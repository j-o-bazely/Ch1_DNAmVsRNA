##########           JB0.1A - Object_Preparation              ##########

###


# This script takes in all the samples to be used in this analysis
# and converts them into objects suitable for analysis with the MethylKit package

# Note: Requires minimum 128Gb to run


###        ###       ###        ###          ###         ###        ###  
#####               


#### 1. ENV PREP                                                                ####

## PACKAGES
library('methylKit')
library('ggplot2')

## PATHS
WorkingDir <- '/data/SBCS-EizaguirreLab/James_B/phd/JB1_DNAm_RNA_comparison/JB1.M_METHYLATION'

objPath <- paste(WorkingDir, '/JB1.M_OBJECTS', sep='')

objOutPath <-paste(objPath, '/JB1.M_MethylKitObjects', sep='')

figOutPath <- paste(WorkingDir, '/JB1.M_FIGURES', sep='')

MetadataPath <- paste(objPath, '/JB1.M_Metadata', sep='')

## custom functions
source(file.path("/data/SBCS-EizaguirreLab/James_B/JB_Island_project/JB_Code/customRfunctions_CY.R"))



#### 2. OBJECT GENERATION
#### Metadata PREP ####
# Load metadata
preMeta <- read.csv(file.path(MetadataPath, 'JB1_UncleanMetadata.csv'))
preMeta <- subset(preMeta, select = c('Female_ID', 'Treatment', 'Nest_ID'))

# Load list of samples
sample_list <- read_lines(file.path(MetadataPath, 'JB1_SampleList.txt'))

# Create empty dataframe for metadata
metadata <- data.frame()

# Loop through preMeta df
for(i in 1:length(preMeta$Nest_ID)){
  # For each row in preMeta, loop through the sample list to see 
  # which hatchling it is for, based on Nest_ID
  for(j in sample_list)
    if(grepl(preMeta$Nest_ID[i], j)){
      # If it's the correct hatchling, add that information to the metadata
      metadata <- rbind(metadata, preMeta[i,])
    }
}

metadata <- metadata %>% 
  rename(Depth = Treatment)

# Add the hatchling IDs to the metadata, leaving us with the correct metadata 
#for each hatchling
metadata$hatchling_ID <- sample_list

# Save the cleaned metadata

write.csv(metadata, file.path(MetadataPath, '/JB1_CleanMetadata.csv'))


#### Global Methylation Prep                                                    ####

## List of paths to sample files
prepMethObj <- function(trt = c('Depth', 'Female_ID')){
  
  # Start timer
  start.time <- Sys.time()
  
  
  # Check user input matches one of two options for trt
  trt <- match.arg(trt)
  
  # Path to sample methylation data
  SamplePathVector <- c()
  for(i in 1:length(metadata$hatchling_ID)){
    SamplePathVector <- append(SamplePathVector, paste(objPath, sprintf('JB1.M_RawReads/JB1_%s.CpG_merged.cov.gz', metadata$hatchling_ID[i]), sep = '/'))
  }
  
  # List of paths to sample methylation datum
  SamplePathList <- list()
  
  for (i in 1:length(SamplePathVector)) {
    SamplePathList[i] <- SamplePathVector[i]
  }
  
  ## Vector of treatment for each sample (deep / shallow)
  if(trt == 'Depth')  {
    treatmentVector <- as.numeric(factor(metadata$Depth))
  }
  
  ## Vector of treatment for each sample (mother_ID)
  if(trt == 'Female_ID')  {
    treatmentVector <- as.numeric(factor(metadata$Female_ID))
  }
  
  ## List of sample IDs
  sampleIdList <- list()
  
  for (i in 1:length(metadata$hatchling_ID)){
    sampleIdList[i] <- metadata$hatchling_ID[i]
  }
  
  # Running MethRead()
  methObj <- methRead(location = SamplePathList,
                      sample.id = sampleIdList,
                      treatment = treatmentVector,
                      context = 'CpG',
                      assembly = 'CarCar_QM_v1_2021_12_Scaff_With_Chr0',
                      pipeline="bismarkCoverage"
  )
  
  ## Coverage & Global Stats ##
  getMethylationStats(methObj[[1]], plot=TRUE)
  
  getCoverageStats(methObj[[1]], plot=TRUE)
  
  
  #### Data Cleaning & Object Merging                                           ####
  
  # Filter object by coverage (low coverage filter and >99.9th percentile in each sample)
  # Filter by 5X
  MethFilter5Cov <-filterByCoverage(methObj,
                                    lo.count=5,
                                    lo.perc=NULL,
                                    hi.count=NULL,
                                    hi.perc=99.9)
  
  MethNormFilter5Cov=normalizeCoverage(MethFilter5Cov)
  
  # Merge samples together into a single united coverage object
  uniteMeth100pc <- methylKit::unite(MethNormFilter5Cov, destrand = FALSE)
  
  
  # Save object
  saveRDS(uniteMeth100pc, file.path(paste('/data/SBCS-EizaguirreLab/James_B/phd/JB1_DNAm_RNA_comparison/JB1.M_METHYLATION/JB1.M_OBJECTS/JB1.M_MethylKitObjects/JB1_UniteCovObj100pc_', trt, '_trt.RDS', sep = '')))
  
  # Remove the uniteMeth object to save memory
  rm(uniteMeth100pc)
  gc()
  
  #### 75% Coverage                                                             ####
  
  # find smallest population to divide by for 75% coverage
  result <- metadata %>%  group_by(get(trt)) %>% summarise(count = n()) %>% arrange(desc(count))
  min <- as.integer(round((as.integer(result[which.min(result$count),][2])/4)*3))
  
  # Run unite() again, with the minimum number of samples per group required to have coverage at a base reduced to 75% for inclusion
  uniteMeth_75pc = methylKit::unite(MethNormFilter5Cov,
                                    min.per.group= (min))
  
  # Convert to a methylBase object
  uniteMeth_75pc = as(uniteMeth_75pc,"methylBase")
  
  # Save object
  saveRDS(uniteMeth_75pc, file.path(paste('/data/SBCS-EizaguirreLab/James_B/phd/JB1_DNAm_RNA_comparison/JB1.M_METHYLATION/JB1.M_OBJECTS/JB1.M_MethylKitObjects/JB1_UniteCovObj75pc_', trt, '_trt.RDS', sep = '')))
  
  rm(uniteMeth_75pc)
  gc()
  # End timer 
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  time.taken
}

# Takes 10 minutes
prepMethObj(trt = 'Depth')

# Takes 9 minutes
prepMethObj(trt = 'Female_ID')


#### Differential Methylation ####

## Depth : 100% coverage ##
# Calculate differential methylation between deep and shallow hatchlings #
diffMethRuntimeFunc<- function(methObj100pc_D){
  start.time <- Sys.time()
  diffMethObj100pc_D <- calculateDiffMeth(methObj100pc_D,
                                     overdispersion = "MN",
                                     adjust="BH")
  saveRDS(diffMethObj100pc_D, file = file.path('/data/SBCS-EizaguirreLab/James_B/phd/JB1_DNAm_RNA_comparison/JB1.M_METHYLATION/JB1.M_OBJECTS/JB1.M_MethylKitObjects/JB1.M_DiffMethObjects/JB1_DiffMeth_UniteCovObj100pc_Depth_trt.RDS'))
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  time.taken
}

diffMethRuntimeFunc(methObj100pc_D)


## Save differential methylation object
saveRDS(diffMethObj100pc_D, file = file.path('/data/SBCS-EizaguirreLab/James_B/phd/JB1_DNAm_RNA_comparison/JB1.M_METHYLATION/JB1.M_OBJECTS/JB1.M_MethylKitObjects/JB1.M_DiffMethObjects/JB1_DiffMeth_UniteCovObj100pc_Depth_trt.RDS'))


## Depth : 75% coverage ##
# Calculate differential methylation between deep and shallow hatchlings #

## Load 75% diff meth object
unitecov75pc_D <- readRDS(file.path('/data/SBCS-EizaguirreLab/James_B/phd/JB1_DNAm_RNA_comparison/JB1.M_METHYLATION/JB1.M_OBJECTS/JB1.M_MethylKitObjects/JB1_UniteCovObj75pc_Depth_trt.RDS'))

### how many sites are there in the whole genome?
n_uniteCov_75pc_D = nrow(unitecov75pc_D)

### main diff meth calculation 
### the automatic p value correction is SLIM and is useful to have so I don't change it here
diffMethObj75pc_D <- calculateDiffMeth(unitecov75pc_D,
                                  overdispersion = "MN",  ##### Note we should include this correction as it looks like the DNAm data are "overdispersed"
                                  chunk.size = 10000)

### this uses the BH correction based on the number of sites in the chromosome
diffMethObj75pc_D$p_adj_BH_n_Chrom = p.adjust(diffMethObj75pc_D$pvalue, method = "BH",
                                         n = n_uniteCov_75pc_D)

### this uses the BH correction based on the number of sites in the whole genome
diffMethObj75pc_D$p_adj_BH_n_WG = p.adjust(diffMethObj75pc_D$pvalue, method = "BH",
                                      n = n_uniteCov_75pc_D )

# Save differential methylation object
saveRDS(diffMethObj75pc_D, file = file.path('/data/SBCS-EizaguirreLab/James_B/phd/JB1_DNAm_RNA_comparison/JB1.M_METHYLATION/JB1.M_OBJECTS/JB1.M_MethylKitObjects/JB1.M_DiffMethObjects/JB1_DiffMeth_UniteCovObj75pc_Depth_trt.RDS'))

