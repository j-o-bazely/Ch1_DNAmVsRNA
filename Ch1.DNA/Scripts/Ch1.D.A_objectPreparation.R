##########           Ch1.DNA_1_ObjectPreparation.r                              ##########

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
WorkingDir <- '/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_DNAmVsRNA/Ch1.DNA'

dataPath <- '/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_dataStorage'

DataOutDir <- '/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_dataStorage/Ch1.M_data/Ch1.M.1_data'

figOutPath <- '/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_DNAmVsRNA/Ch1.DNA/Figures'

MetadataPath <- '/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_DNAmVsRNA/Ch1.DNA/Metadata'

## custom functions
source(file.path("/data/SBCS-EizaguirreLab/James_B/JB_Island_project/JB_Code/customRfunctions_CY.R"))



#### 2. OBJECT GENERATION

## Load metadata
metadata <- read.csv(file = file.path(MetadataPath, 'JB1_CleanMetadata.csv'))

#### Global Methylation Prep                                                    ####

## This function generates a methylRaw object used for global analysis and 
## List of paths to sample files
prepMethObj <- function(trt = c('Depth', 'Female_ID')){
  
  # Start timer
  start.time <- Sys.time()
  
  
  # Check user input matches one of two options for trt
  trt <- match.arg(trt)
  
  # Path to sample methylation data
  SamplePathVector <- c()
  for(i in 1:length(metadata$hatchling_ID)){
    SamplePathVector <- append(SamplePathVector, paste(dataPath, sprintf('Ch1.ReMap_data/Ch1.ReMap.4_data/destrandedCalls/%s/%s.CpG_merged.cov.gz', metadata$hatchling_ID[i], metadata$hatchling_ID[i]), sep = '/'))
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
  methObj <- methRead(location = SamplePathList[1:9],
                      sample.id = sampleIdList[1:9],
                      treatment = treatmentVector[1:9],
                      context = 'CpG',
                      assembly = 'Chang2023',
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
  saveRDS(uniteMeth100pc, file.path(DataOutDir, 'ch1.M.1_uniteMeth100pc.RDS'))
  
  # Remove the uniteMeth object to save memory
  rm(uniteMeth100pc)
  gc()
  
  #### 75% Coverage                                                             ####
  
  # find smallest population to divide by for 75% coverage
  result <- metadata %>%  group_by(get(trt)) %>% summarise(count = n()) %>% arrange(desc(count))
  min <- as.integer(round((as.integer(result[which.min(result$count),][2])/4)*3))
  
  # Run unite() again, with the minimum number of samples per group required to have coverage at a base reduced to 75% for inclusion
  uniteMeth_75pc = methylKit::unite(MethNormFilter5Cov,
                                    min.per.group= 6L)
  
  # Convert to a methylBase object
  uniteMeth_75pc = as(uniteMeth_75pc,"methylBase")
  
  # Save object
  saveRDS(uniteMeth_75pc, file.path(DataOutDir, 'ch1.M.1_uniteMeth75pc.RDS'))
  
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
  diffMethObj100pc_D <- calculateDiffMeth(uniteMeth100pc,
                                          overdispersion = "MN",
                                          adjust="BH",
                                          mc.cores = 5)
  #saveRDS(diffMethObj100pc_D, file = file.path('/data/SBCS-EizaguirreLab/James_B/phd/JB1_DNAm_RNA_comparison/JB1.M_METHYLATION/JB1.M_OBJECTS/JB1.M_MethylKitObjects/JB1.M_DiffMethObjects/JB1_DiffMeth_UniteCovObj100pc_Depth_trt.RDS'))
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  time.taken
}

diffMethRuntimeFunc(uniteMeth100pc)


## Save differential methylation object
saveRDS(diffMethObj100pc_D, file.path(DataOutDir, 'ch1.M.1_diffMeth100pc.RDS'))


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
saveRDS(file.path(DataOutDir, 'ch1.M.1_uniteMeth75pc.RDS'))

