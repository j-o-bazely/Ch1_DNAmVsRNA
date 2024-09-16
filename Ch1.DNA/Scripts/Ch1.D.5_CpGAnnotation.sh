####        Ch1.D.5_CpGAnnotation.R        ####

###        ###       ###        ###          ###         ###        ###   

# This script annotates all identified methylated sites
# with their location in the genome and their associated gene 


###        ###       ###        ###          ###         ###        ###   

####    ENV PREP                                                                ####

## PACKAGES
library('methylKit')
library('ggplot2')
library('GenomicRanges')



## PATHS
WorkingDir <- '/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_DNAmVsRNA/Ch1.DNA'

dataPath <- '/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_dataStorage'

DataOutDir <- '/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_dataStorage/Ch1.M_data/Ch1.M.5_data'

figOutPath <- '/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_DNAmVsRNA/Ch1.DNA/Figures'

MetadataPath <- '/data/SBCS-EizaguirreLab/James_B/cleanPHD/Ch1_DNAmVsRNA/Ch1.DNA/Metadata'

## OBJECTS 
gff3Annot <- rtracklayer::readGFF(file.path(dataPath, 'Ch1_inputData/Ch1_ChangGenome/genomic.gtf'))
gff3Annot_old <- rtracklayer::readGFF(file.path('/data/SBCS-EizaguirreLab/James_B/phd/JB1_DNAm_RNA_comparison/JB1.R_RNA_UPDATE/JB1.R_data/Chang23CanadaGenome/Chang23LoggerheadAnnotation.gtf'))
bedAnnot <- readTranscriptFeatures(file.path(dataPath, 'Ch1_inputData/Ch1_ChangGenome/Chang23LoggerheadAnnotationTEMP.bed12'),
                                   remove.unusual = FALSE,
                                   up.flank = 1500,
                                   down.flank = 500,
                                   unique.prom=FALSE)

## Load object to be annotated
OBJ <- readRDS(file.path(dataPath, 'Ch1.M_data/Ch1.M.1_data/ch1.M.1_uniteMeth75pc.RDS'))

####      OBJECT GENERATION                                                     ####

# If object is of class 'methylBase', convert to a percent methylation matrix
if(class(OBJ) == 'methylBase'){
  OBJ_matrix <- percMethylation(OBJ)
}

# Create initial GRange object
#     GRange object containing location of every CpG in the object      #


# Reassign column names, avoids clashes with genomation later on
OBJ$original_pos <-  OBJ$start
OBJ$original_chrom <- OBJ$chr
colnames(OBJ)

# Vector containing chromosome, start and end location of CpG
myPos = paste(OBJ$chr, OBJ$end, OBJ$original_pos, OBJ$original_chrom)

# Change this vector into a dataframe 
df <- data.frame(chr=sapply(strsplit(myPos, " "), `[`, 1),
                 start=sapply(strsplit(myPos, " "), `[`, 2),
                 end=sapply(strsplit(myPos, " "), `[`, 2),
                 original_pos=sapply(strsplit(myPos, " "), `[`, 3),
                 original_chrom=sapply(strsplit(myPos, " "), `[`, 4))

# Convert this to a GRangesOBJ of every CpG position in the original object
GRangeOBJ <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)


# Use annotateWithGeneParts from genomation, to annotate each of these CpGs
# with features such as distance from a gene and the name of that gene 
A = annotateWithGeneParts(target = as(GRangeOBJ,"GRanges"), feature = bedAnnot)

# Check that the names are promoters, exons, introns and intergenic in A@members before this
GRangeOBJ$featureType = ifelse(A@members[,1]==1, "promoters",
                               ifelse(A@members[,2]==1, "exons",
                                      ifelse(A@members[,3]==1, "introns", "intergenic")))



####    Gene Name Addition                                                      ####


##      Handling genic CpG      ##
# Only genic CpG need to be assigned to genes,
# Subset the GRanges object to remove any 'intergenic' sites
GRangeOBJ1=GRangeOBJ[!GRangeOBJ$featureType %in% "intergenic"]

# Prepare a column in this subset GRanges object to take gene info
GRangeOBJ1$geneInfo <- NA

# Function: Add gene info from BED12 file to the GRanges object
# x = Feature type (e.g. 'exons')
add_geneInfo_genic <- function(x, GRangeOBJ, bedAnnot){  
  
  # findOverlaps() is assigned to ov.
  # This sub-function finds any overlaps between the BED12 file for the given feature type,
  # and sites in the GRanges object for the given feature type.
  ov = GenomicRanges::findOverlaps(
    bedAnnot[[x]],
    GRangeOBJ[GRangeOBJ$featureType %in% x,])
  
  # The 'name' column in the BED12 file, containing the names of the genes,
  # is merged with the 'gene info' column in the correct position in the GRanges object
  mcols(GRangeOBJ[GRangeOBJ$featureType %in% x,])[subjectHits(ov), "geneInfo"] = mcols(bedAnnot[[x]])[queryHits(ov), "name"]
  
  # The GRanges object containing the gene name is now returned
  return(GRangeOBJ)
}

# Run this function on exon, intron, and promoter region CpGs
myGRangeOBJ1 = add_geneInfo_genic("exons", GRangeOBJ1, bedAnnot) # Add gene names for sites on exons
myGRangeOBJ2 = add_geneInfo_genic("introns", myGRangeOBJ1, bedAnnot) # Add gene names for sites on introns
myGRangeOBJ_genic = add_geneInfo_genic("promoters", myGRangeOBJ2, bedAnnot) # Add gene names for sites on promoters


##    Handling intergenic CpG     ##                                                      ##
# Subset for only intergenic CpG 
interGRangeOBJ1 =GRangeOBJ[GRangeOBJ$featureType %in% "intergenic"]

# Rerun annotateWithGeneParts
interA = annotateWithGeneParts(target = as(interGRangeOBJ1,"GRanges"), feature = bedAnnot)

# rows2rm is an index of positions that are intergenic AND 10,000+ bp from a TSS,
# so therefore are considered unrelated to genes
rows2rm = which((interA@dist.to.TSS$dist.to.feature>10000 | interA@dist.to.TSS$dist.to.feature< -10000) &rowSums(interA@members) %in% 0)

if (is_empty(rows2rm)){interGRangeOBJ1 = interGRangeOBJ1} else { interGRangeOBJ1 = interGRangeOBJ1[-rows2rm,] }

# Time to re-annotate this new Granges object minus the truly intergenic regions
# using annotateWithGeneParts!
interB = annotateWithGeneParts(as(interGRangeOBJ1,"GRanges"), bedAnnot)

# We then want to know the genes these CpG are within 10 kb of 
interC = getAssociationWithTSS(interB)

# Now add these associated gene IDs to the initial GRange object subset for intergenic CpG
# So now, every intergenic CpG that wasn't too far away, has an associated gene
interGRangeOBJ1$geneInfo = interC$feature.name

# Shorten the name just for ease of use
my_GRangeOBJ_inter <- interGRangeOBJ1

# merge the intergenic CpG back with the genic CpG for a full cohort
# of CpG associated with genes
finalGRangeOBJ =c(myGRangeOBJ_genic, my_GRangeOBJ_inter)


###     OBJ Subsetting      ###

# Use finalGRangeOBJ to subset OBJ so that only gene-associated CpG remain
# Remove the redundant columns from the GRanges object
finalGRangeOBJ$original_pos <- NULL
finalGRangeOBJ$original_chrom <- NULL

# Use selectByOverlap() to subset the unitecov object to only include CpG related to genes
geneAssocOBJ <- selectByOverlap(OBJ, finalGRangeOBJ)


###   Output Object Generation      ###

# Create matrix from gene-associated CpG OBJ
gaMatrix <- as.data.frame(percMethylation(geneAssocOBJ))

# Add gene name column
gaMatrix$gene <- str_extract(finalGRangeOBJ$geneInfo, "(?<=gene-)[^;]+")
  
# Add feature type column
gaMatrix$region <- finalGRangeOBJ$featureType

## Save matrix  ##
saveRDS(gaMatrix, file.path(DataOutDir, 'Ch1.M.5_percMeth75pc_Annot_Matrix.RDS'))


