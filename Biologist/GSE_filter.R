#############################################################
# Biologist: Elysha Sameth                                  #
#                                                           #
# Filter each MOA for each gene expression technology to    #
# obtain differentialy expressed genes for GSE analysis    #
#############################################################

# Read in the microarray (limma) files
MA_AHR <- read.csv('/projectnb/bf528/users/van-gogh/project_3/Biologist/limma_leflunomide.csv')
MA_DNA <- read.csv('/projectnb/bf528/users/van-gogh/project_3/Biologist/limma_ifosfamide.csv')
MA_PXR <- read.csv('/projectnb/bf528/users/van-gogh/project_3/Biologist/limma_fluconazole.csv')

files <- list(MA_AHR=MA_AHR, MA_DNA=MA_DNA, MA_PXR=MA_PXR)
names <- names(files)

# Function to filter DEG and write probe sets IDs to a file
MA_write_file <- function(df, name) {
  # Filter data by p < 0.05 and abs(FC) > 1.5
  df <- subset(df, P.Value < 0.05 & abs(logFC) > log2(1.5))
  # Get the genes and write it to a file
  # Do not include row name, col name, or quotations (just a list of ids)
  write.table(df$X, paste0('/projectnb/bf528/users/van-gogh/project_3/Biologist/', name, ".txt"), 
              row.names=FALSE, col.names=FALSE, quote=FALSE)
}

# Apply function for each file
lapply(names(files), function(x) MA_write_file(files[[x]], x))


# Read in the RNA-seq (DESeq2) files
RNA_AHR <- read.csv('/projectnb2/bf528/users/van-gogh/project_3/Programmer/deseqResults/deseq_results_AhR.csv')
RNA_DNA <- read.csv('/projectnb2/bf528/users/van-gogh/project_3/Programmer/deseqResults/deseq_results_DNA_Damage.csv')
RNA_PXR <- read.csv('/projectnb2/bf528/users/van-gogh/project_3/Programmer/deseqResults/deseq_results_PXR.csv')

files <- list(RNA_AHR=RNA_AHR, RNA_DNA=RNA_DNA, RNA_PXR=RNA_PXR)
names <- names(files)
# Function to filter DEG and write probe sets/gene names to a file
RNA_write_file <- function(df, name) {
  # Filter data by p < 0.05 and abs(FC) > 1.5
  df <- subset(df, pvalue < 0.05 & abs(log2FoldChange) > log2(1.5))
  # Get the genes and write it to a file
  # Do not include row name, col name, or quotations (just a list of ids)
  write.table(df$X, paste0('/projectnb/bf528/users/van-gogh/project_3/Biologist/', name, ".txt"), 
              row.names=FALSE, col.names=FALSE, quote=FALSE)
}

# Apply function for each file
lapply(names(files), function(x) RNA_write_file(files[[x]], x))