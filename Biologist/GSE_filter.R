#############################################################
# Biologist: Elysha Sameth                                  #
#                                                           #
# Filter each MOA for each gene expression technology to    #
# obtain differentially expressed genes for GSE analysis    #
#############################################################

# Read in the microarray (limma) files
MA_AHR <- read.csv('/projectnb/bf528/users/van-gogh/project_3/Biologist/limma_leflunomide.csv')
MA_DNA <- read.csv('/projectnb/bf528/users/van-gogh/project_3/Biologist/limma_ifosfamide.csv')
MA_PXR <- read.csv('/projectnb/bf528/users/van-gogh/project_3/Biologist/limma_fluconazole.csv')

files <- list(MA_AHR=MA_AHR, MA_DNA=MA_DNA, MA_PXR=MA_PXR)
names <- names(files)

write_file <- function(df, name) {
  # Filter data by adj. p-value
  df <- subset(df, adj.P.Val < 0.05)
  # Get the genes and write it to a file
  write.table(df$X, paste0('/projectnb/bf528/users/van-gogh/project_3/Biologist/', name, ".txt"), 
              row.names=FALSE, col.names=FALSE)
}

# Apply function for each file
lapply(names(files), function(x) write_file(files[[x]], x))