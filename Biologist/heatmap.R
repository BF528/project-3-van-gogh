#############################################################
# Biologist: Elysha Sameth                                  #
#                                                           #
# Create a clustered heatmap of the normalized counts       #
# DESeq2 (RNA-Seq) to assess MOA grouping                   #
#############################################################

library(pheatmap)
library(grid)

# Read in the files
AHR_count <- read.csv('/projectnb2/bf528/users/van-gogh/project_3/Programmer/deseqResults/deseq_norm_counts_AhR.csv', row.names="X")
DNA_count <- read.csv('/projectnb2/bf528/users/van-gogh/project_3/Programmer/deseqResults/deseq_norm_counts_DNA_Damage.csv', row.names="X")
PXR_count <- read.csv('/projectnb2/bf528/users/van-gogh/project_3/Programmer/deseqResults/deseq_norm_counts_PXR.csv', row.names="X")

# Merge all treated samples and controls
merged <- cbind(AHR_count, PXR_count, DNA_count)
# Get the mean of common controls
final_merged <- as.data.frame(
  sapply(unique(names(merged)),
         function(col) rowMeans(merged[names(merged) == col])))

# Create metadata for ColSideColors  
metadata <- data.frame(sample=c("AhR", "AhR", "AhR", "Control", "Control", "Control", 
                                "CAR.PXR", "CAR.PXR", "CAR.PXR", 
                                "DNA", "DNA", "DNA", "Control", "Control", "Control")
                       )
rownames(metadata) <- colnames(final_merged)

# Filter genes that have a variance significantly different from the median variance of all genes using a threshold of p < 0.01
# Perform a two-tailed alternative
df<-ncol(final_merged)-1
chi_sq_lower<-qchisq(0.01/2, df)
chi_sq_upper<-qchisq(1-0.01/2, df)

# Compute the test statistic for each gene: T = (N-1)(s/o)**2 where df=N-1, variance=s, median variance=o
final_merged$variance <- apply(final_merged, 1, var)
final_merged$test_stat <- df*final_merged$variance/median(final_merged$variance)
chi_filter<-subset(final_merged, test_stat > chi_sq_upper | test_stat < chi_sq_lower)
# Remove extra columns
chi_filter<-subset(chi_filter, select = -c(variance, test_stat))

# Have a coefficient of variation > 0.406
variation_filter<-subset(chi_filter, apply(chi_filter, 1, function(x) sd(x)/mean(x)) > 0.406)

# Create the heatmap
setHook('grid.newpage', function() pushViewport(viewport(x=1,y=1,width=0.9,height=0.9,
                                                         name='vp', just=c('right', 'top'))),
        action='prepend')
pheatmap(as.matrix(merged), annotation_col = metadata, scale='row',
         show_rownames=FALSE, main='RNA-Seq Normalized Count Across Samples',
         cluster_rows=FALSE, annotation_names_col = FALSE)
setHook('grid.newpage', NULL, 'replace')
grid.text('Sample', y=-0.07, gp=gpar(fontsize=12))
grid.text('Gene', x=-0.07, rot=90, gp=gpar(fontsize=12))


################################################
# Extra analysis to assess incorrect clustering
library(PCAtools)

# Perform a pca
p <- pca(final_merged, metadata=metadata)
# Biplot(shows PC1 and PC2)
biplot(p, colby='sample',
       colLegendTitle='Sample Type', 
       legendPosition='right',
       title='PCA of RNA-Seq Normalized Counts')
# Pairsplot (first 5 PCs)
pairsplot(p, colby='sample', 
          colkey=c(AhR='chartreuse3', CAR.PXR='coral1', Control='cyan2', 'DNA'='darkorchid1'))
