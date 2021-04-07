### featureCount R Script
### Andrew Gjelsteen BF528 (Spring 2021) - van-gogh group

library(ggplot2)

## Processing Each Count file to combine into a single csv ##
### SRR1177981
SRR1177981 <- read.table("SRR1177981_counts", sep='\t', header = TRUE)
SRR1177981 <- subset(SRR1177981, select = c(1, 7))
colnames(SRR1177981)[2] <- "SRR1177981" #rename column for future purposes
### SRR1177982
SRR1177982 <- read.table("SRR1177982_counts", sep='\t', header = TRUE)
SRR1177982 <- subset(SRR1177982, select = c(1, 7))
colnames(SRR1177982)[2] <- "SRR1177982" #rename column for future purposes
### SRR1177983
SRR1177983 <- read.table("SRR1177983_counts", sep='\t', header = TRUE)
SRR1177983 <- subset(SRR1177983, select = c(1, 7))
colnames(SRR1177983)[2] <- "SRR1177983" #rename column for future purposes
### SRR1178008
SRR1178008 <- read.table("SRR1178008_counts", sep='\t', header = TRUE)
SRR1178008 <- subset(SRR1178008, select = c(1, 7))
colnames(SRR1178008)[2] <- "SRR1178008" #rename column for future purposes
### SRR1178009
SRR1178009 <- read.table("SRR1178009_counts", sep='\t', header = TRUE)
SRR1178009 <- subset(SRR1178009, select = c(1, 7))
colnames(SRR1178009)[2] <- "SRR1178009" #rename column for future purposes
### SRR1178010
SRR1178010 <- read.table("SRR1178010_counts", sep='\t', header = TRUE)
SRR1178010 <- subset(SRR1178010, select = c(1, 7))
colnames(SRR1178010)[2] <- "SRR1178010" #rename column for future purposes
### SRR1178014
SRR1178014 <- read.table("SRR1178014_counts", sep='\t', header = TRUE)
SRR1178014 <- subset(SRR1178014, select = c(1, 7))
colnames(SRR1178014)[2] <- "SRR1178014" #rename column for future purposes
### SRR1178021
SRR1178021 <- read.table("SRR1178021_counts", sep='\t', header = TRUE)
SRR1178021 <- subset(SRR1178021, select = c(1, 7))
colnames(SRR1178021)[2] <- "SRR1178021" #rename column for future purposes
### SRR1178047
SRR1178047 <- read.table("SRR1178047_counts", sep='\t', header = TRUE)
SRR1178047 <- subset(SRR1178047, select = c(1, 7))
colnames(SRR1178047)[2] <- "SRR1178047" #rename column for future purposes

## Merging the files into a single dataframe. Merging by column 'Geneid'
# merge() operates on only two dataframes at a time, so I'll use nested functions

featureMerged <- merge(SRR1177981, 
                       merge(SRR1177982, 
                             merge(SRR1177983, 
                                   merge(SRR1178008, 
                                         merge(SRR1178009, 
                                               merge(SRR1178010, 
                                                     merge(SRR1178014,
                                                           merge(SRR1178021, SRR1178047, by = 'Geneid'),
                                                           by = 'Geneid'),
                                                     by = 'Geneid'),
                                               by = 'Geneid'),
                                         by = 'Geneid'),
                                   by = 'Geneid'),
                             by = 'Geneid'),
                       by = 'Geneid')

write.csv(featureMerged, 'features_Total.csv')

png("boxPlotMultiPanel.png") #Creates .png file to export this as one figure
## Visualizing data with boxplots:

par(mfrow = c(3,3))

# # filter out rows that have any zeros
# cnts <- subset(cnts,rowSums(cnts==0)==0)
# # change Geneid from column into rownames
# cnts_named <- cnts[,-1]
# rownames(cnts_named) <- cnts[,1]
# cnts <- cnts_named

boxplot(log10(SRR1177981$SRR1177981), main = 'SRR1177981 distribution', ylab = 'log(counts)', cex.lab = 1.1)
boxplot(log10(SRR1177982$SRR1177982), main = 'SRR1177982 distribution', ylab = 'log(counts)', cex.lab = 1.1)
boxplot(log10(SRR1177983$SRR1177983), main = 'SRR1177983 distribution', ylab = 'log(counts)', cex.lab = 1.1)
boxplot(log10(SRR1178008$SRR1178008), main = 'SRR1178008 distribution', ylab = 'log(counts)', cex.lab = 1.1)
boxplot(log10(SRR1178009$SRR1178009), main = 'SRR1178009 distribution', ylab = 'log(counts)', cex.lab = 1.1)
boxplot(log10(SRR1178010$SRR1178010), main = 'SRR1178010 distribution', ylab = 'log(counts)', cex.lab = 1.1)
boxplot(log10(SRR1178014$SRR1178014), main = 'SRR1178014 distribution', ylab = 'log(counts)', cex.lab = 1.1)
boxplot(log10(SRR1178021$SRR1178021), main = 'SRR1178021 distribution', ylab = 'log(counts)', cex.lab = 1.1)
boxplot(log10(SRR1178047$SRR1178047), main = 'SRR1178047 distribution', ylab = 'log(counts)', cex.lab = 1.1)

dev.off()

rn4 <- read.csv('group_3_rna_info.csv')
controlCounts <- read.csv('control_counts.csv')
controlCounts <- subset(controlCounts, select = c("Geneid", "SRR1178035", "SRR1178074", "SRR1178075")) # subsetting out the control samples

DESeqData <- merge(featureMerged, controlCounts, by = "Geneid")

### Installing bioconductor locally to access
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")


BiocManager::install("apeglm")
#########

browseVignettes("DESeq2")

library(DESeq2)

# load counts
cnts <- DESeqData

# filter out rows that have any zeros
cnts <- subset(cnts,rowSums(cnts==0)==0)
# change Geneid from column into rownames
cnts_named <- cnts[,-1]
rownames(cnts_named) <- cnts[,1]
cnts <- cnts_named


### Subsetting three count groups:
cnts1 <- subset(cnts, select = c("SRR1178035", "SRR1178074","SRR1178075", "SRR1177981", "SRR1177982", "SRR1177983"))
cnts2 <- subset(cnts, select = c("SRR1178035", "SRR1178074","SRR1178075", "SRR1178008", "SRR1178009", "SRR1178010"))
cnts3 <- subset(cnts, select = c("SRR1178035", "SRR1178074","SRR1178075", "SRR1178014", "SRR1178021", "SRR1178047"))

### Setting Geneid as the row names
info2 <- info[,-1]
rownames(info2) <- info[,1]
info <- info2


# sample information
info <- read.csv("group_EX_rna_info.csv")
info2 <- info[,-1]
rownames(info2) <- info[,1]
info <- info2


##### create the DESeq object (First) #####
dds1 <- DESeqDataSetFromMatrix(
  countData = cnts1,
  colData = info,
  design= ~ mode_of_action
)

# relevel mode_of_action as factor
dds1$mode_of_action <- relevel(dds1$mode_of_action, ref='Control')

# run DESeq
dds1 <- DESeq(dds1)
res1 <- results(dds1, contrast=c('mode_of_action','DNA_Damage','Control'))
res1 <- lfcShrink(dds1, coef=2)

# write out DE results
write.csv(res1,'deseq_results_DNA_Damage.csv')

# write out matrix of normalized counts
write.csv(counts(dds1,normalized=TRUE),'deseq_norm_counts_DNA_Damage.csv')

##### create the DESeq object (Second) #####
dds2 <- DESeqDataSetFromMatrix(
  countData = cnts2,
  colData = info,
  design= ~ mode_of_action
)

# relevel mode_of_action as factor
dds2$mode_of_action <- relevel(dds2$mode_of_action, ref='Control')

# run DESeq
dds2 <- DESeq(dds2)
res2 <- results(dds2, contrast=c('mode_of_action','AhR','Control'))
res2 <- lfcShrink(dds2, coef=2)

# write out DE results
write.csv(res2,'deseq_results_AhR.csv')

# write out matrix of normalized counts
write.csv(counts(dds2,normalized=TRUE),'deseq_norm_counts_AhR.csv')
 
##### create the DESeq object (Third) #####
dds3 <- DESeqDataSetFromMatrix(
  countData = cnts3,
  colData = info,
  design= ~ mode_of_action
)

# relevel mode_of_action as factor
dds3$mode_of_action <- relevel(dds3$mode_of_action, ref='Control')

# run DESeq
dds3 <- DESeq(dds3)
res3 <- results(dds3, contrast=c('mode_of_action','CAR/PXR','Control'))
res3 <- lfcShrink(dds3, coef=2)

# write out DE results
write.csv(res3,'deseq_results_PXR.csv')

# write out matrix of normalized counts
write.csv(counts(dds3,normalized=TRUE),'deseq_norm_counts_PXR.csv')

##### Import the csv files and filter by p-value #####

deseq1 <- read.csv("deseq_results_DNA_damage.csv", na.strings=FALSE)
deseq1 <- deseq1[deseq1$padj < 0.05, ]
deseq1 <- na.omit(deseq1)

length(deseq1$log2FoldChange[deseq1$pvalue < 0.05])

deseq2 <- read.csv("deseq_results_AhR.csv", na.strings=FALSE)
deseq2 <- deseq2[deseq2$padj < 0.05, ]
deseq2 <- na.omit(deseq2)
length(deseq2$padj)
length(deseq2$log2FoldChange[deseq1$pvalue < 0.05])

deseq3 <- read.csv("deseq_results_PXR.csv", na.strings=FALSE)
deseq3 <- deseq3[deseq3$padj < 0.05, ]
deseq3 <- na.omit(deseq3)
length(deseq3$padj)
length(deseq3$log2FoldChange[deseq3$pvalue < 0.05])
######################

## Create histograms of fold change values from the significant DE genes for each analysis.
## Also create scatter plots of fold change vs nominal p-value.

png("histogram.png")
par(mfrow = c(1, 3))
hist(deseq1$log2FoldChange[deseq1$pvalue < 0.05], main = "DE genes - Ifosfamide", xlab = 'Log Fold Change', ylab = 'count', col = 'steelblue2')
hist(deseq2$log2FoldChange[deseq2$pvalue < 0.05], main = "DE genes - Leflunomide", xlab = 'Log Fold Change', ylab = 'count', col = 'steelblue2')
hist(deseq3$log2FoldChange[deseq3$pvalue < 0.05], main = "DE genes - Fluconazole", xlab = 'Log Fold Change', ylab = 'count', col = 'steelblue2')
dev.off()

### Scatterplots for these 
library(ggplot2)
install.packages(c("hrbrthemes"))
library(hrbrthemes)
library(devtools)

## Getting -log(pvalue) column for each of the 3 dataframes.
deseq1$logpValue <- -log10(as.numeric(deseq1[,5]))
deseq2$logpValue <- -log10(as.numeric(deseq2[,5]))
deseq3$logpValue <- -log10(as.numeric(deseq3[,5]))

deseq1$color = "black"
deseq1$color[deseq1$logpValue >= 1.3] = "steelblue3"
deseq1$color[abs(deseq1$log2FoldChange) > 1.5] = "darkgreen"
deseq1$color[abs(deseq1$log2FoldChange) > 1.5 & deseq1$logpValue > 1.3] = 'red'

deseq2$color = "black"
deseq2$color[deseq2$logpValue > 1.3] = "steelblue3"
deseq2$color[deseq2$log2FoldChange > 1.5] = "darkgreen"
deseq2$color[deseq2$log2FoldChange < -1.5] = "darkgreen"
deseq2$color[deseq2$log2FoldChange > 1.5 & deseq2$logpValue > 1.3] = 'red'
deseq2$color[deseq2$log2FoldChange < -1.5 & deseq2$logpValue > 1.3] = 'red'

deseq3$color = "black"
deseq3$color[deseq3$logpValue > 1.3] = "steelblue3"
deseq3$color[deseq3$log2FoldChange > 1.5 | deseq3$log2FoldChange < -1.5] = "darkgreen"
deseq3$color[deseq3$log2FoldChange > 1.5 && deseq3$logpValue > 1.3] = 'red'
deseq3$color[deseq3$log2FoldChange < -1.5 && deseq3$logpValue > 1.3] = 'red'

png('Scatterplots.png')

par(mfrow = c(2, 2), cex.main = 0.85, cex.lab = 0.85)
## First Plot
plot(logpValue ~ log2FoldChange, data = deseq1, xlab = expression('log'[2]*'fold change'), 
  ylab = expression('-log'[10]*'(p)'), main = 'Fold Change vs Nominal p-value - Fluconazole',
    pch = 20, col = deseq1$color, cex = 0.95) 

abline(h = 1.3, lty = 2)
abline(v = c(-1.5,1.5), lty = c(2,2))
###Second plot
plot(logpValue ~ log2FoldChange, data = deseq2, xlab = expression('log'[2]*'fold change'), 
     ylab = expression('-log'[10]*'(p)'), main = 'Fold Change vs Nominal p-value - Leflunomide',
     pch = 20, col = deseq2$color, cex = 0.95)
abline(h = 1.3, lty = 2)
abline(v = c(-1.5,1.5), lty = c(2,2))

### Third Plot
plot(logpValue ~ log2FoldChange, data = deseq3, xlab = expression('log'[2]*'fold change'), 
     ylab = expression('-log'[10]*'(p)'), main = 'Fold Change vs Nominal p-value - Ifosfamide',
     pch = 20, col = deseq3$color, cex = 0.95)
abline(h = 1.3, lty = 2)
abline(v = c(-1.5,1.5), lty = c(2,2))

# legend("bottom", xpd = TRUE, horiz = TRUE, inset = c(0, 0), bty = "n", c("NS",expression('-log'[2]*'FC'),'p-value', expression('-log'[2]*'FC and p-value')), 
#        pch=c(20, 20, 20,20),col=c("black","steelblue3","darkgreen","red"))

dev.off()

