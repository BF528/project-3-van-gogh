# Monil Gandhi - Analyst

#File Description - Run limma on microarray data to find the DE genes and write the csv to the analyst folder

# install packages and load libraries
install.packages("pacman")
require(pacman)  
library(pacman)  
library(reshape2)
library(grid)
library(gridExtra)
install.packages("limma")
library(limma)

#run pacman for downloading the necessary libraries
pacman::p_load(pacman, dplyr, GGally, ggplot2, ggthemes, 
               ggvis, httr, lubridate, plotly, rio, rmarkdown, shiny, 
               stringr, tidyr, ggpubr, grid, gridExtra) 

#import the limma package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")

#import the group 3 toxo info csv
samples <- read.csv("/project/bf528/project_3/groups/group_3_mic_info.csv")

#import the rma table
rma <- read.table('/projectnb/bf528/project_3/samples/liver-normalization-rma.txt',
                  sep='\t',
                  as.is=TRUE,
                  header=TRUE,
                  row.names=1,
)

# subset the full expression matrix to just those in this comparison
rma.subset <- rma[paste0('X',samples$array_id)]

#LEFLUNOMIDE
rma.subset.leflmide <- rma[paste0('X',samples$array_id[samples$chemical=='LEFLUNOMIDE' | (samples$chemical=='Control' & samples$vehicle=='CORN_OIL_100_%')])]

#call design to create a matrix
design <- model.matrix(
  ~factor(
    samples$chemical[samples$vehicle=='CORN_OIL_100_%'],
    levels=c('Control','LEFLUNOMIDE')
  )
)

#assign the column name
colnames(design) <- c('Intercept','LEFLUNOMIDE')

# run limma
fit <- lmFit(rma.subset.leflmide, design)
fit <- eBayes(fit)
t <- topTable(fit, coef=2, n=nrow(rma.subset.leflmide), adjust='BH')

# t_adj_pval <- subset(t, adj.P.Val < 0.05)

#get the DE genes on P.Value < 0.05 & abs(logFC) > log2(1.5)
t_adj_pval <- subset(t, P.Value < 0.05 & abs(logFC) > log2(1.5))

#print the DE genes in LEFLUNOMIDE
print(nrow(t_adj_pval))

#write the csv to the folder
write.csv(t,'/projectnb/bf528/users/van-gogh/project_3/Analyst/limma_leflunomide.csv')



#FLUCONAZOLE
rma.subset.fluzole <- rma[paste0('X',samples$array_id[samples$chemical=='FLUCONAZOLE' | (samples$chemical=='Control' & samples$vehicle=='CORN_OIL_100_%')])]

#call design to create a matrix
design <- model.matrix(
  ~factor(
    samples$chemical[samples$vehicle=='CORN_OIL_100_%'],
    levels=c('Control','FLUCONAZOLE')
  )
)

#assign the column name
colnames(design) <- c('Intercept','FLUCONAZOLE')

# run limma
fit <- lmFit(rma.subset.fluzole, design)
fit <- eBayes(fit)
t <- topTable(fit, coef=2, n=nrow(rma.subset.fluzole), adjust='BH')

# t_adj_pval <- subset(t, adj.P.Val < 0.05)

#get the DE genes on P.Value < 0.05 & abs(logFC) > log2(1.5)
t_adj_pval <- subset(t, P.Value < 0.05 & abs(logFC) > log2(1.5))

#print the DE genes in FLUCONAZOLE
print(nrow(t_adj_pval))

#write the csv to the folder
write.csv(t,'/projectnb/bf528/users/van-gogh/project_3/Analyst/limma_fluconazole.csv')


#IFOSFAMIDE
rma.subset.isomide <- rma[paste0('X',samples$array_id[samples$chemical=='IFOSFAMIDE' | (samples$chemical=='Control' & samples$vehicle=='SALINE_100_%')])]

#call design to create a matrix
design <- model.matrix(
  ~factor(
    samples$chemical[samples$vehicle=='SALINE_100_%'],
    levels=c('Control','IFOSFAMIDE')
  )
)

#assign the column name
colnames(design) <- c('Intercept','IFOSFAMIDE')

# run limma
fit <- lmFit(rma.subset.isomide, design)
fit <- eBayes(fit)
t <- topTable(fit, coef=2, n=nrow(rma.subset.isomide), adjust='BH')

# t_adj_pval <- subset(t, adj.P.Val < 0.05)

#get the DE genes on P.Value < 0.05 & abs(logFC) > log2(1.5)
t_adj_pval <- subset(t, P.Value < 0.05 & abs(logFC) > log2(1.5))

#print the DE genes in IFOSFAMIDE
print(nrow(t_adj_pval))

#write the csv to the folder
write.csv(t,'/projectnb/bf528/users/van-gogh/project_3/Analyst/limma_ifosfamide.csv')





