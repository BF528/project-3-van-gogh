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


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")


samples <- read.csv("/project/bf528/project_3/groups/group_3_mic_info.csv")

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

design <- model.matrix(
  ~factor(
    samples$chemical[samples$vehicle=='CORN_OIL_100_%'],
    levels=c('Control','LEFLUNOMIDE')
  )
)

colnames(design) <- c('Intercept','LEFLUNOMIDE')

# run limma
fit <- lmFit(rma.subset.leflmide, design)
fit <- eBayes(fit)
t <- topTable(fit, coef=2, n=nrow(rma.subset.leflmide), adjust='BH')

# t_adj_pval <- subset(t, adj.P.Val < 0.05)

t_adj_pval <- subset(t, P.Value < 0.05 & abs(logFC) > log2(1.5))

print(nrow(t_adj_pval))

write.csv(t,'/projectnb/bf528/users/van-gogh/project_3/Analyst/limma_leflunomide.csv')



#FLUCONAZOLE
rma.subset.fluzole <- rma[paste0('X',samples$array_id[samples$chemical=='FLUCONAZOLE' | (samples$chemical=='Control' & samples$vehicle=='CORN_OIL_100_%')])]

design <- model.matrix(
  ~factor(
    samples$chemical[samples$vehicle=='CORN_OIL_100_%'],
    levels=c('Control','FLUCONAZOLE')
  )
)

colnames(design) <- c('Intercept','FLUCONAZOLE')

# run limma
fit <- lmFit(rma.subset.fluzole, design)
fit <- eBayes(fit)
t <- topTable(fit, coef=2, n=nrow(rma.subset.fluzole), adjust='BH')

# t_adj_pval <- subset(t, adj.P.Val < 0.05)

t_adj_pval <- subset(t, P.Value < 0.05 & abs(logFC) > log2(1.5))


print(nrow(t_adj_pval))
write.csv(t,'/projectnb/bf528/users/van-gogh/project_3/Analyst/limma_fluconazole.csv')


#IFOSFAMIDE
rma.subset.isomide <- rma[paste0('X',samples$array_id[samples$chemical=='IFOSFAMIDE' | (samples$chemical=='Control' & samples$vehicle=='SALINE_100_%')])]

design <- model.matrix(
  ~factor(
    samples$chemical[samples$vehicle=='SALINE_100_%'],
    levels=c('Control','IFOSFAMIDE')
  )
)

colnames(design) <- c('Intercept','IFOSFAMIDE')

# run limma
fit <- lmFit(rma.subset.isomide, design)
fit <- eBayes(fit)
t <- topTable(fit, coef=2, n=nrow(rma.subset.isomide), adjust='BH')

# t_adj_pval <- subset(t, adj.P.Val < 0.05)

t_adj_pval <- subset(t, P.Value < 0.05 & abs(logFC) > log2(1.5))


print(nrow(t_adj_pval))
write.csv(t,'/projectnb/bf528/users/van-gogh/project_3/Analyst/limma_ifosfamide.csv')





