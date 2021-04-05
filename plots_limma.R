# Monil Gandhi - Analyst

#File Description - Import the csv from the limma analysis and plot the histograms and volcano plots for the DE genes

# install packages and load libraries
install.packages("pacman")
require(pacman)  
library(pacman)  
library(reshape2)
library(grid)
library(gridExtra)
install.packages("limma")
library(tidyverse)
library(EnhancedVolcano)
library(limma)

#run pacman for downloading the necessary libraries
pacman::p_load(pacman, dplyr, GGally, ggplot2, ggthemes, 
               ggvis, httr, lubridate, plotly, rio, rmarkdown, shiny, 
               stringr, tidyr, ggpubr, grid, gridExtra, hrbrthemes, EnhancedVolcano) 

# load limma data
limma_leflunomide <- read.csv("/projectnb/bf528/users/van-gogh/project_3/Analyst/limma_leflunomide.csv")
limma_ifosfamide <- read.csv("/projectnb/bf528/users/van-gogh/project_3/Analyst/limma_ifosfamide.csv")
limma_fluconazole <- read.csv("/projectnb/bf528/users/van-gogh/project_3/Analyst/limma_fluconazole.csv")


# get the DE gens with P.Value < 0.05 & abs(logFC) > log2(1.5) for each of the three samples

# histogram_leflunomide <- limma_leflunomide[limma_leflunomide$adj.P.Val < 0.05,]
histogram_leflunomide <- limma_leflunomide[limma_leflunomide$P.Val < 0.05 & abs(limma_leflunomide$logFC) > log2(1.5),]
# histogram_fluconazole <- limma_fluconazole[limma_fluconazole$adj.P.Val < 0.05,]
histogram_fluconazole <- limma_fluconazole[limma_fluconazole$P.Val < 0.05 & abs(limma_fluconazole$logFC) > log2(1.5),]
# histogram_ifosfamide <- limma_ifosfamide[limma_ifosfamide$adj.P.Val < 0.05,]
histogram_ifosfamide <- limma_ifosfamide[limma_ifosfamide$P.Val < 0.05 & abs(limma_ifosfamide$logFC) > log2(1.5),]

#Plot the histograms for the DE genes for each of the three samples
l <- ggplot(histogram_leflunomide, aes(x=logFC))+
  geom_histogram(color="darkblue", fill="lightblue") + labs(title="DE genes - Leflunomide", x="Log Fold Change", y="Count") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_text(size=10), axis.title.y=element_text(size=10))

f <- ggplot(histogram_fluconazole, aes(x=logFC))+
  geom_histogram(color="darkblue", fill="lightblue") + labs(title="DE genes - Fluconazole", x="Log Fold Change", y="Count") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_text(size=10), axis.title.y=element_text(size=10))

i <- ggplot(histogram_ifosfamide, aes(x=logFC))+
  geom_histogram(color="darkblue", fill="lightblue") + labs(title="DE genes - Ifosfamide", x="Log Fold Change", y="Count") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_text(size=10), axis.title.y=element_text(size=10))

# Plot the histograms on the same plane
histogram_plane <- ggarrange(l, f, i, 
                      labels = c("1.", "2.", "3."), 
                      ncol = 3, nrow = 1)

# ggsave("histograms_part5.png", plot = histogram_plane, width=7, height=3, units="in")

#Plot the volcano plot for Leflunomide
row.names(limma_leflunomide) <- limma_leflunomide$X

limma_leflunomide_volcano <- subset(limma_leflunomide,select= c("logFC","P.Value"))

l <- EnhancedVolcano(limma_leflunomide_volcano,
                lab = rownames(limma_leflunomide_volcano),
                x = 'logFC',
                y = 'P.Value',
                title = 'Fold Change vs Nominal p-value - Leflunomide',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                pointSize = 1.5,
                labSize = 4.0,
                xlim = c(-5,9), widthConnectors = 0.5)


#Plot the volcano plot for Fluconazole
row.names(limma_fluconazole) <- limma_fluconazole$X

limma_fluconazole_volcano <- subset(limma_fluconazole,select= c("logFC","P.Value"))

f <- EnhancedVolcano(limma_fluconazole_volcano,
                lab = rownames(limma_fluconazole_volcano),
                x = 'logFC',
                y = 'P.Value',
                title = 'Fold Change vs Nominal p-value - Fluconazole',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                pointSize = 1.5,
                labSize = 4.0,
                xlim = c(-4,4), widthConnectors = 0.5)

#Plot the volcano plot for Ifosfamide
row.names(limma_ifosfamide) <- limma_ifosfamide$X

limma_ifosfamide_volcano <- subset(limma_ifosfamide,select= c("logFC","P.Value"))

i <- EnhancedVolcano(limma_ifosfamide_volcano,
                lab = rownames(limma_ifosfamide_volcano),
                x = 'logFC',
                y = 'P.Value',
                title = 'Fold Change vs Nominal p-value - Ifosfamide',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                pointSize = 1.5,
                labSize = 4.0,
                xlim = c(-3,3), widthConnectors = 0.5)


# volcano_plane <- ggarrange(l, f, i, 
#                       labels = c("1.", "2.", "3."), 
#                       ncol = 3, nrow = 1)

# ggsave("volcano_leflunomide.png", plot = l, width=7, height=3, units="in")
# ggsave("volcano_fluconazole.png", plot = f, width=7, height=3, units="in")
# ggsave("volcano_ifosfamide.png", plot = i, width=7, height=3, units="in")



