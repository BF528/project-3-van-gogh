# Monil Gandhi - Analyst

#File Description - Concordance analysis between RNA-seq and microarray

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

#load deseq data
deseq_leflunomide <- read.csv("/projectnb/bf528/users/van-gogh/project_3/Programmer/deseq_results2.csv")
deseq_fluconazole <- read.csv("/projectnb/bf528/users/van-gogh/project_3/Programmer/deseq_results3.csv")
deseq_ifosfamide <- read.csv("/projectnb/bf528/users/van-gogh/project_3/Programmer/deseq_results1.csv")

example_limma_microarray <- read.csv("/project/bf528/project_3/results/example_limma_results.csv")

#load rna seq  data
example_rna_deseq <- read.csv("/project/bf528/project_3/results/example_deseq_results.csv")


map_affy <- read.csv("/project/bf528/project_3/refseq_affy_map.csv")

#get the probe ids
example_limma_id <- example_limma_microarray$X
lefl_probe_id <- limma_leflunomide$X
fluco_probe_id <- limma_fluconazole$X
ifo_probe_id <- limma_ifosfamide$X

#map the probes to refseqs
assign_refseq <- function(map_affy, sample_df, probe_ids_df){
  for (num in 1:length(probe_ids_df))
  {
    refseq_ids <- list(map_affy$REFSEQ[which(probe_ids_df[num] == map_affy$PROBEID)]) 
    sample_df$REFSEQID[num] <- refseq_ids
  }
  return(sample_df)
}


example_limma_microarray <- assign_refseq(map_affy, example_limma_microarray, example_limma_id)
limma_leflunomide <- assign_refseq(map_affy, limma_leflunomide,lefl_probe_id)
limma_fluconazole <- assign_refseq(map_affy, limma_fluconazole,fluco_probe_id)
limma_ifosfamide <- assign_refseq(map_affy, limma_ifosfamide,ifo_probe_id)

#subset the df
subset_col_df <- function(sample_df){
  empty_sample_df <- sample_df[FALSE, names(sample_df) %in% c("X", "REFSEQID", "logFC", "AveExpr")]
  return(empty_sample_df)
}

subset_example_limma <- as.data.frame(subset_col_df(example_limma_microarray))
subset_leflunomide <- as.data.frame(subset_col_df(limma_leflunomide))
subset_fluconazole <- as.data.frame(subset_col_df(limma_fluconazole))
subset_ifosfamide <- as.data.frame(subset_col_df(limma_ifosfamide))

#get the significant probes from the microarray dataset
significant_probes_genes <- function(sample_df, sample_empty_df){
  count <- 0
  for (record in 1:nrow(sample_df)){
    # if((sample_df[record, "P.Value"] < 0.05) & (abs(sample_df[record, "logFC"]) > 1.5))
    if((sample_df[record, "P.Value"] < 0.05) & (abs(sample_df[record, "logFC"]) > log2(1.5)))
    {
     if(length(sample_df[record, "REFSEQID"][[1]])  > 0)
       {
       for(refseq_id in sample_df[record, "REFSEQID"][[1]])
         {
         new_row <- data.frame(sample_df[record, "X"], refseq_id, sample_df[record, "logFC"], sample_df[record, "AveExpr"]) 
         names(new_row)<-c("X","REFSEQID", "logFC", "AveExpr") 
         sample_empty_df <- rbind(sample_empty_df, new_row)
         }
      } 
    }
    count <- count + 1
    
  }
  return(sample_empty_df)
}

subset_example_limma <- significant_probes_genes(example_limma_microarray, subset_example_limma)
subset_leflunomide <- significant_probes_genes(limma_leflunomide, subset_leflunomide)
subset_fluconazole <- significant_probes_genes(limma_fluconazole, subset_fluconazole)
subset_ifosfamide <- significant_probes_genes(limma_ifosfamide, subset_ifosfamide)


#final dataframes of the microarray analysis
subset_example_limma_grouped <- subset_example_limma %>% group_by(REFSEQID) %>% summarise(logFC=median(logFC), AveExpr=median(AveExpr))

subset_leflunomide_grouped  <- subset_leflunomide %>% group_by(REFSEQID) %>% summarise(logFC=median(logFC), AveExpr=median(AveExpr))

subset_fluconazole_grouped  <- subset_fluconazole %>% group_by(REFSEQID) %>% summarise(logFC=median(logFC), AveExpr=median(AveExpr))

subset_ifosfamide_grouped  <- subset_ifosfamide %>% group_by(REFSEQID) %>% summarise(logFC=median(logFC), AveExpr=median(AveExpr))

#get the significant genes from the rna-seq
significant_deseq_genes <- function(sample_deseq_df){
  # return(subset(sample_deseq_df, ((pvalue < 0.05) & (abs(log2FoldChange) > 1.5))))
  return(subset(sample_deseq_df, ((pvalue < 0.05) & (abs(log2FoldChange) > log2(1.5)))))
}

sig_deseq_leflunomide <- significant_deseq_genes(deseq_leflunomide)
sig_deseq_fluconazole <-  significant_deseq_genes(deseq_fluconazole)
sig_deseq_ifosfamide <-  significant_deseq_genes(deseq_ifosfamide)
example_deseq_significant <- significant_deseq_genes(example_rna_deseq)

#calculate the concordance
calculate_concordance <- function(df_deseq_rna, df_limma_microarray){
  counter = 0
  for(num in 1:nrow(df_deseq_rna)){
    matching_record <- df_limma_microarray[which(df_limma_microarray$REFSEQID == df_deseq_rna[num, "X"]), ]
    if(nrow(matching_record) != 0){
      if(sign(df_deseq_rna[num, "log2FoldChange"]) == sign(matching_record[["logFC"]])){
        counter <- counter + 1
      }
    }
  }

  return ((2 * counter)/(nrow(df_deseq_rna) + nrow(df_limma_microarray)))
}

# #Example set plots and concordance calculation
# value_example_concordance <- calculate_concordance(example_deseq_significant, subset_example_limma_grouped)
# 
# #Example plot deseq
# df_example_plot_deseq <- data.frame("Concordance_DEG" = value_example_concordance, "Treatment_Effect" = nrow(example_deseq_significant))
# 
# ggplot(df_example_plot_deseq, aes(x=Treatment_Effect, y=Concordance_DEG)) + geom_point()
# 
# #Example plot limma_microarray
# df_example_plot_microarray <- data.frame("Concordance_DEG" = value_example_concordance, "Treatment_Effect" = nrow(subset_example_limma_grouped))
# 
# ggplot(df_example_plot_microarray, aes(x=Treatment_Effect, y=Concordance_DEG)) + geom_point()



#leflunomide
value_concordance_leflunomide <- calculate_concordance(sig_deseq_leflunomide, subset_leflunomide_grouped)

#fluconazole
value_concordance_fluconazole <- calculate_concordance(sig_deseq_fluconazole, subset_fluconazole_grouped)

#ifosfamide
value_concordance_ifosfamide <- calculate_concordance(sig_deseq_ifosfamide, subset_ifosfamide_grouped)


#scatter plot for deseq
labels_chemicals <- c("LEF", "FLU", "IFO")

deseq_df <- data.frame("Concordance_DEG"=c(value_concordance_leflunomide, value_concordance_fluconazole, value_concordance_ifosfamide), 
                       "Treatment_Effect"=c(nrow(sig_deseq_leflunomide), nrow(sig_deseq_fluconazole), nrow(sig_deseq_ifosfamide)))

plot_deseq <- ggplot(deseq_df, aes(x=Treatment_Effect, y=Concordance_DEG)) +
  geom_point(color="black", fill=c('green', 'yellow', 'blue'), shape=21, size=3) +
  stat_smooth(method = "lm", linetype="dashed", se=FALSE, color="black") + 
  labs(title="RNA-Seq Plot", x="Treatment Effect Size - RNAseq", y="Concordance DEG") +
  stat_cor(aes(label=paste(..rr.label..))) +
  geom_text(label=labels_chemicals, hjust = 0, nudge_x=3.0, nudge_y=-0.005, size=4.0)



#scatter plot for microarray
labels_chemicals <- c("LEF", "FLU", "IFO")

microarray_df <- data.frame("Concordance_DEG"=c(value_concordance_leflunomide, value_concordance_fluconazole, value_concordance_ifosfamide), 
                       "Treatment_Effect"=c(nrow(subset_leflunomide_grouped), nrow(subset_fluconazole_grouped), nrow(subset_ifosfamide_grouped)))

microarray_df
plot_microarray <- ggplot(microarray_df, aes(x=Treatment_Effect, y=Concordance_DEG)) +
  geom_point(color="black", fill=c('green', 'yellow', 'blue'), shape=21, size=3) +
  stat_smooth(method = "lm", linetype="dashed", se=FALSE, color="black") + 
  labs(title="Microarray Plot", x="Treatment Effect Size - Microarray", y="Concordance DEG") +
  stat_cor(aes(label=paste(..rr.label..))) +
  geom_text(label=labels_chemicals, hjust = 0, nudge_x=3.0, nudge_y=-0.005, size=4.0)



plot_plane <- ggarrange(plot_deseq, plot_microarray,
                     labels = c("A.", "B."),
                     ncol = 2, nrow = 1)


#function call for the median of rna-seq
get_median_deseq_baseMean <- function(df_sample){
  median_val <- median(df_sample$baseMean)
  above_mean_results <- subset(df_sample, baseMean > median_val)
  below_mean_results <- subset(df_sample, baseMean < median_val)
  to_return <- list(above_mean_results,below_mean_results)
  return(to_return)
}

# #example deseq median results
# list_deseq_median_results <- get_median_deseq_baseMean(example_deseq_significant)
# df_deseq_above_median <- list_deseq_median_results[[1]]
# df_deseq_below_median <- list_deseq_median_results[[2]]

#leflunomide
list_deseq_median_results_lef <- get_median_deseq_baseMean(sig_deseq_leflunomide)
df_deseq_above_median_lef <- list_deseq_median_results_lef[[1]]
df_deseq_below_median_lef <- list_deseq_median_results_lef[[2]]

#fluconazole
list_deseq_median_results_flu <- get_median_deseq_baseMean(sig_deseq_fluconazole)
df_deseq_above_median_flu <- list_deseq_median_results_flu[[1]]
df_deseq_below_median_flu <- list_deseq_median_results_flu[[2]]

#ifosfamide
list_deseq_median_results_ifo <- get_median_deseq_baseMean(sig_deseq_ifosfamide)
df_deseq_above_median_ifo <- list_deseq_median_results_ifo[[1]]
df_deseq_below_median_ifo <- list_deseq_median_results_ifo[[2]]

#function call for the median of microarray
get_median_limma_AveExpr <- function(df_sample){
  median_val <- median(df_sample$AveExpr)
  above_mean_results <- subset(df_sample, AveExpr > median_val)
  below_mean_results <- subset(df_sample, AveExpr < median_val)
  to_return <- list(above_mean_results,below_mean_results)
  return(to_return)
}

# #example limma median results
# list_limma_median_results <- get_median_limma_AveExpr(subset_example_limma_grouped)
# df_limma_above_median <- list_limma_median_results[[1]]
# df_limma_below_median <- list_limma_median_results[[2]]

#leflunomide
list_limma_median_results_lef <- get_median_limma_AveExpr(subset_leflunomide_grouped)
df_limma_above_median_lef <- list_limma_median_results_lef[[1]]
df_limma_below_median_lef <- list_limma_median_results_lef[[2]]

#fluconazole
list_limma_median_results_flu <- get_median_limma_AveExpr(subset_fluconazole_grouped)
df_limma_above_median_flu <- list_limma_median_results_flu[[1]]
df_limma_below_median_flu <- list_limma_median_results_flu[[2]]

#ifosfamide
list_limma_median_results_ifo <- get_median_limma_AveExpr(subset_ifosfamide_grouped)
df_limma_above_median_ifo <- list_limma_median_results_ifo[[1]]
df_limma_below_median_ifo <- list_limma_median_results_ifo[[2]]

#concordance below group
#leflunomide
value_concordance_leflunomide_below <- calculate_concordance(df_deseq_below_median_lef, df_limma_below_median_lef)

#fluconazole
value_concordance_fluconazole_below <- calculate_concordance(df_deseq_below_median_flu, df_limma_below_median_flu)

#ifosfamide
value_concordance_ifosfamide_below <- calculate_concordance(df_deseq_below_median_ifo, df_limma_below_median_ifo)

#concordance above groups
#leflunomide
value_concordance_leflunomide_above <- calculate_concordance(df_deseq_above_median_lef, df_limma_above_median_lef)

#fluconazole
value_concordance_fluconazole_above <- calculate_concordance(df_deseq_above_median_flu, df_limma_above_median_flu)

#ifosfamide
value_concordance_ifosfamide_above <- calculate_concordance(df_deseq_above_median_ifo, df_limma_above_median_ifo)

# barplot 
df_barplot <- data.frame("Concordance"=c(value_concordance_leflunomide_above, value_concordance_leflunomide_below, value_concordance_leflunomide,
                                             value_concordance_fluconazole_above , value_concordance_fluconazole_below, value_concordance_fluconazole,
                                             value_concordance_ifosfamide_above, value_concordance_ifosfamide_below, value_concordance_ifosfamide),
                  "Median"=c("Above", "Below", "Total", "Above", "Below", "Total", "Above", "Below", "Total"),
                  "Sample"=c("Leflunomide","Leflunomide","Leflunomide",
                               "Fluconazole","Fluconazole","Fluconazole",
                               "Ifosfamide","Ifosfamide","Ifosfamide"))

df_barplot

ggplot(df_barplot, aes(fill=Median, x=Sample, y=Concordance)) +
  geom_bar(stat="identity", position="dodge") +
  labs(title="Above - Below Median Concordance")










