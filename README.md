# Project Description

Wang et al present a large study comparing microarray and RNA-Seq gene expression data from a set of toxicological treatments with known mechanism of action (MOA) measured in rat liver. The goal of the study was to characterize the concordance of differential gene expression across platforms, test and compare how effective each platform is at detecting expected pathway-level effects based on each treatment’s MOA, and assess the MOA prediction accuracy of each platform using a test set.

In this project, we will reproduce the results from Figure 2a and Figure 3b+c, as well as compare the pathway enrichment results reported in the paper for toxgroup 3. In doing so, we will:
* Align short reads to the rat genome using STAR and quantify expression using a read counting strategy
* Perform differential expression analysis of RNA-Seq data with DESeq2
* Perform differential expression analysis of pre-normalized microarray expression data using limma
* Map between Affymetrix and refSeq identifier systems

Wang, Charles, Binsheng Gong, Pierre R. Bushel, Jean Thierry-Mieg, Danielle Thierry-Mieg, Joshua Xu, Hong Fang, et al. 2014. “A comprehensive study design reveals treatment- and transcript abundance–dependent concordance between RNA-seq and microarray data” Nature Biotechnology 32 (9): 926–32. PMID: 4243706

# Contributors
* Data Curator: Lindsay Wang (@LindsayW007) 
* Programmer: Andrew Gjelsteen (@agjelste)
* Analyst: Monil Gandhi (@gandhimonil9823)
* Biologist: Elysha Sameth (@esameth)

# Repository Contents
## Data Curator
### STAR.qsub
* Dependencies: STAR aligner
* Execution: `qsub STAR.qsub`
* Outputs: STAR outputs with `.bam` file and corresponding alignment statistics

### multiqc.qsub
* Dependencies: multiqx
* Execution: `qsub multiqc.qsub`
* Inputs: STAR outptus
* Outputs: multiqc outputs

## Programmer
### featureCounts.qsub
* Dependencies: multiqc, python 2.7
* Execution: 'qsub featureCounts.qsub'
* Inputs: 9 bam files
* Outputs: 9 count files, featureCounts .html with plots and statistics

### ProgrammerScript.R
* Dependencies: R, ggplot2, DESeq2
* Inputs: 9 count files, corresponding to each bam input
* Outputs: 3 csv files for each chemical sample, histogram plot for counts, scatterplot for counts

## Analyst
### analyst_limma.R	
* Dependencies: R
* Inputs: Two files - group_3_toxo_info.csv and liver-normalization-rma.txt
* Outputs: 3 csv files for each chemical sample 

### concordance.R	
* Dependencies: R
* Inputs: Three files (for each MOA) from `limma` and `DESeq2`
* Outputs: Histograms and Barplots

## Biologist
### GSE_filter.R
* Dependencies: R
* Inputs: Three files (for each MOA) from `limma` and `DESeq2`
* Outputs: List of probe IDs/gene names that are differentialy expressed using `| log2FC | > log2(1.5)` and `unadjusted p < 0.05` for GSEA

### heatmap.R
* Dependencies: R, BiocManager
* Inputs: Three files (for each MOA) of the normalized counts from `DESeq2`
* Outputs: Heatmap-based hierarchical clustering of the MOAs along with a PCA of the samples
