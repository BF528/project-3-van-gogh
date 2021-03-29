# Project Description

A brief description of what this repository is for and what it contains

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
