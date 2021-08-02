# Supplementary material for *Mining conserved and divergent signals in 5' splicing site sequences* paper. 
In this repository you will find donor sequences data (9-basis long donor site sequences, and single and two-site nucleotide frequencies), the code implementing the regularized maximum entropy models (cpp and R scripts), and R-scripts to visualize and summarize two-site interaction patterns.


## Directory Guide

**data/**: Directory containing input data to fed the C++ program. 
* **data/freqs/** contains organism subfolders with single-site and two-site observed frequency input files. 
* **data/seqs/** 9 nucleotide length donor sequences in Rdata files named after each species.

**code_standalone/**: C++ code and program used to get a family of regularized maximum entropy models from observed single-site and two-site frequencies. In a typical run, models are obtained sequentially increasing the regularization level.

**code_compileNrun/**: R-scripts and C++ code (same as above) to infer multiple fitting models for different organisms at once.

**runs/**: standard output directory where fitting runs will be stored (other locations may be specified within the C++ code). Inside each organism's subfolder the **error.txt** file tracks the evolution of fitting error as the run went on. For each organism, results for different &gamma; values are stored at different sub-folders. 

**scripts/**: R scripts used to create circos diagrams and to get Table 1 of the paper used to analyse conserved signals. These scripts are intended to be used using the base directory as their working directory. This is to say, it's recommended to use **Mining5ss_SM** as the working directory for these scripts.

