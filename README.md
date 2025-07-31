## BGM
This repository contains the implementation of the method proposed in the paper "Bi-Gaussian Mirrors for False Discovery Rate Control." The repository is organized as follows:


## Simulations/
This folder includes the code for reproducing the simulation results presented in Figures 2–4 of the paper:

"Simu_Figure#.R": Main scripts for generating the results shown in Figure #.

"functions.R": Supporting utility functions used across the simulation scripts.


## Real data/
This folder contains two subfolders, each corresponding to a real-world dataset analyzed in the paper:

"HIV_code/" contains code for analyzing the HIV drug resistance dataset. Its scripts "HIV_#.R" reproduce the results in Figures 5–7 of the paper. Original data of this example are publicly available at: https://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/. For convenience, we also provide the dataset as "#_DATA.txt" files.

"scRNAseq_code/" contains code for analyzing the single-cell RNA sequencing (scRNA-seq) dataset. The script "scRNA_seq.R" reproduces results reported in Table 3 of the paper. Utility functions are included in "functions.R". Original data of this analysis are available at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141834. Due to file size limits, the raw ".txt" data are not included in this repository.


## Notation
Please note that the reproduced empirical results may differ slightly from those reported in the paper. This variation is expected due to the inherent randomness in the sampling process and the non-deterministic nature of the algorithm.


## Function Documentation
All utility functions are documented via in-line annotations in their respective files, with clear usage instructions and explanations of key implementation steps.
