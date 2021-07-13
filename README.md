## Meta-analysis of TCGA copy number data
In this repository are the scripts used for the analysis performed in the manuscript: Somatic Copy Number Alterations in Human Cancers: An Analysis of Publically Available Data From the Cancer Genome Atlas
L. Harbers et al, Frontiers in Oncology (2021)

### Obtaining TCGA copy number data
TCGA copy number data can be obtained by running the `getTCGA_data.R` script. Set the cache (where you want to save the data) and specify where you want to save the clinical data and it will automatically download all the available copy number data and associated clinical information. 

### Running the analysis and plotting
To then run the analysis and generate the plots that are shown in the manuscript you will have to run `TCGA-cna-plots.R`. All the required data files are available in the `data/` directory of this repository. Except the TCGA copy number data, which is too large to include. But this can be downloaded with the above mentioned script. 
