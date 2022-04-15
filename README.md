# AmrPlusPlus_SNP

## Introduction
The **AMRPlusPlus_SNP** repository contains multiple progams used to either extract SNP info or to verify ARGs for resistant-conferring SNPs. Contained further down below are descriptions to each folder contained in this repo.

## Lists of Folders in main
### *AMR_hits* Folder
Collection of csv files that contain the count matrix to all **AMRPlusPlus** hits from the SAM files in the *SAM_files* folder.

### *extracted_SNP_files* Folder
All SNP information extracted from **MetaMARC** in `MetamarcInfoExtraction` and **KARGVA** in `KargvaInfoExtracion`.

### *fluro_filter* Folder
Contains the `fluro_filter` program that filters SAM files for Fluoroquinolones ARGs.  
Is currently deprecated.

### *KargvaInfoExtraction* Folder
Contains the `KargvaInfoExtracion` program that extracts SNP information retrieved from the database in [the **KARGVA** repository](https://github.com/DataIntellSystLab/KARGVA).  
Currently only contains a copy of the database.
