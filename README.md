# AmrPlusPlus_SNP

## Introduction
The **AMRPlusPlus_SNP** repository contains multiple progams used to either extract SNP info or to verify ARGs for resistant-conferring SNPs. Contained further down below are descriptions to each folder contained in this repo.

## Lists of Folders in Main
### *AMR_hits* Folder
Collection of csv files that contain the count matrix to all **AMRPlusPlus** hits from the SAM files in the *SAM_files* folder.

### *extracted_SNP_files* Folder
All SNP information extracted from **MetaMARC** in `MetamarcInfoExtraction` and **KARGVA** in `KargvaInfoExtracion`.

### *fluro_filter* Folder
Contains the `fluro_filter` program that filters SAM files for Fluoroquinolones ARGs.  
Is currently deprecated.

### *KargvaInfoExtraction* Folder
Contains the `KargvaInfoExtracion` program that extracts SNP information retrieved from the database in the [KARGVA](https://github.com/DataIntellSystLab/KARGVA) repository.  
Currently only contains a copy of the database.

### *mapping_files* Folder
Contains two csv files: the first file is used to map the MEGARes v1 database headers (used by **MetaMARC**) to headers used by external databases, and the second file is used to map the MEGARes v2 database headers (used by **AMRPlusPlus_SNP**) to headers used by external databases.

### *MetamarcInfoExtraction* Folder
Contains the `MetamarcInfoExtraction` program that extracts SNP information retrieved from the database in the [MetaMARC](https://github.com/lakinsm/meta-marc) repository and saves it in the *extracted_SNP_files* folder.  
Also contains the *metamarc_files* folder which includes the folowing files: 
- "mmarc_model_members.csv": lists mmarc models and corresponding ARGs; ARGs are represented by their MEGARes v1 database header
- "mmarc_snpsearch_metadata2.txt": lists of SNP information for each mmarc model 
- "mmarc_snpsearch_metadata2_modified.txt": another version of the previous text file that has been modified for use by the `MetamarcInfoExtraction` program

### *SAM_files* Folder
All SAM files provided by the FDA for SNP verification

### *SNP_Verification* Folder
Contains the `SNP_Verification` program that verifies SAM files (such as the ones in the *SAM_files* folder) for resistance-conferring SNPs.  
Uses all of the SNP information saved in the *extracted_SNP_files* folder.  
Is currently in the Testing phase

### *Test* Folder
Where output to the `SNP_Verification` program is currently saved
