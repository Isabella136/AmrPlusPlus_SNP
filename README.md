# AmrPlusPlus_SNP

## Introduction
The **AMRPlusPlus_SNP** repository contains multiple progams used to either extract SNP info or to verify ARGs for resistant-conferring SNPs. 

## SNP_Verification Instructions

### For Mac

1. Download the latest version of Python 3 by using this link: https://www.python.org/downloads/macos/
    - Follow the instructions given by the installer
2. Open terminal
3. Install pysam by typing `pip install pysam` into terminal
4. Install the AMRPlusPlus_SNP repo by typing `gh repo clone https://github.com/Isabella136/AmrPlusPlus_SNP`
5. Run the SNP_Verification.py program by typing `python3 ./SNP_Verification.py -i <inputFile> -o <outputFolder>`

### For Linux

1. Open terminal
2. To download the latest version of Python3, type `sudo apt-get update` press enter, then type `sudo apt-get install python3.10 python3-pip`
3. Install pysam by typing `pip install pysam` into terminal
4. Install the AMRPlusPlus_SNP repo by typing `gh repo clone https://github.com/Isabella136/AmrPlusPlus_SNP`
5. Run the SNP_Verification.py program by typing `python3 ./SNP_Verification.py -i <inputFile> -o <outputFolder>`

### For Windows

Currently, pysam is difficult to install on Windows, therefore it is recommended to instead use a linux environment via command prompt to run the SNP_Verification program

### List of Arguments
-c: conditions for redistribution  
-h: help  
-i: input file  
-o: output folder  
-w: warranty disclaimer  
--mt_and_wt: true by default, used in case of insertion leading to presence of both mt and wt; if true, mark as resistant; if false, mark as susceptible

## Lists of Folders in Main
### *AMR_hits* Folder
Collection of csv files that contain the count matrix to all **AMRPlusPlus** hits from the SAM files in the *SAM_files* folder.

### *extracted_SNP_files* Folder
All SNP information extracted from **MetaMARC** in `MetamarcInfoExtraction` and **KARGVA** in `KargvaInfoExtracion`.

### *InfoExtraction* Folder
Contains the following three folders:

#### *KargvaInfoExtraction* Folder
Contains the `KargvaInfoExtracion` program that extracts SNP information retrieved from the database in the [KARGVA](https://github.com/DataIntellSystLab/KARGVA) repository and saves it in the *extracted_SNP_files* folder. 

#### *LiteratureInfoExtraction* Folder
Contains the `LiteratureInfoExtracion` program that extracts SNP information retrieved from the literature (which was previously compiled into the "SNPinfo_literature.csv" file) and saves it in the *extracted_SNP_files* folder.  

#### *MetamarcInfoExtraction* Folder
Contains the `MetamarcInfoExtraction` program that extracts SNP information retrieved from the database in the [MetaMARC](https://github.com/lakinsm/meta-marc) repository and saves it in the *extracted_SNP_files* folder.  
Also contains the *metamarc_files* folder which includes the folowing files: 
- "mmarc_model_members.csv": lists mmarc models and corresponding ARGs; ARGs are represented by their MEGARes v1 database header
- "mmarc_snpsearch_metadata2.txt": lists of SNP information for each mmarc model 
- "mmarc_snpsearch_metadata2_modified.txt": another version of the previous text file that has been modified for use by the `MetamarcInfoExtraction` program

### *mapping_files* Folder
Contains two csv files: the first file is used to map the MEGARes v1 database headers (used by **MetaMARC**) to headers used by external databases, and the second file is used to map the MEGARes v2 database headers (used by **AMRPlusPlus_SNP**) to headers used by external databases.

### *SAM_files* Folder
All SAM files provided by the FDA for SNP verification

### *SNP_Verification* Folder
Contains the `SNP_Verification` program that verifies SAM files (such as the ones in the *SAM_files* folder) for resistance-conferring SNPs.  
Uses all of the SNP information saved in the *extracted_SNP_files* folder.  
Is currently in the Testing phase

### *Test* Folder
Where output to the `SNP_Verification` program is currently saved
