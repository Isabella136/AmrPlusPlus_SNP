# AmrPlusPlus_SNP

## Introduction
The **AMRPlusPlus_SNP** repository contains multiple progams used to either extract SNP info or to verify ARGs for resistant-conferring SNPs. 

## SNP_Verification Installation Instructions

### For Mac and Linux

1. If not already done so, install latest version of Python3 and of the pysam and pandas modules
2. Open terminal
3. Install the AMRPlusPlus_SNP repo by typing `gh repo clone https://github.com/Isabella136/AmrPlusPlus_SNP`
4. Run the SNP_Verification.py program by typing `python3 SNP_Verification.py`

### For Windows

Currently, pysam is difficult to install on Windows, therefore it is recommended to instead use a linux environment via command prompt to run the SNP_Verification program

## List of Arguments
-a: amrplusplus; is either true or false  
-c: config file; if this argument is used, must be the first listed  
-h: help  
-i: input file  
-l: license disclaimer  
-o: output folder  
-r: conditions for redistribution  

--mt_and_wt:            true by default, used in case of insertion leading to presence of both mt and wt; if true, mark as resistant; if false, mark as susceptible  
--detailed_output:      false by default, determines whether a more detailed output will be given; can be either false, 'all', or include a list of accessions seperated by commas  
--count_matrix:         count matrix that will be updated if amrplusplus is true  
--count_matrix_final:   the file where the updated count matrix will be found if amrplusplus is true

## Lists of Folders in Main

### *data* Folder
Contains all SNP information extracted from `SNPInfoExtraction` (see `in-depth` branch), the megares database, and sample input/output.

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

### *SNP_Verification_Process* and *SNP_Verification_Tools* Folder
Contains the `SNP_Verification` program that verifies SAM files for resistance-conferring SNPs.  
Uses all of the SNP information saved in the *data* folder.  

## More Information
- [Genes Issues](https://github.com/Isabella136/AmrPlusPlus_SNP/blob/main/genes-issues.md)
- [Usage](https://github.com/Isabella136/AmrPlusPlus_SNP/blob/main/usage.md)
- [SNP Info Guide](https://github.com/Isabella136/AmrPlusPlus_SNP/blob/main/data/SNPInfoGuide.md)
