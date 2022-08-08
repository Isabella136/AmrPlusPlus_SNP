# AmrPlusPlus_SNP

## Introduction
The **AMRPlusPlus_SNP** repository contains multiple progams used to either extract SNP info or to verify ARGs for resistant-conferring SNPs. 

## SNP_Verification Instructions

### For Mac and Linux

1. If not already done so, install latest version of Python3 and of the pysam module
2. Open terminal
3. Install the AMRPlusPlus_SNP repo by typing `gh repo clone https://github.com/Isabella136/AmrPlusPlus_SNP`
4. Run the SNP_Verification.py program by typing `python3 ./SNP_Verification.py -i <inputFile> -o <outputFolder>`

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

### *extracted_SNP_files* Folder
All SNP information extracted from `SNPInfoExtraction` (see `in-depth` branch).

### *SNP_Verification_Process* and *SNP_Verification_Tools* Folder
Contains the `SNP_Verification` program that verifies SAM files (such as the ones in the *SAM_files* folder) for resistance-conferring SNPs.  
Uses all of the SNP information saved in the *extracted_SNP_files* folder.  
Is currently in the Testing phase

