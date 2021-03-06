# SNP_Verification Instructions

## For Mac

1. Download the latest version of Python 3 by using this link: https://www.python.org/downloads/macos/
    - Follow the instructions given by the installer
2. Open terminal
3. Install pysam by typing `pip install pysam` into terminal
4. Install the AMRPlusPlus_SNP repo by typing `gh repo clone https://github.com/Isabella136/AmrPlusPlus_SNP`
5. Run the SNP_Verification.py program by typing `python3 ./SNP_Verification.py -i <inputFile> -o <outputFolder>`

## For Linux

1. Open terminal
2. To download the latest version of Python3, type `sudo apt-get update` press enter, then type `sudo apt-get install python3.10 python3-pip`
3. Install pysam by typing `pip install pysam` into terminal
4. Install the AMRPlusPlus_SNP repo by typing `gh repo clone https://github.com/Isabella136/AmrPlusPlus_SNP`
5. Run the SNP_Verification.py program by typing `python3 ./SNP_Verification.py -i <inputFile> -o <outputFolder>`

## For Windows

Currently, pysam is difficult to install on Windows, therefore it is recommended to instead use a linux environment via command prompt to run the SNP_Verification program

## List of Arguments
-c: conditions for redistribution  
-h: help  
-i: input file  
-o: output folder  
-w: warranty disclaimer  
--mt_and_wt: true by default, used in case of insertion leading to presence of both mt and wt; if true, mark as resistant; if false, mark as susceptible
        
