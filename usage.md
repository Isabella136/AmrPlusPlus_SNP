# As a Stand-Alone Program
## Default Testing
When testing the program to see if it was properly installed, simply typing into command line `python3 SNP_Verification.py` from the program's directory should yield results into a *test_output* folder which will be located in *data*.  
If outside the program's directory, the command `python3 /path/to/AmrPlusPlus_SNP/SNP_Verification.py -c /path/to/AmrPlusPlus_SNP/config.ini` should be given instead.  
The files located in *data/test_output/S1_test.amr.alignment.dedup* should have the same outputted information as the files in *data/sample_output/S1_test.amr.alignment.dedup*

# As a Module to AMR++ v3
## Default Testing
The newest version of the [AMR++ bioinformatics pipeline](https://github.com/Microbial-Ecology-Group/AMRplusplus) makes use of the SNP_Verification program.  
When testing the AMR++ pipeline, SNP_Verification is called upon three times, once for each test input file.  
The command used by nextflow is the following: `python3 SNP_Verification.py -c -a -i ${sam} -o ${sample_id}_SNPs`
