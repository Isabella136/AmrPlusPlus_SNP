# As a Stand-Alone Program
## Default Testing
When testing the program to see if it was properly installed, simply running the program with no argument from the program's directory should yield results into a `test_output/` folder which will be located in `data/`:  
```
python3 SNP_Verification.py
```
If outside the program's directory, the `-c` argument should be used to specifiy the config file:
```
python3 path/to/AmrPlusPlus_SNP/SNP_Verification.py -c path/to/AmrPlusPlus_SNP/config.ini
``` 
The files located in `data/test_output/S1_test.amr.alignment.dedup/` should have the same outputted information as the files in `data/sample_output/S1_test.amr.alignment.dedup/`
## Specifying SAM File Input
The `-i` argument is used when changing the SAM file input from what is specified in the default config file:
```
python3 SNP_Verification.py -i path/to/input.sam
```
## Specifying Output Folder
The `-o` argument is used when changing the output folder from what is specified in the default config file:
```
python3 SNP_Verification.py -o path/to/outputFolder
```
## Specifying Config File
The `-c` argument is used when applying a different config file:
```
python3 SNP_Verification.py -c newConfig.ini
```
## Getting Gene-Specific Detailed Output
By default, the SNP_Verification produces a generalized output for each of the five gene types (Normal, Intrinsic, Hypersusceptible, Frameshift, and Suppressible frameshift). To have a more detailed output for one or more SNPConfirmation genes identified in the SAM file input, type the appropriate command:
```
python3 SNP_Verification.py --detailed_output all             #For a detailed output for all SNPConfirmation genes
python3 SNP_Verification.py --detailed_output MEG_XX          #For a detailed output for MEG_XX only 
                                                              #Will also limit the generalized output to MEG_XX only
python3 SNP_Verification.py --detailed_output MEG_XX,MEG_YY   #For a detailed output for MEG_XX and MEG_YY
                                                              #Will also limit the generalized output to MEG_XX and MEG_YY
```
## Changing Verification Rules for Wild-Type-next-to-Mutant Scenarios
By default, if in the case of an insertion both the wild-type and the mutant are present next to each other, the SNP_Verification program will mark the read as resistant. To change this, use the `--mt_and_wt` argument.
```
python3 SNP_Verification.py --mt_and_wt false
```

# As a Module to AMR++ v3
## Default Testing
The newest version of the [AMR++ bioinformatics pipeline](https://github.com/Microbial-Ecology-Group/AMRplusplus) makes use of the SNP_Verification program.  
When testing the AMR++ pipeline, SNP_Verification is called upon three times, once for each test input file.  
The command used by nextflow is listed below, where 
- the `-a` argument specify that the program is used by the AMR++ pipeline and therefore should recalculate the count matrix outputted by the `resistomeresults` process
- the `-i` argument calls upon one of three SAM file that will be used as input (all three of which are located in the `data\raw\` folder in the AMR++ v3 pipeline)
- the `-o` argument declares the directory where the SNP_Verification-specific output will be saved
- the `--count_matrix` argument will specify the count matrix that needs to be updated with the SNP_Verification output.
Both the output folder and the count matrix should be found in `test_results/ResistomeResults`  
```
python3 $amrsnp/SNP_Verification.py -c $amrsnp/config.ini -a -i ${sam} -o ${sample_id}_SNPs --count_matrix ${raw_count_matrix_snp}
```
## Specifying Sample Input
When used in the AMR++ ppeline, the SAM input to the SNP_Verification program is produced by the `bwa_align` process and is based on the samples specified in the `--reads` parameter when running `main_AMR++.nf`:
```
nextflow run main_AMR++.nf -profile conda --pipeline standard_AMR --reads "path/to/your/reads/*_R{1,2}.fastq.gz"
```
## Specifying Path to `ResistomeResults/`
When used in the AMR++ pipeline, you can specify the path to `ResistomeResults` (and other outputs given by AMR++) by specifying the `--output` parameter when running `main_AMR++.nf`:
```
nextflow run main_AMR++.nf -profile conda --pipeline standard_AMR --output "path/to/output/"
```
