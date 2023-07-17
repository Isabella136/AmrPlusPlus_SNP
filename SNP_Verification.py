#!/usr/bin/env python3

#   AMRPlusPlus_SNP_Verification
#   Copyright (C) 2022  Nathalie Bonin
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see https://www.gnu.org/licenses/.

from SNP_Verification_Tools import Gene
from SNP_Verification_Processes import verify
from Bio import SeqIO
import pysam, sys, getopt, os, configparser, csv
import pandas as pd
import multiprocessing
from concurrent.futures import ProcessPoolExecutor as ppe

# from argparse import ArgumentParser

# def parse_args():
#     parser = ArgumentParser()
#     parser.add_argument('--amr', default='true', choices=['true', 'false'])
#     parser.add_argument('--bam', required=True)
#     parser.add_argument('--outdir', required=True)
#     parser.add_argument('--mt_and_wt', default='true', choices=['true', 'false'])
#     parser.add_argument('--detailed_output', default='false', choices=['true', 'false'])
#     parser.add_argument('--count_matrix') #?
#     parser.add_argument('--count_matrix_final') #?
#     return parser.parse_args()

def parse_config():
# Define Command Line Arguments and Read Config
    configFile = "config.ini"
    config = configparser.ConfigParser()
    arguments_string = """

        -a: amrplusplus; is either 'true' or 'false'
        -c: config file; if this argument is used, must be the first listed
        -h: help
        -i: BAM input file
        -l: license disclaimer
        -o: main output folder
        -r: conditions for redistribution
        -t: threads for multiprocessing
        
        --mt_and_wt:            true by default, used in case of insertion leading to presence of both mt and wt; if true, mark as resistant; if false, mark as susceptible
        --detailed_output:      false by default, determines whether a more detailed output will be given; can be either 'false', 'all', or include a list of accessions seperated by commas
        --count_matrix:         count matrix that will be updated if amrplusplus is true
        --count_matrix_final:   the file where the updated count matrix will be found if amrplusplus is true
        
        """
    argList = []
    try:
        options, args = getopt.getopt(sys.argv[1:], "hlra:c:i:o:t:", ["mt_and_wt=", "detailed_output=", "count_matrix=", "count_matrix_final="])
    except getopt.GetoptError:
        print("ERROR - this is the list of arguments recognized by the program:{}".format(arguments_string))
        sys.exit(-1)
    countMatrixFinal = False
    for i, (opt, arg) in enumerate(options):
        if opt == "-h":
            print("List of arguments:{}".format(arguments_string))
            sys.exit()
        elif opt == "-r":
            print("""
            
            AMRPlusPlus_SNP_Verification
            Copyright (C) 2022  Nathalie Bonin
            
            This program is free software: you can redistribute it and/or modify
            it under the terms of the GNU General Public License as published by
            the Free Software Foundation, either version 3 of the License, or
            (at your option) any later version.
            
            """)
            sys.exit()
        elif opt == "-l":
            print("""
            
            AMRPlusPlus_SNP_Verification
            Copyright (C) 2022  Nathalie Bonin
            
            This program is distributed in the hope that it will be useful,
            but WITHOUT ANY WARRANTY; without even the implied warranty of
            MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
            GNU General Public License for more details.
            
            """)
            sys.exit()
        elif opt == "-c":
            if i > 0:
                print("ERROR: config file must be the first argument listed")
                sys.exit(-1)
            configFile = arg
            config.read(configFile)
        elif opt == "-i":
            if i == 0:
                config.read(configFile)
            config['FULL_FILE_NAMES']['SAMPLE'] = os.path.basename(arg).split('.')[0]
            config['SOURCE_FILES']['BAM_INPUT'] = arg
        elif opt == "-o":
            if i == 0:
                config.read(configFile)
            config['FOLDERS']['MAIN_OUTPUT_FOLDER'] = arg
        elif opt == "-a":
            if i == 0:
                config.read(configFile)
            if arg not in ['true', 'false']:
                print("ERROR: '-a' argument can only be either 'true' or 'false'")
                sys.exit(-1)
            config['SETTINGS']['AMRPLUSPLUS'] = arg
        elif opt == "-t":
            if i == 0:
                config.read(configFile)
            config['SETTINGS']['THREADS'] = arg
        elif opt == "--mt_and_wt":
            if i == 0:
                config.read(configFile)
            if arg not in ['true', 'false']:
                print("ERROR: mt_and_wt argument can only be either 'true' or 'false'")
                sys.exit(-1)
            config['SETTINGS']['MT_AND_WT'] = arg
        elif opt =="--detailed_output":
            if i == 0:
                config.read(configFile)
            if arg != "false":
                config['SETTINGS']['DETAILED'] = "true"
                if arg != "all":
                    argList = arg.split(',')
        elif opt == "--count_matrix":
            if i == 0:
                config.read(configFile)
            if not(countMatrixFinal):
                config['OUTPUT_FILES']['COUNT_MATRIX_FINAL'] = arg
            config['SOURCE_FILES']['COUNT_MATRIX'] = arg
        elif opt == "--count_matrix_final":
            if i == 0:
                config.read(configFile)
            config['OUTPUT_FILES']['COUNT_MATRIX_FINAL'] = arg
            countMatrixFinal = True
        i += 1
    if len(options) == 0:
        config.read(configFile)

    # Add new information to config buffer
    config['FOLDERS']['SAMPLE_OUTPUT'] = (config['FOLDERS']['MAIN_OUTPUT_FOLDER'] + 
                                        config['FULL_FILE_NAMES']['SAMPLE']) + '/'
    config['FOLDERS']['SAMPLE_DETAILED_OUTPUT'] = (config['FOLDERS']['SAMPLE_OUTPUT'] + 
                                                config['FOLDERS']['DETAILED_FOLDER'])
    config['FOLDERS']['TEMP'] = os.path.dirname(config['TEMP_FILES']['TEMP_BAM_SORTED'])

    config['FULL_FILE_NAMES']['COUNT_MATRIX_FINAL'] = (config['FOLDERS']['MAIN_OUTPUT_FOLDER'] + 
                                                    config['OUTPUT_FILES']['COUNT_MATRIX_FINAL'])
    config['FULL_FILE_NAMES']['NTYPE_OUTPUT'] = (config['FOLDERS']['SAMPLE_OUTPUT'] + 
                                                config['OUTPUT_FILES']['NORMAL_TYPE_OUTPUT'])
    config['FULL_FILE_NAMES']['FTYPE_OUTPUT'] = (config['FOLDERS']['SAMPLE_OUTPUT'] + 
                                                config['OUTPUT_FILES']['FRAMESHIFT_TYPE_OUTPUT'])
    config['FULL_FILE_NAMES']['HTYPE_OUTPUT'] = (config['FOLDERS']['SAMPLE_OUTPUT'] + 
                                                config['OUTPUT_FILES']['HYPERSUSCEPTIBLE_TYPE_OUTPUT'])
    config['FULL_FILE_NAMES']['STYPE_OUTPUT'] = (config['FOLDERS']['SAMPLE_OUTPUT'] + 
                                                config['OUTPUT_FILES']['SUPPRESSIBLE_TYPE_OUTPUT'])
    config['FULL_FILE_NAMES']['ITYPE_OUTPUT'] = (config['FOLDERS']['SAMPLE_OUTPUT'] + 
                                                config['OUTPUT_FILES']['INTRINSIC_TYPE_OUTPUT'])
    return (config, argList)

def dir_check(config):
    # Verify existance of folders 
    if not(os.path.exists(config['FOLDERS']['SAMPLE_DETAILED_OUTPUT'])):
        os.makedirs(config['FOLDERS']['SAMPLE_DETAILED_OUTPUT'])
    if not(os.path.exists(config['FOLDERS']['TEMP'])):
        os.makedirs(config['FOLDERS']['TEMP'])

def parse_snp_info(config):
    gene_dict = dict()
    # Get Variants Info
    for gene in SeqIO.parse(config['SOURCE_FILES']['SNP_INFO_FASTA'], 'fasta'):
        # Find index of last pipe ('|') before variant list
        index = -1
        for pipe in range(0, 5):
            index = gene.name[index+1:].find('|') + index + 1
        name = gene.name[:index]
        variants = gene.name[index+1:]

        # Store information required for creating Gene object
        gene_dict[name + "|RequiresSNPConfirmation"] = [name, gene.seq, variants]
    return gene_dict

def iterate(process_vars):

    gene_name = process_vars[0]
    gene_object = process_vars[1]
    config = process_vars[2]
    argList = process_vars[3]
    lock = process_vars[4]

    with pysam.AlignmentFile(config['TEMP_FILES']['TEMP_BAM_SORTED'], "r") as samfile:
        alignment_iterator = list(samfile.fetch(reference=gene_name))

    # Create Gene object
    if 'Must:' in gene_object[2]:
        if 'Nuc:' in gene_object[2]:            gene_variant = Gene.IntrinsicrRNA(gene_object[0],gene_object[1],gene_object[2])
        else:                                   gene_variant = Gene.IntrinsicProtein(gene_object[0],gene_object[1],gene_object[2])
    elif 'Hyper:' in gene_object[2]:            gene_variant = Gene.Hypersusceptible(gene_object[0],gene_object[1],gene_object[2])
    elif 'FS-' in gene_object[2]:
        if 'suppression' in gene_object[2]:     gene_variant = Gene.Suppressible(gene_object[0],gene_object[1],gene_object[2])
        elif 'MEG_6142' in gene_name:           gene_variant = Gene.NormalProtein(gene_object[0],gene_object[1],gene_object[2])
        else:                                   gene_variant = Gene.Frameshift(gene_object[0],gene_object[1],gene_object[2])
    else:
        if 'Nuc:' in gene_object[2]:            gene_variant = Gene.NormalrRNA(gene_object[0],gene_object[1],gene_object[2])
        else:                                   gene_variant = Gene.NormalProtein(gene_object[0],gene_object[1],gene_object[2])

    # Go through all alignments to gene
    DEBUGGING_MODE = config.getboolean('SETTINGS', 'DEBUGGING_MODE')
    if DEBUGGING_MODE:
        lock.acquire()
    for read in alignment_iterator:
        if (read.cigarstring == None):
            continue
        elif (len(argList) != 0) and (gene_variant.getName().split("|")[0] not in argList):
            continue
        verify(read, gene_variant, config)
        gene_variant.resetForNextRead()

    if DEBUGGING_MODE:
        lock.release()
    
    return {gene_name : gene_variant}

def process_genes(config, argList, gene_dict):
    processes = list()

    pysam.sort("-o", config['TEMP_FILES']['TEMP_BAM_SORTED'], config['SOURCE_FILES']['BAM_INPUT'])
    pysam.index(config['TEMP_FILES']['TEMP_BAM_SORTED'], config['TEMP_FILES']['TEMP_BAM_SORTED']+'.bai')

    m = multiprocessing.Manager()
    lock = m.Lock()

    for gene_name, gene_object in gene_dict.items():
        processes.append((gene_name, gene_object, config, argList, lock))

    with ppe(int(config['SETTINGS']['THREADS'])) as p:
        results = p.map(iterate, processes)

    return results


# Function that appends gene.getOutputInfo() to output file
def appendGeneOutputInfo(name, output_info, csvwriter, countMatrix, config):
    new_row = [name]
    new_row.extend(output_info)
    csvwriter.writerow(new_row)

    # Update count matrix if previously found in AMR++
    if config.getboolean('SETTINGS', 'AMRPLUSPLUS'):
        if name in list(countMatrix['gene_accession'].values):
            index = countMatrix['gene_accession'][countMatrix['gene_accession']==name].index[0]
            prevResCount = countMatrix.loc[index, config['FULL_FILE_NAMES']['SAMPLE']]
            if prevResCount != 0:
                newCount = output_info[1]
                newCount += (output_info[2] + output_info[5]) if len(output_info) == 10 else 0  #If gene is intrinsic 
                                                                                                #Adds 'some' and 'acquired' counts                 
                countMatrix.loc[index, config['FULL_FILE_NAMES']['SAMPLE']] = newCount


def create_output(config, argList, gene_variant_dict):
    # Create output files and write headers
    with (open(config['FULL_FILE_NAMES']['NTYPE_OUTPUT'], "w") as outputN, 
        open(config['FULL_FILE_NAMES']['FTYPE_OUTPUT'], "w") as outputF, 
        open(config['FULL_FILE_NAMES']['HTYPE_OUTPUT'], "w") as outputH,
        open(config['FULL_FILE_NAMES']['STYPE_OUTPUT'], "w") as outputS,
        open(config['FULL_FILE_NAMES']['ITYPE_OUTPUT'], "w") as outputI):

        nf_header = ["Gene","Number of reads", "Resistant", "Missense",
                    "Insertion", "Deletion", "Previously recorded nonsense",
                    "N-tuple", "Nonstop", "12+bp indel", "12+ bp frameshift",
                    "Newly found nonsense", "Frameshift till end"]

        h_header = nf_header.copy()
        h_header.extend(['Hypersusceptible mutations + resistance-conferring mutations'])

        s_header = nf_header.copy()
        s_header.pop()
        s_header.extend(["Frameshift at end", "Suppressible frameshift at res 531", 
                        "Frameshift at res 531 that is not suppressible"])

        i_header = ["Gene", "Number of reads", "All", "Some", "None", "Mutations", "Acquired", 
                    "12+bp indel", "12+ bp frameshift", "Nonsense", "Frameshift till end"]

        nwriter = csv.writer(outputN, delimiter = ',')
        fwriter = csv.writer(outputF, delimiter= ',')
        hwriter = csv.writer(outputH, delimiter= ',')
        swriter = csv.writer(outputS, delimiter= ',')
        iwriter = csv.writer(outputI, delimiter= ',')

        nwriter.writerow(nf_header)
        fwriter.writerow(nf_header)
        hwriter.writerow(h_header)
        swriter.writerow(s_header)
        iwriter.writerow(i_header)

        # Retrieve count matrix data
        countMatrix = pd.read_csv(config['SOURCE_FILES']['COUNT_MATRIX']) if config.getboolean('SETTINGS', 'AMRPLUSPLUS') else None

        # Run through results
        for name, gene in gene_variant_dict.items():
            if (len(argList) != 0) and (gene.getName().split("|")[0] not in argList):
                continue
            
            # Add to correct output
            tag = gene.getGeneTag()
            if   tag == 'N': appendGeneOutputInfo(name, gene.getOutputInfo(), nwriter, countMatrix, config)
            elif tag == 'F': appendGeneOutputInfo(name, gene.getOutputInfo(), fwriter, countMatrix, config)
            elif tag == 'H': appendGeneOutputInfo(name, gene.getOutputInfo(), hwriter, countMatrix, config)
            elif tag == 'S': appendGeneOutputInfo(name, gene.getOutputInfo(), swriter, countMatrix, config)
            else:            appendGeneOutputInfo(name, gene.getOutputInfo(), iwriter, countMatrix, config)

            # Print more detailed output if requested
            if config.getboolean('SETTINGS', 'DETAILED') and (gene.getOutputInfo()[0] > 0):
                with open(config['FOLDERS']['SAMPLE_DETAILED_OUTPUT'] + name + ".csv", "w") as detailedOutput:
                    gene.writeAdditionalInfo(detailedOutput)
            
            gene.clearOutputInfo()

    # Print count matrix
    if config.getboolean('SETTINGS', 'AMRPLUSPLUS'): countMatrix.to_csv(config['FULL_FILE_NAMES']['COUNT_MATRIX_FINAL'], index=False)
    sys.exit(0)


def main():
    config, argList = parse_config()
    dir_check(config)
    gene_dict = parse_snp_info(config)
    results = process_genes(config, argList, gene_dict)

    gene_variant_dict = dict()
    [gene_variant_dict.update(r) for r in results]

    create_output(config, argList, gene_variant_dict)

if __name__ == '__main__':
    main()
