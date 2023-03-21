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
from SNP_Verification_Tools import gene_dict
from Bio import SeqIO
from multiprocessing import Process
import SNP_Verification_Tools
import pysam, sys, getopt, os, configparser, shutil, numpy, csv, multiprocessing
import pandas as pd



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
    
    --mt_and_wt:            true by default, used in case of insertion leading to presence of both mt and wt; if true, mark as resistant; if false, mark as susceptible
    --detailed_output:      false by default, determines whether a more detailed output will be given; can be either 'false', 'all', or include a list of accessions seperated by commas
    --count_matrix:         count matrix that will be updated if amrplusplus is true
    --count_matrix_final:   the file where the updated count matrix will be found if amrplusplus is true
    
    """
argList = []
try:
    options, args = getopt.getopt(sys.argv[1:], "hlrac:i:o:", ["mt_and_wt=", "detailed_output=", "count_matrix=", "count_matrix_final="])
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
        config['SETTINGS']['AMRPLUSPLUS'] = 'true'
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
                                      config['FULL_FILE_NAMES']['SAMPLE'])
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

# Verify existance of folders 
if not(os.path.exists(config['FOLDERS']['SAMPLE_DETAILED_OUTPUT'])):
    os.makedirs(config['FOLDERS']['SAMPLE_DETAILED_OUTPUT'])
if not(os.path.exists(config['FOLDERS']['TEMP'])):
    os.makedirs(config['FOLDERS']['TEMP'])

# Get Variants Info
for gene in SeqIO.parse(config['SOURCE_FILES']['SNP_INFO_FASTA'], 'fasta'):
    # Find index of last pipe ('|') before variant list
    index = -1
    for pipe in range(0, 5):
        index = gene.name[index+1:].find('|') + index + 1
    name = gene.name[:index]
    variants = gene.name[index+1:]

    # Create Gene object, add to dictionary
    if 'Must:' in variants:
        if 'Nuc:' in variants:          gene_dict[name + "|RequiresSNPConfirmation"] = Gene.IntrinsicrRNA(name, gene.seq, variants)
        else:                           gene_dict[name + "|RequiresSNPConfirmation"] = Gene.IntrinsicProtein(name, gene.seq, variants)
    elif 'Hyper:' in variants:          gene_dict[name + "|RequiresSNPConfirmation"] = Gene.Hypersusceptible(name, gene.seq, variants)
    elif 'FS-' in variants:
        if 'suppression' in variants:   gene_dict[name + "|RequiresSNPConfirmation"] = Gene.Suppressible(name, gene.seq, variants)
        elif 'MEG_6142' in name:        gene_dict[name + "|RequiresSNPConfirmation"] = Gene.NormalProtein(name, gene.seq, variants)
        else:                           gene_dict[name + "|RequiresSNPConfirmation"] = Gene.Frameshift(name, gene.seq, variants)
    else:
        if 'Nuc:' in variants:          gene_dict[name + "|RequiresSNPConfirmation"] = Gene.NormalrRNA(name, gene.seq, variants)
        else:                           gene_dict[name + "|RequiresSNPConfirmation"] = Gene.NormalProtein(name, gene.seq, variants)

# Define processes for the analysis of BAM file sorted by MEGARes reference
#pysam.sort("-o", config['TEMP_FILES']['TEMP_BAM_SORTED'], config['SOURCE_FILES']['BAM_INPUT'])
#pysam.index(config['TEMP_FILES']['TEMP_BAM_SORTED'], config['TEMP_FILES']['TEMP_BAM_SORTED']+'.bai')
with pysam.AlignmentFile(config['TEMP_FILES']['TEMP_BAM_SORTED'], "r") as samfile:
    def iterate(alignment_iterator, gene):
        for read in alignment_iterator:
            if (read.cigarstring == None):
                continue
            elif (len(argList) != 0) and (gene.getName().split("|")[0] not in argList):
                continue
            verify(read, gene, config)
            gene.resetForNextRead()
    processes = list()
    for gene in gene_dict:
        iter = samfile.fetch(reference=gene)
        processes.append(Process(target=iterate, args=(list(iter),gene_dict[gene])))
    for startIndex in range(0, len(processes), 12):
        endIndex = (startIndex + 12) if len(processes) - startIndex >= 12 else len(processes)
        for currentIndex in range(startIndex, endIndex):
            processes[currentIndex].start()
        for currentIndex in range(startIndex, endIndex):
            processes[currentIndex].join()

# Function that appends gene.getOutputInfo() to output file
def appendGeneOutputInfo(name, output_info, csvwriter, countMatrix):
    csvwriter.writerow([name].extend(output_info))

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


# Create output files and write headers
with (open(config['FULL_FILE_NAMES']['NTYPE_OUTPUT'], "w") as outputN, 
      open(config['FULL_FILE_NAMES']['FTYPE_OUTPUT'], "w") as outputF, 
      open(config['FULL_FILE_NAMES']['HTYPE_OUTPUT'], "w") as outputH,
      open(config['FULL_FILE_NAMES']['STYPE_OUTPUT'], "w") as outputS,
      open(config['FULL_FILE_NAMES']['ITYPE_OUTPUT'], "w") as outputI):

    nf_header = ["Gene,Number of reads", "Resistant", "Missense",
                 "Insertion", "Deletion", "Previously recorded nonsense",
                 "N-tuple", "Nonstop", "12+bp indel", "12+ bp frameshift",
                 "Newly found nonsense", "Frameshift till end"]

    h_header = nf_header.copy()
    h_header = h_header.append('Hypersusceptible mutations + resistance-conferring mutations')

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
    countMatrix = pd.read_csv(config['FULL_FILE_NAMES']['COUNT_MATRIX']) if config.getboolean('SETTINGS', 'AMRPLUSPLUS') else None

    # Run through results
    for name, gene in gene_dict.items():
        if (len(argList) != 0) and (gene.getName().split("|")[0] not in argList):
            continue
        
        # Add to correct output
        tag = gene.getGeneTag()
        if   tag == 'N': appendGeneOutputInfo(name, gene.getOutputInfo(), nwriter, countMatrix)
        elif tag == 'F': appendGeneOutputInfo(name, gene.getOutputInfo(), fwriter, countMatrix)
        elif tag == 'H': appendGeneOutputInfo(name, gene.getOutputInfo(), hwriter, countMatrix)
        elif tag == 'S': appendGeneOutputInfo(name, gene.getOutputInfo(), swriter, countMatrix)
        else:            appendGeneOutputInfo(name, gene.getOutputInfo(), iwriter, countMatrix)

        # Print more detailed output if requested
        if config.getboolean('SETTINGS', 'DETAILED') and (gene.getOutputInfo()[0] > 0):
            with open(config['FOLDERS']['SAMPLE_DETAILED_OUTPUT'] + name + ".csv", "w") as detailedOutput:
                gene.writeAdditionalInfo(detailedOutput)
        
        gene.clearOutputInfo()

# Print count matrix
if config.getboolean('SETTINGS', 'AMRPLUSPLUS'): countMatrix.to_csv(config['FULL_FILE_NAMES']['COUNT_MATRIX_FINAL'], index=False)
sys.exit(0)