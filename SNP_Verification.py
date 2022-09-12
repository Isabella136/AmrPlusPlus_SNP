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

from SNP_Verification_Tools.Gene import Gene
from SNP_Verification_Processes import verify
from SNP_Verification_Tools import geneDict
import SNP_Verification_Tools
import pysam, sys, getopt, os, configparser, shutil, numpy
import pandas as pd

argList = []
configFile = "config.ini"
config = configparser.ConfigParser()

try:
    options, args = getopt.getopt(sys.argv[1:], "hlrac:i:o:", ["mt_and_wt=", "detailed_output=", "count_matrix=", "count_matrix_final="])
except getopt.GetoptError:
    print("""
        ERROR - this is the list of arguments recognized by the program:

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
        
        """)
    sys.exit(-1)

i = 0
confChanged = False
output = False
countMatrixFinal = False
for opt, arg in options:
    if opt == "-h":
        print("""
        List of arguments:

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
        
        """)
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
        confChanged = True
    elif opt == "-i":
        if i == 0:
            config.read(configFile)
        config['SOURCE_FILES']['SAM_INPUT'] = arg
        if output:
            temp = config['OUTPUT_FILES']['OUTPUT_FOLDER'][:config['OUTPUT_FILES']['OUTPUT_FOLDER'].rfind('/')]
            config['OUTPUT_FILES']['OUTPUT_FOLDER'] = temp + "/" + config['SOURCE_FILES']['SAM_INPUT'].split('/')[-1].split('.')[0]
        else:
            config['OUTPUT_FILES']['OUTPUT_FOLDER'] = config['OUTPUT_FILES']['OUTPUT_FOLDER'] + "/" + config['SOURCE_FILES']['SAM_INPUT'].split('/')[-1].split('.')[0]
    elif opt == "-o":
        if i == 0:
            config.read(configFile)
        config['OUTPUT_FILES']['OUTPUT_FOLDER'] = arg + "/" + config['SOURCE_FILES']['SAM_INPUT'].split('/')[-1].split('.')[0] 
        output = True
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
    elif opt == "--count_matrix=":
        if i == 0:
            config.read(configFile)
        if not(countMatrixFinal):
            config['OUTPUT_FILES']['COUNT_MATRIX_FINAL'] = arg
        config['SOURCE_FILES']['COUNT_MATRIX'] = arg
    elif opt == "--count_matrix_final=":
        if i == 0:
            config.read(configFile)
        config['OUTPUT_FILES']['COUNT_MATRIX_FINAL'] = arg
        countMatrixFinal = True
    i += 1

if (i == 0) or ((i == 1) and confChanged):
    config.read(configFile)
else:
    newConfigFile = configFile[:configFile.rfind('.')]
    underscore = newConfigFile.rfind('_0')
    if (underscore != -1) and newConfigFile[underscore+2:].isnumeric:
        newConfigFile = newConfigFile[:underscore+2] + str(int(newConfigFile[underscore+2:]) + 1)
    else:
        newConfigFile = newConfigFile + '_01'
    newConfigFile = newConfigFile + '.ini'

    with open(newConfigFile, 'w') as configFileTemp:
        config.write(configFileTemp)

if not(os.path.exists(config['OUTPUT_FILES']['OUTPUT_FOLDER'])):
    os.makedirs(config['OUTPUT_FILES']['OUTPUT_FOLDER'])
if config.getboolean('SETTINGS', 'AMRPLUSPLUS'):
    if config['OUTPUT_FILES']['COUNT_MATRIX_FINAL'].find('/') != -1:
        if not(os.path.exists(config['OUTPUT_FILES']['COUNT_MATRIX_FINAL'][:config['OUTPUT_FILES']['COUNT_MATRIX_FINAL'].find('/')])):
            os.makedirs(config['OUTPUT_FILES']['COUNT_MATRIX_FINAL'][:config['OUTPUT_FILES']['COUNT_MATRIX_FINAL'].find('/')])
            with open(config['OUTPUT_FILES']['COUNT_MATRIX_FINAL'], 'w') as fp: pass
        elif not(os.path.exists(config['OUTPUT_FILES']['COUNT_MATRIX_FINAL'])):
            with open(config['OUTPUT_FILES']['COUNT_MATRIX_FINAL'], 'w') as fp: pass
    elif not(os.path.exists(config['OUTPUT_FILES']['COUNT_MATRIX_FINAL'])):
        with open(config['OUTPUT_FILES']['COUNT_MATRIX_FINAL'], 'w') as fp: pass
if config.getboolean('SETTINGS', 'DETAILED') and not(os.path.exists(config['OUTPUT_FILES']['OUTPUT_FOLDER'] + config['OUTPUT_FILES']['DETAILED_FOLDER'])):
    os.mkdir(config['OUTPUT_FILES']['OUTPUT_FOLDER'] + config['OUTPUT_FILES']['DETAILED_FOLDER'])
if not(os.path.exists(config['TEMP_FILES']['TEMP_SAM_SORTED'][:config['TEMP_FILES']['TEMP_SAM_SORTED'].rfind('/')])):
    os.mkdir(config['TEMP_FILES']['TEMP_SAM_SORTED'][:config['TEMP_FILES']['TEMP_SAM_SORTED'].rfind('/')])

SNPinfo = open(config['SOURCE_FILES']['SNP_INFO_FASTA'], "rt")
isSequence = False
name = ""
snp = ""
sequence = ""
for line in SNPinfo:
    if isSequence:
        sequence = line
        geneDict.update({name + "|RequiresSNPConfirmation":Gene(name, sequence[:-1], snp)})
        isSequence = False
    else:
        temp = 0
        for i in range(0, 5):
            temp = line[temp+1:].find('|') + temp + 1
        name = line[1:temp]
        snp = line[temp+1:len(line)-1]
        isSequence = True
SNPinfo.close()

    #Analyze SAM file
pysam.sort("-o", config['TEMP_FILES']['TEMP_SAM_SORTED'], config['SOURCE_FILES']['SAM_INPUT'])
samfile = pysam.AlignmentFile(config['TEMP_FILES']['TEMP_SAM_SORTED'], "r")
iter = samfile.fetch()
for read in iter:
    gene = geneDict.get(read.reference_name, False)
    if (gene == False):
        continue
    elif (read.cigarstring == None):
        continue
    elif (len(argList) != 0) and (gene.getName().split("|")[0] not in argList):
        continue
    verify(read, gene, config)
    gene.resetForNextRead()
samfile.close() 

#Function that appends gene.getOutputInfo() to output file
def appendGeneOutputInfo(name, outputInfo, file, countMatrix):
    file.write("\n" + name)
    for info in outputInfo.values():
        file.write("," + str(info))
    if type(countMatrix) == pd.DataFrame:
        if name in list(countMatrix['gene_accession'].values):
            index = countMatrix['gene_accession'][countMatrix['gene_accession']==name].index[0]
            prevResCount = countMatrix.loc[index, config['SOURCE_FILES']['SAM_INPUT'].split('/')[-1].split('.')[0]]
            if prevResCount != 0:
                newCount = outputInfo[1]
                if len(outputInfo) == 10:                                   #If gene is intrinsic
                    newCount += (outputInfo[2] + outputInfo[5])             #Adds 'some' and 'acquired' counts
                countMatrix.loc[index, config['SOURCE_FILES']['SAM_INPUT'].split('/')[-1].split('.')[0]] = newCount


#Create output files and write headers
outputN = open(config['OUTPUT_FILES']['OUTPUT_FOLDER'] + config['OUTPUT_FILES']['NORMAL_TYPE_OUTPUT'], "w")
outputN.write("Gene,Number of reads,Resistant,Missense,Insertion,Deletion,Previously recorded nonsense,N-tuple,Nonstop,12+bp indel,12+ bp frameshift,Newly found nonsense,Frameshift till end")
outputN.close()
outputF = open(config['OUTPUT_FILES']['OUTPUT_FOLDER'] + config['OUTPUT_FILES']['FRAMESHIFT_TYPE_OUTPUT'], "w")
outputF.write("Gene,Number of reads,Resistant,Missense,Insertion,Deletion,Previously recorded nonsense,N-tuple,Nonstop,12+bp indel,12+ bp frameshift,Newly found nonsense,Frameshift till end")
outputF.close()
outputH = open(config['OUTPUT_FILES']['OUTPUT_FOLDER'] + config['OUTPUT_FILES']['HYPERSUSCEPTIBLE_TYPE_OUTPUT'], "w")
outputH.write("Gene,Number of reads,Resistant,Missense,Insertion,Deletion,Previously recorded nonsense,N-tuple,Nonstop,12+bp indel,12+ bp frameshift,Newly found nonsense,Frameshift till end,Hypersusceptible mutations + resistance-conferring mutations")
outputH.close()
outputS = open(config['OUTPUT_FILES']['OUTPUT_FOLDER'] + config['OUTPUT_FILES']['SUPPRESSIBLE_TYPE_OUTPUT'], "w")
outputS.write("Gene,Number of reads,Resistant,Missense,Insertion,Deletion,Previously recorded nonsense,N-tuple,Nonstop,12+bp indel,12+ bp frameshift,Newly found nonsense,Frameshift at end ,Suppressible frameshift at res 531,Frameshift at res 531 that is not suppressible")
outputS.close()
outputI = open(config['OUTPUT_FILES']['OUTPUT_FOLDER'] + config['OUTPUT_FILES']['INTRINSIC_TYPE_OUTPUT'], "w")
outputI.write("Gene,Number of reads,All,Some,None,Mutations,Acquired,12+bp indel,12+ bp frameshift,Nonsense,Frameshift till end")
outputI.close()

countMatrix = None
if config.getboolean('SETTINGS', 'AMRPLUSPLUS'):
    if config['OUTPUT_FILES']['COUNT_MATRIX_FINAL'] != config['SOURCE_FILES']['COUNT_MATRIX']:
        shutil.copyfile(config['SOURCE_FILES']['COUNT_MATRIX'], config['OUTPUT_FILES']['COUNT_MATRIX_FINAL'])
    countMatrix = pd.read_csv(config['OUTPUT_FILES']['COUNT_MATRIX_FINAL'])

for name, gene in geneDict.items():
    if (len(argList) != 0) and (gene.getName().split("|")[0] not in argList):
        continue
    tag = gene.getGeneTag()
    if tag == 'N':
        outputN = open(config['OUTPUT_FILES']['OUTPUT_FOLDER'] + config['OUTPUT_FILES']['NORMAL_TYPE_OUTPUT'], "a")
        appendGeneOutputInfo(name, gene.getOutputInfo(), outputN, countMatrix)
        outputN.close()
    elif tag == 'F':
        outputF = open(config['OUTPUT_FILES']['OUTPUT_FOLDER'] + config['OUTPUT_FILES']['FRAMESHIFT_TYPE_OUTPUT'], "a")
        appendGeneOutputInfo(name, gene.getOutputInfo(), outputF, countMatrix)
        outputF.close()
    elif tag == 'H':
        outputH = open(config['OUTPUT_FILES']['OUTPUT_FOLDER'] + config['OUTPUT_FILES']['HYPERSUSCEPTIBLE_TYPE_OUTPUT'], "a")
        appendGeneOutputInfo(name, gene.getOutputInfo(), outputH, countMatrix)
        outputH.close()
    elif tag == 'S':
        outputS = open(config['OUTPUT_FILES']['OUTPUT_FOLDER'] + config['OUTPUT_FILES']['SUPPRESSIBLE_TYPE_OUTPUT'], "a")
        appendGeneOutputInfo(name, gene.getOutputInfo(), outputS, countMatrix)
        outputS.close()
    else:
        outputI = open(config['OUTPUT_FILES']['OUTPUT_FOLDER'] + config['OUTPUT_FILES']['INTRINSIC_TYPE_OUTPUT'], "a")
        appendGeneOutputInfo(name, gene.getOutputInfo(), outputI, countMatrix)
        outputI.close()
    if config.getboolean('SETTINGS', 'DETAILED') and (gene.getOutputInfo()[0] > 0):
        detailedOutput = open(config['OUTPUT_FILES']['OUTPUT_FOLDER'] + config['OUTPUT_FILES']['DETAILED_FOLDER'] + "/" + name + ".csv", "w")
        gene.writeAdditionalInfo(detailedOutput)
        detailedOutput.close()
    
    gene.clearOutputInfo()
if type(countMatrix) == pd.DataFrame:
    countMatrix.to_csv(config['OUTPUT_FILES']['COUNT_MATRIX_FINAL'], index=False)
sys.exit(0)