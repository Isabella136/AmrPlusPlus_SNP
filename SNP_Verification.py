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
import pysam, sys, getopt, os, configparser

argList = []
configFile = "config.ini"
config = configparser.ConfigParser()
inputFile = ""
outputFolder = ""
mt_and_wt = ""
detailed_output = ""

try:
    options, args = getopt.getopt(sys.argv[1:], "hlrci:o:", ["mt_and_wt=", "detailed_output="])
except getopt.GetoptError:
    print("SNP_Verification.py -i <inputFile> -o <outputFolder>")
    sys.exit(-1)
for opt, arg in options:
    if opt == "-h":
        print("List of arguments:\n\n\n-c:config.file\n-h: help\n-i: input file\n-l: license disclaimer\n-o: output folder\n-r: conditions for redistribution\n\n--mt_and_wt: true by default, used in case of insertion leading to presence of both mt and wt; if true, mark as resistant; if false, mark as susceptible\n--detailed_output: false by default, determines whether a more detailed output will be given\n\n")
        sys.exit()
    elif opt == "-r":
        print("\n\nAMRPlusPlus_SNP_Verification\nCopyright (C) 2022  Nathalie Bonin\n\nThis program is free software: you can redistribute it and/or modify\nit under the terms of the GNU General Public License as published by\nthe Free Software Foundation, either version 3 of the License, or\n(at your option) any later version.\n\n")
        sys.exit()
    elif opt == "-l":
        print("\n\nAMRPlusPlus_SNP_Verification\nCopyright (C) 2022  Nathalie Bonin\n\nThis program is distributed in the hope that it will be useful,\nbut WITHOUT ANY WARRANTY; without even the implied warranty of\nMERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\nGNU General Public License for more details.\n\n")
        sys.exit()
    elif opt == "-c":
        configFile = arg
    elif opt == "-i":
        inputFile = arg
    elif opt == "-o":
        outputFolder = arg
    elif opt == "--mt_and_wt":
        mt_and_wt = arg
    elif opt =="--detailed_output":
        if arg != "false":
            detailed_output = "true"
            if arg != "All":
                argList = arg.split(',')

config.read(configFile)
if len(inputFile) > 0:
    config['SOURCE_FILES']['SAM_INPUT'] = inputFile
if len(outputFolder) > 0:
    config['OUTPUT_FILES']['OUTPUT_FOLDER'] = outputFolder + "/" + config['SOURCE_FILES']['SAM_INPUT'].split('/')[-1][:-4]
if len(mt_and_wt) > 0:
    config['SETTINGS']['MT_AND_WT'] = mt_and_wt
if len(detailed_output) > 0:
    config['SETTINGS']['DETAILED'] = detailed_output

with open(configFile, 'w') as configFile:
    config.write(configFile)

if not(os.path.exists(config['OUTPUT_FILES']['OUTPUT_FOLDER'])):
    os.mkdir(config['OUTPUT_FILES']['OUTPUT_FOLDER'])
if config.getboolean('SETTINGS', 'DETAILED') and not(os.path.exists(config['OUTPUT_FILES']['OUTPUT_FOLDER'] + config['OUTPUT_FILES']['DETAILED_FOLDER'])):
    os.mkdir(config['OUTPUT_FILES']['OUTPUT_FOLDER'] + config['OUTPUT_FILES']['DETAILED_FOLDER'])

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
    elif (len(argList) != 0) and (gene.split("|")[0] not in argList):
        continue
    verify(read, gene, config)
    gene.resetForNextRead()
samfile.close() 

#Function that appends gene.getOutputInfo() to output file
def appendGeneOutputInfo(name, outputInfo, file):
    file.write("\n" + name)
    for info in outputInfo.values():
        file.write("," + str(info))


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

for name, gene in geneDict.items():
    if (len(argList) != 0) and (gene.split("|")[0] not in argList):
        continue
    tag = gene.getGeneTag()
    if tag == 'N':
        outputN = open(config['OUTPUT_FILES']['OUTPUT_FOLDER'] + config['OUTPUT_FILES']['NORMAL_TYPE_OUTPUT'], "a")
        appendGeneOutputInfo(name, gene.getOutputInfo(), outputN)
        outputN.close()
    elif tag == 'F':
        outputF = open(config['OUTPUT_FILES']['OUTPUT_FOLDER'] + config['OUTPUT_FILES']['FRAMESHIFT_TYPE_OUTPUT'], "a")
        appendGeneOutputInfo(name, gene.getOutputInfo(), outputF)
        outputF.close()
    elif tag == 'H':
        outputH = open(config['OUTPUT_FILES']['OUTPUT_FOLDER'] + config['OUTPUT_FILES']['HYPERSUSCEPTIBLE_TYPE_OUTPUT'], "a")
        appendGeneOutputInfo(name, gene.getOutputInfo(), outputH)
        outputH.close()
    elif tag == 'S':
        outputS = open(config['OUTPUT_FILES']['OUTPUT_FOLDER'] + config['OUTPUT_FILES']['SUPPRESSIBLE_TYPE_OUTPUT'], "a")
        appendGeneOutputInfo(name, gene.getOutputInfo(), outputS)
        outputS.close()
    else:
        outputI = open(config['OUTPUT_FILES']['OUTPUT_FOLDER'] + config['OUTPUT_FILES']['INTRINSIC_TYPE_OUTPUT'], "a")
        appendGeneOutputInfo(name, gene.getOutputInfo(), outputI)
        outputI.close()
    if config.getboolean('SETTINGS', 'DETAILED') and (gene.getOutputInfo()[0] > 0):
        detailedOutput = open(config['OUTPUT_FILES']['OUTPUT_FOLDER'] + config['OUTPUT_FILES']['DETAILED_FOLDER'] + "/" + name + ".csv", "w")
        gene.writeAdditionalInfo(detailedOutput)
        detailedOutput.close()
    
    gene.clearOutputInfo()

sys.exit(0)