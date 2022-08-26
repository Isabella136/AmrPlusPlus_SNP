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
import pysam, sys, getopt, os

snpInfoPrint = lambda a, b, c : a + ": " + b + " resitant reads out of " + c + " total reads\n"

inputFile = []
outputFolder = ""
argList = []
try:
    options, args = getopt.getopt(sys.argv[1:], "hlci:o:", ["mt_and_wt=", "detailed_output="])
except getopt.GetoptError:
    print("SNP_Verification.py -i <inputFile> -o <outputFolder>")
    sys.exit(-1)
for opt, arg in options:
    if opt == "-h":
        print("List of arguments:\n\n\n-c: conditions for redistribution\n-h: help\n-i: input file\n-o: output folder\n-l: license disclaimer\n\n--mt_and_wt: true by default, used in case of insertion leading to presence of both mt and wt; if true, mark as resistant; if false, mark as susceptible\n--detailed_output: false by default, determines whether a more detailed output will be given\n\n")
        sys.exit()
    elif opt == "-c":
        print("\n\nAMRPlusPlus_SNP_Verification\nCopyright (C) 2022  Nathalie Bonin\n\nThis program is free software: you can redistribute it and/or modify\nit under the terms of the GNU General Public License as published by\nthe Free Software Foundation, either version 3 of the License, or\n(at your option) any later version.\n\n")
        sys.exit()
    elif opt == "-l":
        print("\n\nAMRPlusPlus_SNP_Verification\nCopyright (C) 2022  Nathalie Bonin\n\nThis program is distributed in the hope that it will be useful,\nbut WITHOUT ANY WARRANTY; without even the implied warranty of\nMERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\nGNU General Public License for more details.\n\n")
        sys.exit()
    elif opt == "-i":
        inputFile.append(arg)
    elif opt == "-o":
        outputFolder = arg
    elif opt == "--mt_and_wt":
        if arg == "False":
            SNP_Verification_Tools.mt_and_wt = False
    elif opt =="--detailed_output":
        if arg != "None":
            SNP_Verification_Tools.detailed = True
            if arg != "All":
                argList = arg.split(',')
if len(inputFile) == 0:
    inputFile.append("Test/Frameshift.sam")
    inputFile.append("Test/Deletion4.sam")
    outputFolder = "Sample_Output"
    # print("SNP_Verification.py -i <inputFile> -o <outputFolder>")
    # sys.exit(-1)
if not(os.path.exists(outputFolder)):
    os.mkdir(outputFolder)

for file in inputFile:
    if not(os.path.exists(outputFolder + "/" + file.split('/')[-1][:-4])):
        os.mkdir(outputFolder + "/" + file.split('/')[-1][:-4])
    if SNP_Verification_Tools.detailed and not(os.path.exists(outputFolder + "/" + file.split('/')[-1][:-4] + "/Detailed_Output")):
        os.mkdir(outputFolder + "/" + file.split('/')[-1][:-4] + "/Detailed_Output")

SNPinfo = open("extracted_SNP_files/SNPinfo.fasta", "rt")
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

for filename in inputFile:
    #Analyze SAM file
    fullOutputPath = outputFolder + "/" + filename.split('/')[-1][:-4]
    pysam.sort("-o", fullOutputPath + "/Sorted_" + filename[filename.rfind("/")+1:], filename)
    samfile = pysam.AlignmentFile(fullOutputPath + "/Sorted_" + filename[filename.rfind("/")+1:], "r")
    iter = samfile.fetch()
    for read in iter:
        gene = geneDict.get(read.reference_name, False)
        if (gene == False):
            continue
        elif (read.cigarstring == None):
            continue
        elif (len(argList) != 0) and (gene.split("|")[0] not in argList):
            continue
        verify(read, gene)
        gene.resetForNextRead()
    samfile.close() 

    #Function that appends gene.getOutputInfo() to output file
    def appendGeneOutputInfo(name, outputInfo, file):
        file.write("\n" + name)
        for info in outputInfo.values():
            file.write("," + str(info))


    #Create output files and write headers
    outputN = open(fullOutputPath + "/Normal_Type_Genes_.csv", "w")
    outputN.write("Gene,Number of reads,Resistant,Missense,Insertion,Deletion,Previously recorded nonsense,N-tuple,Nonstop,12+bp indel,12+ bp frameshift,Newly found nonsense,Frameshift till end")
    outputN.close()
    outputF = open(fullOutputPath + "/Frameshift_Type_Genes.csv", "w")
    outputF.write("Gene,Number of reads,Resistant,Missense,Insertion,Deletion,Previously recorded nonsense,N-tuple,Nonstop,12+bp indel,12+ bp frameshift,Newly found nonsense,Frameshift till end")
    outputF.close()
    outputH = open(fullOutputPath + "/Hypersusceptible_Mutations_Type_Genes.csv", "w")
    outputH.write("Gene,Number of reads,Resistant,Missense,Insertion,Deletion,Previously recorded nonsense,N-tuple,Nonstop,12+bp indel,12+ bp frameshift,Newly found nonsense,Frameshift till end,Hypersusceptible mutations + resistance-conferring mutations")
    outputH.close()
    outputS = open(fullOutputPath + "/Suppressible_Frameshift_Type_Genes.csv", "w")
    outputS.write("Gene,Number of reads,Resistant,Missense,Insertion,Deletion,Previously recorded nonsense,N-tuple,Nonstop,12+bp indel,12+ bp frameshift,Newly found nonsense,Frameshift at end ,Suppressible frameshift at res 531,Frameshift at res 531 that is not suppressible")
    outputS.close()
    outputI = open(fullOutputPath + "/Intrinsic_Resistance_Genes.csv", "w")
    outputI.write("Gene,Number of reads,All,Some,None,Mutations,Acquired,12+bp indel,12+ bp frameshift,Nonsense,Frameshift till end")
    outputI.close()

    for name, gene in geneDict.items():
        if (len(argList) != 0) and (gene.split("|")[0] not in argList):
            continue
        tag = gene.getGeneTag()
        if tag == 'N':
            outputN = open(fullOutputPath + "/Normal_Type_Genes_.csv", "a")
            appendGeneOutputInfo(name, gene.getOutputInfo(), outputN)
            outputN.close()
        elif tag == 'F':
            outputF = open(fullOutputPath + "/Frameshift_Type_Genes.csv", "a")
            appendGeneOutputInfo(name, gene.getOutputInfo(), outputF)
            outputF.close()
        elif tag == 'H':
            outputH = open(fullOutputPath + "/Hypersusceptible_Mutations_Type_Genes.csv", "a")
            appendGeneOutputInfo(name, gene.getOutputInfo(), outputH)
            outputH.close()
        elif tag == 'S':
            outputS = open(fullOutputPath + "/Suppressible_Frameshift_Type_Genes.csv", "a")
            appendGeneOutputInfo(name, gene.getOutputInfo(), outputS)
            outputS.close()
        else:
            outputI = open(fullOutputPath + "/Intrinsic_Resistance_Genes.csv", "a")
            appendGeneOutputInfo(name, gene.getOutputInfo(), outputI)
            outputI.close()
        if SNP_Verification_Tools.detailed and (gene.getOutputInfo()[0] > 0):
            detailedOutput = open(fullOutputPath + "/Detailed_Output/" + name + ".csv", "w")
            gene.writeAdditionalInfo(detailedOutput)
            detailedOutput.close()
        
        gene.clearOutputInfo()

sys.exit(0)