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
from SNP_Verification_Tools.SNP import SNP
from SNP_Verification_Processes import verify
import pysam
import sys, getopt

mt_and_wt = True #used in case of insertion leading to presence of both mt and wt; if true, mark as resistant; if false, mark as susceptible

geneDict = {}
argInfoDict = {}
intrinsicInfoDict = {}
frameshiftInfoDict = {}
meg_3180InfoDict = {}
snpInfoPrint = lambda a, b, c : a + ": " + b + " resitant reads out of " + c + " total reads\n"

inputFile = []
outputFolder = ""
try:
    options, args = getopt.getopt(sys.argv[1:], "hwci:o:", ["mt_and_wt="])
except getopt.GetoptError:
    print("SNP_Verification.py -i <inputFile> -o <outputFolder>")
    sys.exit(-1)
for opt, arg in options:
    if opt == "-h":
        print("List of arguments:\n\n\n-c: conditions for redistribution\n-h: help\n-i: input file\n-o: output folder\n-w: warranty disclaimer\n\n--mt_and_wt: true by default, used in case of insertion leading to presence of both mt and wt; if true, mark as resistant; if false, mark as susceptible\n\n")
        sys.exit()
    elif opt == "-c":
        print("\n\nAMRPlusPlus_SNP_Verification\nCopyright (C) 2022  Nathalie Bonin\n\nThis program is free software: you can redistribute it and/or modify\nit under the terms of the GNU General Public License as published by\nthe Free Software Foundation, either version 3 of the License, or\n(at your option) any later version.\n\n")
        sys.exit()
    elif opt == "-w":
        print("\n\nAMRPlusPlus_SNP_Verification\nCopyright (C) 2022  Nathalie Bonin\n\nThis program is distributed in the hope that it will be useful,\nbut WITHOUT ANY WARRANTY; without even the implied warranty of\nMERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\nGNU General Public License for more details.\n\n")
        sys.exit()
    elif opt == "-i":
        inputFile.append(arg)
    elif opt == "-o":
        outputFolder = arg
    elif opt == "--mt_and_wt":
        if arg == "False":
            mt_and_wt = False
if len(inputFile) == 0:
    print("SNP_Verification.py -i <inputFile> -o <outputFolder>")
    sys.exit(-1)

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

for name in inputFile:

    #Analyze SAM file
    pysam.sort("-o", outputFolder + "/Sorted_" + name[name.rfind("/")+1:], name)
    samfile = pysam.AlignmentFile(outputFolder + "/Sorted_" + name[name.rfind("/")+1:], "r")
    iter = samfile.fetch()
    for read in iter:
        gene = geneDict.get(read.reference_name, False)
        if (gene == False):
            continue
        elif (read.cigarstring == None):
            continue
        verify(read, gene, argInfoDict, intrinsicInfoDict, frameshiftInfoDict, meg_3180InfoDict, mt_and_wt)
    samfile.close() 

    #Output SNP Info
    output = open(outputFolder + "/" + name[name.rfind("/")+1:] + ".txt", "w")
    for argName in argInfoDict:
        output.write(snpInfoPrint(argName, str(argInfoDict[argName][0]), str(argInfoDict[argName][1])))
    argInfoDict = {
    }
    output.close()

    #Output Intrinsic Resistance Info
    if len(intrinsicInfoDict) != 0:
        output = open(outputFolder + "/intrinsic_resistance_" + name[name.rfind("/")+1:] + ".csv", "w")
        output.write("ARG,All nuc/aa required present,Some nuc/aa positions in query,No nuc/aa positions in query,Mutants in nuc/aa positions\n")
        for argName in intrinsicInfoDict:
            output.write(argName + ",")
            if "All" in intrinsicInfoDict[argName]:
                for query in intrinsicInfoDict[argName]["All"][:-1]:
                    output.write(query + ";")
                output.write(intrinsicInfoDict[argName]["All"][-1] + ",")
            else: output.write(",")
            if "Some" in intrinsicInfoDict[argName]:
                for query in intrinsicInfoDict[argName]["Some"][:-1]:
                    output.write(query + ";")
                output.write(intrinsicInfoDict[argName]["Some"][-1] + ",")
            else: output.write(",")
            if "NA" in intrinsicInfoDict[argName]:
                for query in intrinsicInfoDict[argName]["NA"][:-1]:
                    output.write(query + ";")
                output.write(intrinsicInfoDict[argName]["NA"][-1] + ",")
            else: output.write(",")
            if "Mutant" in intrinsicInfoDict[argName]:
                for query in intrinsicInfoDict[argName]["Mutant"][:-1]:
                    output.write(query + ";")
                output.write(intrinsicInfoDict[argName]["Mutant"][-1] + "\n")
            else: output.write("\n")
        intrinsicInfoDict = {}
        output.close()

    #Output Info on Insertion/Deletion Difference
    if len(frameshiftInfoDict) != 0:
        output = open(outputFolder + "/insertion_deletion_frameshift_" + name[name.rfind("/")+1:] + ".csv", "w")
        output.write("ARG,Reads with a frameshift at the end of their sequence\n")
        for argName in frameshiftInfoDict:
            output.write(argName + ",")
            for query in frameshiftInfoDict[argName][:-1]:
                output.write(query + ",")
            output.write(frameshiftInfoDict[argName][-1] + "\n")
        frameshiftInfoDict = {}
        output.close()

sys.exit(0)