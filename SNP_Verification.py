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

from SNP_Verification.Gene import Gene
from SNP_Verification.SNP import SNP
from SNP_Verification import dnaTranslate
import pysam
import sys, getopt

mt_and_wt = True #used in case of insertion leading to presence of both mt and wt; if true, mark as resistant; if false, mark as susceptible

geneDict = {

}

argInfoDict = {

}

intrinsicArgInfoDict = {}

intrinsic = ["MEG_1628", "MEG_3979", "MEG_3983", "MEG_4279", "MEG_4280", "MEG_4281", "MEG_4282", "MEG_6092", "MEG_6093"]

snpInfoPrint = lambda a, b, c : a + ": " + b + " resitant reads out of " + c + " total reads\n"

def intrinsicResistant(name, queryName, missing, messageType, aaOrNu):
    if name not in intrinsicArgInfoDict:
        intrinsicArgInfoDict.update({name:list()})
    if missing == None:
        if messageType == None:
            intrinsicArgInfoDict[name].append("Based on the sequence given, " + queryName + " contains the " + aaOrNu + " required for intrinsic resistance")
        elif messageType == "All":
            intrinsicArgInfoDict[name].append(queryName + " contains ALL " + aaOrNu + " required for intrinsic resistance")
        else:
            intrinsicArgInfoDict[name].append("The sequence given for " + queryName + " does not include the position where the " + aaOrNu + " required for intrinsic resistance are located")
    elif messageType == "MEG_6093":
        intrinsicArgInfoDict[name].append("The sequence given for " + queryName + " does not have the following amino acids required for resistance: " + str(missing[0]) + "; nor does it have this alternative group of amino acids that can also induce resistance: " + str(missing[1]))
    else:
        intrinsicArgInfoDict[name].append("The sequence given for " + queryName + " does not have the following " + aaOrNu + " required for intrinsic resistance:" + str(missing))
def resistant(name, increment):
    argInfo = argInfoDict.get(name, False)
    if (argInfo == False):
        argInfoDict.update({name:(increment, 1)})
    else:
        temp = list(argInfo)
        temp[0] += increment
        temp[1] += 1
        argInfo = tuple(temp)
        argInfoDict.update({name:argInfo})

def disregard(name):
    resistant(name, 0)

def extendCigar(cigar):
    toReturn = ""
    count = 0
    for c in cigar:
        if c.isdigit():
            count = count * 10 + int(c)
        else:
            while count > 0:
                toReturn = toReturn + c
                count -= 1
    return toReturn

def mapCigarToAlignment(cigar, aligned_pair, rRna):
    alignment_map = []
    queryLength = 0
    refLength = 0
    index = 0
    shift = 0
    prevShift = 0
    splitIndex = -1
    translatable = False
    for op in cigar:
        if op == "S":
            index += 1
        elif (op == "M") | (op == "X") | (op == "="):
            prevShift = shift
            if not(translatable) and not(rRna):
                if splitIndex < 0:
                    splitIndex = aligned_pair[index][1] % 3
                else:
                    splitIndex += 1
                if (splitIndex % 3) != 0:
                    index += 1
                else:
                    translatable = True
                    queryLength += 1
                    refLength += 1
                    alignment_map.append(("M", aligned_pair[index][0], aligned_pair[index][1], shift, prevShift))
                    index += 1
            else:
                if not(translatable):
                    translatable = True
                queryLength += 1
                refLength += 1
                alignment_map.append(("M", aligned_pair[index][0], aligned_pair[index][1], shift, prevShift))
                index += 1
        elif op == "I":
            prevShift = shift
            shift += 1
            if not(translatable) and not(rRna):
                if splitIndex < 0:
                    temp = aligned_pair[index][1]
                    i = 1
                    while temp == None:
                        temp = aligned_pair[index+i][1]
                        i += 1
                    splitIndex = temp % 3
                else:
                    splitIndex += 1
                if (splitIndex % 3) != 0:
                    index += 1
                else:
                    translatable = True
                    queryLength += 1
                    alignment_map.append((op, aligned_pair[index][0], aligned_pair[index][1], shift, prevShift))
                    index += 1
            else:
                if not(translatable):
                    translatable = True
                queryLength += 1
                alignment_map.append((op, aligned_pair[index][0], aligned_pair[index][1], shift, prevShift))
                index += 1
        elif op == "D":
            prevShift = shift
            shift -= 1
            if not(translatable) and not(rRna):
                if ((splitIndex % 3) != 2) and (splitIndex != -1):
                    index += 1
                else:
                    translatable = True
                    refLength += 1
                    alignment_map.append((op, aligned_pair[index][0], aligned_pair[index][1], shift, prevShift))
                    index += 1
            else:
                if not(translatable):
                    translatable = True
                refLength += 1
                alignment_map.append((op, aligned_pair[index][0], aligned_pair[index][1], shift, prevShift))
                index += 1
        else:
            continue
    if not(rRna):
        while (queryLength % 3) != 0:
            poped = alignment_map.pop()
            if poped[0] != "D" : queryLength -= 1
    return alignment_map

def transformNtAlignmentMap(nt_alignment_map):
    new_nt_alignment_map = {}
    ntQueryIndex = 0
    for nt in nt_alignment_map:
        if nt[0] == "I":
            continue
        new_nt_alignment_map.update({nt[2]:(ntQueryIndex,)})
        ntQueryIndex += 1
    return new_nt_alignment_map

def aaAlignment(nt_alignment_map):
    aa_alignment_map = {

    }
    ntRefIndex = nt_alignment_map[0][2]
    i = 1
    while ntRefIndex == None:
        ntRefIndex = nt_alignment_map[i][2]
        i+=1
    ntQueryIndex = 0
    aaQueryIndex = 0
    insertCount = 0
    prevAaShift = None    
    firstAlignment = None
    thirdAlignment = 0
    mapIndex = 0
    inbetween = None
    deleteCount = 0
    hasDeletion = False
    def addToMapDeletion(index, toAdd):
        aa = aa_alignment_map.get(index-1, False)
        if aa == False:
            aa_alignment_map.update({index-1:(inbetween, )})
        else:
            temp = list(aa)
            temp.append(inbetween)
            aa = tuple(temp)
            aa_alignment_map.update({index-1:aa})
        addToMap(index, toAdd)
    def addToMapInbetween(index, toAdd):
        aa = aa_alignment_map.get(index, False)
        if aa == False:
            if (inbetween != None):
                aa_alignment_map.update({index:(inbetween, toAdd, )})
            else:
                aa_alignment_map.update({index:(toAdd,)})
        else:
            temp = list(aa)
            if (inbetween != None):
                temp.append(inbetween)
            temp.append(toAdd)
            aa = tuple(temp)
            aa_alignment_map.update({index:aa})
    def addToMap(index, toAdd):
        aa = aa_alignment_map.get(index, False)
        if aa == False:
            aa_alignment_map.update({index:(toAdd,)})
        else:
            temp = list(aa)
            temp.append(toAdd)
            aa = tuple(temp)
            aa_alignment_map.update({index:aa})
    for nt in nt_alignment_map:
        if (ntQueryIndex % 3) == 0:
            if deleteCount == 0: #1/2D3M
                hasDeletion = False
        if (ntRefIndex % 3) == 0:
            if nt[0] == "I":
                insertCount += 1
                deleteCount = 0
                if (nt[3] % 3) == 0:  
                    if insertCount == 3: #3I
                        if inbetween == None:
                            addToMap(int(ntRefIndex/3)-1, thirdAlignment)
                        inbetween = thirdAlignment
                    else: #1I...1M2I or 2I...2M1I
                        addToMapInbetween(int(ntRefIndex/3)-1, thirdAlignment)
                        inbetween = None
                    insertCount = 0
                if prevAaShift == None:
                    prevAaShift = nt[4]
                ntQueryIndex += 1
            elif nt[0] == "D":
                ntRefIndex += 1
                if prevAaShift == None:
                    prevAaShift = nt[4]
                insertCount = 0
                deleteCount += 1
                if deleteCount == 3:
                    addToMap(int(ntRefIndex/3)-1, '-')
                    inbetween = '-'
                hasDeletion = True
            else: #nt[0] == "M"
                deleteCount = 0
                ntQueryIndex += 1
                ntRefIndex += 1
                insertCount = 0
                if prevAaShift == None:
                    prevAaShift = nt[4]
                if (ntQueryIndex % 3) == 0: #2I1M3M
                    if (prevAaShift == 0):
                        inbetween = None
                    elif inbetween == None:
                        if not(hasDeletion):
                            addToMap(int(ntRefIndex/3)-1, thirdAlignment)
                        if nt[3] < 0:
                            addToMap(int(ntRefIndex/3), thirdAlignment)
                        inbetween = thirdAlignment
                    else:
                        inbetween = thirdAlignment
                firstAlignment = aaQueryIndex
        elif (ntRefIndex % 3) == 1:
            if nt[0] == "I":
                deleteCount = 0
                ntQueryIndex += 1
                if (ntQueryIndex % 3) == 0: #1I2M2M1I
                    if prevAaShift == 0:
                        inbetween = None
                    elif inbetween == None:
                        addToMap(int(ntRefIndex/3)-1, thirdAlignment)
                        inbetween = thirdAlignment
                    else:
                        inbetween = thirdAlignment
            elif nt[0] == "D":
                ntRefIndex += 1
                deleteCount += 1
                if deleteCount == 3:
                    addToMap(int(ntRefIndex/3)-1, '-')
                    inbetween = '-'
                hasDeletion = True
            else: #nt[0] == "M"
                deleteCount = 0
                ntQueryIndex += 1
                ntRefIndex += 1
                if firstAlignment == None:
                    firstAlignment = aaQueryIndex
                if prevAaShift == None:
                    prevAaShift = nt[4]
                if (ntQueryIndex % 3) == 0: #1I2M3M
                    if (prevAaShift == 0):
                        inbetween = None
                    elif inbetween == None:
                        if not(hasDeletion):
                            addToMap(int(ntRefIndex/3)-1, thirdAlignment)
                        if nt[3] < 0:
                            addToMap(int(ntRefIndex/3), thirdAlignment)
                        inbetween = thirdAlignment
                    else:
                        inbetween = thirdAlignment
        else: #(ntRefIndex % 3) == 2
            if nt[0] == "I":
                deleteCount = 0
                ntQueryIndex += 1
            elif nt[0] == "D":
                hasDeletion = True
                deleteCount += 1
                ntRefIndex += 1
                if firstAlignment != None:
                    if (prevAaShift % 3) == 0: #1M2D, 2M1D
                        addToMap(int(nt[2]/3), thirdAlignment)
                else: #3D
                    addToMap(int(nt[2]/3), '-')
                    if (prevAaShift % 3) != 0:
                        inbetween = '-'
                prevAaShift = None
                firstAlignment = None
            else: #nt[0] == "M"
                deleteCount = 0
                ntQueryIndex += 1
                ntRefIndex += 1
                if firstAlignment == None:
                    firstAlignment = aaQueryIndex
                if prevAaShift == None:
                    prevAaShift = nt[4]
                if firstAlignment == thirdAlignment:
                    if ((prevAaShift % 3) == 0) | ((nt[4] % 3) == 0): #3M, 1D1I2M, 1D/2I...2D1M, 2D/1I...1D2M, 2D3M, 1D3M
                        if inbetween != None:
                            if hasDeletion:
                                if inbetween != '-':
                                    inbetween = thirdAlignment
                                    addToMapDeletion(int(nt[2]/3), thirdAlignment)
                                else:
                                    addToMapInbetween(int(nt[2]/3), thirdAlignment)
                            else:
                                addToMapInbetween(int(nt[2]/3), thirdAlignment)
                            inbetween = None
                        else:
                            addToMap(int(nt[2]/3), thirdAlignment)
                elif (prevAaShift % 3) == 0:
                    if (nt[4] % 3) == 0: #3Mand3xI
                        if inbetween != None:
                            aa_alignment_map.update({int(nt[2]/3):(inbetween, firstAlignment, thirdAlignment)})
                            inbetween = None
                        else:
                            aa_alignment_map.update({int(nt[2]/3):(firstAlignment, thirdAlignment)})
                    else: #1/2Iand3M
                        addToMapInbetween(int(nt[2]/3), firstAlignment)
                        inbetween = None
                elif (nt[4] % 3) == 0: #2D/1I...1M1D1M
                    if inbetween != None:
                        if hasDeletion:
                            if inbetween != '-':
                                inbetween = thirdAlignment
                                addToMapDeletion(int(nt[2]/3), thirdAlignment)
                            else:
                                addToMapInbetween(int(nt[2]/3), thirdAlignment)
                        else:
                            addToMapInbetween(int(nt[2]/3), thirdAlignment)
                        inbetween = None
                    else:
                        addToMap(int(nt[2]/3), thirdAlignment)
                prevAaShift = None
                firstAlignment = None

        aaQueryIndex = int(ntQueryIndex / 3)
        thirdAlignment = aaQueryIndex
        mapIndex += 1
    return aa_alignment_map    

def verifyNonsense(stopLocation, aa_alignment_map, gene, name):
    SNPInfo = gene.condensedNonInfo()
    if (len(SNPInfo) == 0):
        disregard(name)
    else:
        res = 0
        for snp in SNPInfo:
            stopLocTemp = aa_alignment_map.get(snp[1]-1,False)
            if stopLocTemp != False:
                if stopLocTemp[0] == stopLocation:
                    res = 1
                    break
        resistant(name, res)

def verifyMultiple(aa_alignment_map, gene, name, aaQuerySequence):
    SNPInfo = gene.condensedMultInfo()
    if (len(SNPInfo) == 0):
        disregard(name)
    else:
        resBool = True
        for snpMult in SNPInfo:
            for snp in snpMult:
                #codons that count for one deletion mutation but aren't fully deleted in read
                #if >1, snp is disregarded
                codonNotFullyDeleted = 0 
                if (aa_alignment_map.get(snp[1]-1,False)) == False:
                    resBool = False
                    break
                for mt in snp[2]: 
                    delMut = 0
                    misMut = 0
                    actualResBool = False
                    for queryIndex in tuple(aa_alignment_map[snp[1]-1]):
                        if queryIndex == '-':
                            if mt == queryIndex:
                                delMut += 1
                                actualResBool = True
                                continue
                        elif mt == aaQuerySequence[queryIndex]:
                            actualResBool = True
                            break
                        else:
                            misMut += 1
                            resBool = False
                    if (delMut > 0) & (misMut > 0):
                        codonNotFullyDeleted += 1
                    if actualResBool: 
                        resBool = actualResBool
                        break
                if not(resBool): break
            if resBool: 
                if codonNotFullyDeleted < 2:
                    resistant(name, 1)
                    break
                else:
                    resBool = False
        if not(resBool):
            disregard(name)

def verify(read, gene):
    name = gene.getName()
    rRna = gene.rRna()
    if not(rRna):
        cigarOpCount = read.get_cigar_stats()[0].tolist()
        #if (insertions - deletions) %3 != 0, disregard
        if ((cigarOpCount[1] - cigarOpCount[2]) % 3) != 0: 
            disregard(name) 
            return None
    cigar = extendCigar(read.cigarstring)
    aligned_pair = read.get_aligned_pairs()
    nt_alignment_map = mapCigarToAlignment(cigar, aligned_pair, rRna)
    querySeq = read.query_sequence
    mapOfInterest = transformNtAlignmentMap(nt_alignment_map)
    start = nt_alignment_map[0][1]
    i = 1
    while start == None:
        start = nt_alignment_map[i][1]
        i += 1
    end = nt_alignment_map[len(nt_alignment_map)-1][1]
    i = 2
    while end == None:
        end = nt_alignment_map[len(nt_alignment_map)-i][1]
        i += 1
    trimmedQuerySequence = querySeq[start:end+1]
    seqOfInterest = trimmedQuerySequence
    if not(rRna):
        aaQuerySequence = dnaTranslate(trimmedQuerySequence)
        stopLocation = aaQuerySequence.find('*') 
        aa_alignment_map = aaAlignment(nt_alignment_map)
        #if nonsense mutations
        if (stopLocation != -1) & ((stopLocation < len(aaQuerySequence)-1) | (read.reference_end < gene.ntSeqLength())): 
            verifyNonsense(stopLocation, aa_alignment_map, gene, name)
        seqOfInterest = aaQuerySequence
        mapOfInterest = aa_alignment_map
    else :
        begin = mapOfInterest[0][0]+1
        end = mapOfInterest[list(mapOfInterest.keys())[-1]][0]+1
        mustList, all = gene.getFirstMustBetweenParams(begin, end)
        missingList = []
        if mustList != None:
            for must in mustList:
                listInMissingList = []
                hasAllMust = True
                while(must.getPos() < end):
                    hasMust = False
                    for queryIndex in mapOfInterest[must.getPos()-1]:
                        if queryIndex == None:
                            continue
                        if must.getWt() == seqOfInterest[queryIndex]:
                            hasMust = True
                            continue
                    if not(hasMust): 
                        hasAllMust = False
                        listInMissingList.append(must.condensedInfo())
                    must = must.getNext()
                    if must == None:
                        all = all & True
                        break
                if hasAllMust:
                    missingList = None
                    break
                elif name == "MEG_6093":
                    missingList.append(listInMissingList)
                else:
                    missingList = listInMissingList
            # def intrinsicResistant(name, queryName, missing, messageType, aaOrNu):
            messageType = None
            if name == "MEG_6093": messageType = "MEG_6093"
            if all: messageType = "All"
            intrinsicResistant(name, read.query_name, missingList, messageType, gene.aaOrNu())
            if missingList == None:
                resistant(name, 1)
        elif name in intrinsic:
            intrinsicResistant(name, read.query_name, None, "NA", gene.aaOrNu())
        SNPInfo = gene.condensedRegDelInfo()
        res = 0
        for snp in SNPInfo:
            if (mapOfInterest.get(snp[1]-1,False)) == False:
                continue
            for mt in snp[2]: 
                for queryIndex in tuple(mapOfInterest[snp[1]-1]):
                    if queryIndex == '-': continue
                    if mt_and_wt:
                        if mt == seqOfInterest[queryIndex]:
                            res = 1
                            break
                    else:
                        if mt == seqOfInterest[queryIndex]:
                            res = 1
                        elif snp[0] == seqOfInterest[queryIndex]:
                            res = -1
                            break
                if res == 1: break
                elif res == -1:
                    res = 0
                    break
            if res == 1: break
        if res == 0 : #take care of Mult
            verifyMultiple(mapOfInterest, gene, name, seqOfInterest)
        else: resistant(name, res)

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
    output = open(outputFolder + "/" + name[name.rfind("/")+1:] + ".txt", "w")
    pysam.sort("-o", outputFolder + "/Sorted_" + name[name.rfind("/")+1:] + ".sam", name)
    samfile = pysam.AlignmentFile(outputFolder + "/Sorted_" + name[name.rfind("/")+1:], "r")
    iter = samfile.fetch()
    for read in iter:
        gene = geneDict.get(read.reference_name, False)
        if (gene == False):
            continue
        elif (read.cigarstring == None):
            continue
        verify(read, gene)
    samfile.close() 
    for name in argInfoDict:
        output.write(snpInfoPrint(name, str(argInfoDict[name][0]), str(argInfoDict[name][1])))
    argInfoDict = {
    }
    output.close()
sys.exit(0)