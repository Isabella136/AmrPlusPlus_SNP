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

mt_and_wt = True #used in case of insertion leading to presence of both mt and wt; if true, mark as resistant; if false, mark as susceptible

geneDict = {

}

argInfoDict = {

}

snpInfoPrint = lambda a, b, c : a + ": " + b + " resitant reads out of " + c + " total reads\n"

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

def mapCigarToAlignment(cigar, aligned_pair):
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
            if not(translatable):
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
                queryLength += 1
                refLength += 1
                alignment_map.append(("M", aligned_pair[index][0], aligned_pair[index][1], shift, prevShift))
                index += 1
        elif op == "I":
            prevShift = shift
            shift += 1
            if not(translatable):
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
                queryLength += 1
                alignment_map.append((op, aligned_pair[index][0], aligned_pair[index][1], shift, prevShift))
                index += 1
        elif op == "D":
            prevShift = shift
            shift -= 1
            if not(translatable):
                if ((splitIndex % 3) != 2) and (splitIndex != -1):
                    index += 1
                else:
                    translatable = True
                    refLength += 1
                    alignment_map.append((op, aligned_pair[index][0], aligned_pair[index][1], shift, prevShift))
                    index += 1
            else:
                refLength += 1
                alignment_map.append((op, aligned_pair[index][0], aligned_pair[index][1], shift, prevShift))
                index += 1
        else:
            continue
    while (queryLength % 3) != 0:
        poped = alignment_map.pop()
        if poped[0] != "D" : queryLength -= 1
    return alignment_map



def aaAlignment(nt_alignment_map):
    aa_alignment_map = {

    }
    def addToMapDeletion():
    def addToMapInbetween(index, toAdd):
        aa = aa_alignment_map.get(index, False)
        if aa == False:
            aa_alignment_map.update({index:(toAdd,)})
        else:
            temp = list(aa)
            temp.append(toAdd)
            aa = tuple(temp)
            aa_alignment_map.update({index:aa})
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
    for nt in nt_alignment_map:
        if (ntQueryIndex % 3) == 0:
            if deleteCount > 0: #1/2D3M
                hasDeletion = False
        if (ntRefIndex % 3) == 0:
            if nt[0] == "I":
                insertCount += 1
                deleteCount = 0
                if (nt[3] % 3) == 0:  
                    if insertCount == 3: #3I
                        if inbetween == None:
                            addToMapInbetween(int(ntRefIndex/3)-1, thirdAlignment)
                        inbetween = thirdAlignment
                    else: #1I...1M2I or 2I...2M1I
                        aa = aa_alignment_map.get(int(ntRefIndex/3)-1, False)
                        if aa == False:
                            if (inbetween != None):
                                aa_alignment_map.update({int(ntRefIndex/3)-1:(inbetween, thirdAlignment, )})
                            else:
                                aa_alignment_map.update({int(ntRefIndex/3)-1:(thirdAlignment,)})
                        else:
                            temp = list(aa)
                            if (inbetween != None):
                                temp.append(inbetween)
                            temp.append(thirdAlignment)
                            aa = tuple(temp)
                            aa_alignment_map.update({int(ntRefIndex/3)-1:aa})
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
                    addToMapInbetween(int(ntRefIndex/3)-1, '-')
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
                            addToMapInbetween(int(ntRefIndex/3)-1, thirdAlignment)
                        if nt[3] < 0:
                            addToMapInbetween(int(ntRefIndex/3), thirdAlignment)
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
                        addToMapInbetween(int(ntRefIndex/3)-1, thirdAlignment)
                        inbetween = thirdAlignment
                    else:
                        inbetween = thirdAlignment
            elif nt[0] == "D":
                ntRefIndex += 1
                deleteCount += 1
                if deleteCount == 3:
                    addToMapInbetween(int(ntRefIndex/3)-1, '-')
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
                            addToMapInbetween(int(ntRefIndex/3)-1, thirdAlignment)
                        if nt[3] < 0:
                            addToMapInbetween(int(ntRefIndex/3), thirdAlignment)
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
                        if inbetween != None:
                            aa_alignment_map.update({int(nt[2]/3):(inbetween, thirdAlignment,)})
                            inbetween = None
                        else:
                            aa_alignment_map.update({int(nt[2]/3):(thirdAlignment,)})
                else: #3D
                    if (prevAaShift % 3) == 0:
                        if inbetween != None:
                            aa_alignment_map.update({int(nt[2]/3):('-', inbetween, )})
                            inbetween = None
                        else:
                            aa_alignment_map.update({int(nt[2]/3):('-',)})
                    else:
                        aa_alignment_map.update({int(nt[2]/3):('-',)})
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
                        if (inbetween != None) & (inbetween != thirdAlignment):
                            if hasDeletion:
                                if inbetween != '-':
                                    inbetween = thirdAlignment
                                aa = aa_alignment_map.get(int(nt[2]/3) - 1, False)
                                if aa == False:
                                    aa_alignment_map.update({int(nt[2]/3) - 1:(inbetween,)})
                                else:
                                    temp = list(aa)
                                    temp.append(inbetween)
                                    aa = tuple(temp)
                                    aa_alignment_map.update({int(nt[2]/3) - 1:aa})
                                aa_alignment_map.update({int(nt[2]/3):(thirdAlignment, )})
                            else:
                                aa_alignment_map.update({int(nt[2]/3):(inbetween, thirdAlignment,)})
                            inbetween = None
                        else:
                            aa_alignment_map.update({int(nt[2]/3):(thirdAlignment,)})
                elif (prevAaShift % 3) == 0:
                    if (nt[4] % 3) == 0: #3Mand3xI
                        if inbetween != None:
                            aa_alignment_map.update({int(nt[2]/3):(inbetween, firstAlignment, thirdAlignment)})
                            inbetween = None
                        else:
                            aa_alignment_map.update({int(nt[2]/3):(firstAlignment, thirdAlignment)})
                    else: #1/2Iand3M
                        if inbetween != None:
                            aa_alignment_map.update({int(nt[2]/3):(inbetween, firstAlignment,)})
                            inbetween = None
                        else:
                            aa_alignment_map.update({int(nt[2]/3):(firstAlignment,)})
                elif (nt[4] % 3) == 0: #2D/1I...1M1D1M
                    if (inbetween != None) & (inbetween != thirdAlignment):
                        if hasDeletion:
                            if inbetween != '-':
                                inbetween = thirdAlignment
                            aa = aa_alignment_map.get(int(nt[2]/3) - 1, False)
                            if aa == False:
                                aa_alignment_map.update({int(nt[2]/3) - 1:(inbetween,)})
                            else:
                                temp = list(aa)
                                temp.append(inbetween)
                                aa = tuple(temp)
                                aa_alignment_map.update({int(nt[2]/3) - 1:aa})
                            aa_alignment_map.update({int(nt[2]/3):(thirdAlignment, )})
                        else:
                            aa_alignment_map.update({int(nt[2]/3):(inbetween, thirdAlignment,)})
                        inbetween = None
                    else:
                        aa_alignment_map.update({int(nt[2]/3):(thirdAlignment,)})  
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
                if (aa_alignment_map.get(snp[1]-1,False)) == False:
                    resBool = False
                    continue
                for mt in snp[2]: 
                    for queryIndex in tuple(aa_alignment_map[snp[1]-1]):
                        if queryIndex == '-':
                            if mt == queryIndex:
                                resBool = True
                                break
                        elif mt == aaQuerySequence[queryIndex]:
                            resBool = True
                            break
                        else:
                            resBool = False
                    if resBool: break
                if not(resBool): break
            if resBool: 
                resistant(name, 1)
                break
        if not(resBool):
            disregard(name)
def verify(read, gene):
    name = gene.getName()
    cigarOpCount = read.get_cigar_stats()[0].tolist()
    #if (insertions - deletions) %3 != 0, disregard
    if ((cigarOpCount[1] - cigarOpCount[2]) % 3) != 0: disregard(name) 
    else :
        cigar = extendCigar(read.cigarstring)
        aligned_pair = read.get_aligned_pairs()
        nt_alignment_map = mapCigarToAlignment(cigar, aligned_pair)
        querySeq = read.query_sequence
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
        aaQuerySequence = dnaTranslate(trimmedQuerySequence)
        stopLocation = aaQuerySequence.find('*') 
        aa_alignment_map = aaAlignment(nt_alignment_map)
        #if nonsense mutations
        if (stopLocation != -1) & ((stopLocation < len(aaQuerySequence)-1) | (read.reference_end < gene.ntSeqLength())): 
            verifyNonsense(stopLocation, aa_alignment_map, gene, name)
        else :
            SNPInfo = gene.condensedRegDelInfo()
            res = 0
            for snp in SNPInfo:
                if (aa_alignment_map.get(snp[1]-1,False)) == False:
                    continue
                for mt in snp[2]: 
                    for queryIndex in tuple(aa_alignment_map[snp[1]-1]):
                        if queryIndex == '-': continue
                        if mt == aaQuerySequence[queryIndex]:
                            res = 1
                            break
                    if res == 1: break
                if res == 1: break
            if res == 0 : #take care of Mult
                verifyMultiple(aa_alignment_map, gene, name, aaQuerySequence)
            else: resistant(name, res)


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
output = open("Test/insertion_and_deletion_tests_output.txt", "w")
fileName = ["Insertion1", "Insertion3", "Insertion4_1", "Insertion4_2", "Insertion2_1", "Insertion2_2", "Insertion2_3", "Insertion5", "Deletion1", "Deletion2", "Insertion6", "Deletion3_1", "Deletion3_2", "InsertionDeletion1_1", "InsertionDeletion1_2", "InsertionDeletion2_1", "InsertionDeletion2_2", "InsertionDeletion3", "InsertionDeletion4", "Deletion4", "Deletion5"]#, "Test", "Test2", "P_BPW_50_2_filtered"]
#outputName = ["Test/test_output.txt", "Test/test2_output.txt", "Test/P_BPW_50_2_filtered_output.txt"]
fileNameIndex = 0
outputNameIndex = 0 
for name in fileName:
    output.write(name + "\n")
    pysam.sort("-o", "Test/Sorted_" + name + ".sam", "Test/" + name + ".sam")
    
    samfile = pysam.AlignmentFile("Test/Sorted_" + name + ".sam", "r")

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
    #if (fileNameIndex >= 20) & (fileNameIndex < 23) :
    #    output.close()
    #    output = open(outputName[outputNameIndex], "w")
    #    outputNameIndex += 1
    #fileNameIndex += 1
output.close()