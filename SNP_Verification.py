from SNP_Verification.Gene import Gene
from SNP_Verification.SNP import SNP
from SNP_Verification import dnaTranslate
import pysam

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
    translatable = False
    for op in cigar:
        if op == "S":
            index += 1
        elif (op == "M") | (op == "X") | (op == "="):
            prevShift = shift
            if not(translatable):
                if (aligned_pair[index][1] % 3) != 0:
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
                index += 1
            else:
                queryLength += 1
                alignment_map.append((op, aligned_pair[index][0], aligned_pair[index][1], shift, prevShift))
                index += 1
        elif op == "D":
            prevShift = shift
            shift -= 1
            if not(translatable):
                index += 1
            else:
                refLength += 1
                alignment_map.append((op, aligned_pair[index][0], aligned_pair[index][1], shift, prevShift))
                index += 1
        else:
            continue
    while (refLength % 3) != 0:
        poped = alignment_map.pop()
        if poped[0] != "I" : refLength -= 1
    return alignment_map

def aaAlignment(nt_alignment_map):
    aa_alignment_map = {

    }
    insertCount = 0
    ntRefIndex = nt_alignment_map[0][2]
    ntQueryIndex = 0
    aaQueryIndex = 0
    prevAaShift = None    
    firstAlignment = None
    thirdAlignment = 0
    mapIndex = 0
    for nt in nt_alignment_map:
        if (ntRefIndex % 3) == 0:
            if nt[0] == "I":
                insertCount += 1
                if (nt[3] % 3) == 0: 
                    if insertCount == 3: insertCount = 0 #Can only happen if inserts are consecutive
                    else:
                        aa = aa_alignment_map.get(int(ntRefIndex/3-1), False)
                        if aa == False:
                            aa_alignment_map.update({int(ntRefIndex/3)-1:(thirdAlignment,)})
                        else:
                            temp = list(aa)
                            temp.append(thirdAlignment)
                            aa = tuple(temp)
                            aa_alignment_map.update({int(ntRefIndex/3)-1:aa})
                        insertCount = 0
                if prevAaShift == None:
                    prevAaShift = nt[4]
                ntQueryIndex += 1
            elif nt[0] == "D":
                insertCount = 0
                firstAlignment = aaQueryIndex
                ntRefIndex += 1
            else: #nt[0] == "M"
                insertCount = 0
                ntQueryIndex += 1
                ntRefIndex += 1
                if prevAaShift == None:
                    prevAaShift = nt[4]
                firstAlignment = aaQueryIndex
        elif (ntRefIndex % 3) == 1:
            if nt[0] == "I":
                ntQueryIndex += 1
            elif nt[0] == "D":
                ntRefIndex += 1
            else: #nt[0] == "M"
                ntQueryIndex += 1
                ntRefIndex += 1
                if prevAaShift == None: #preceded by "D"
                    prevAaShift = nt[4]
                    firstAlignment = aaQueryIndex
        else: #(ntRefIndex % 3) == 2
            if nt[0] == "I":
                ntQueryIndex += 1
            elif nt[0] == "D":
                ntRefIndex += 1
                if prevAaShift != None:
                    if (prevAaShift % 3) == 0: #1M2D, 2M1D
                        aa_alignment_map.update({int(nt[2]/3):(thirdAlignment,)})
                prevAaShift = None
                firstAlignment = None
            else: #nt[0] == "M"
                ntQueryIndex += 1
                ntRefIndex += 1
                thirdAlignment = aaQueryIndex
                if prevAaShift == None: #preceded by 2 "D"
                    prevAaShift = nt[4]
                    firstAlignment = aaQueryIndex
                if firstAlignment == thirdAlignment:
                    if (prevAaShift % 3) == 0: #3M, 1D1I2M, 1D/2I...2D1M, 2D/1I...1D2M
                        aa_alignment_map.update({int(nt[2]/3):(thirdAlignment,)})
                elif (prevAaShift % 3) == 0:
                    if (nt[4] % 3) == 0: #3Mand3xI
                        aa_alignment_map.update({int(nt[2]/3):(firstAlignment, thirdAlignment)})
                    else: #3Mand1/2I
                        aa_alignment_map.update({int(nt[2]/3):(firstAlignment,)})
                elif (nt[4] % 3) == 0: #2D/1I...1M1D1M
                    aa_alignment_map.update({int(nt[2]/3):(thirdAlignment,)})    
                prevAaShift = None
                firstAlignment = None
        aaQueryIndex = int(ntQueryIndex / 3)
        mapIndex += 1
    return aa_alignment_map    


def verify(read, gene):
    name = gene.getName()
    cigarOpCount = read.get_cigar_stats()[0].tolist()
    #if insertions != deletions, disregard
    if cigarOpCount[1] != cigarOpCount[2]: disregard(name) 
    else :
        cigar = extendCigar(read.cigarstring)
        aligned_pair = read.get_aligned_pairs()
        nt_alignment_map = mapCigarToAlignment(cigar, aligned_pair)
        querySeq = read.query_sequence
        trimmedQuerySequence = querySeq[nt_alignment_map[0][1]:nt_alignment_map[len(nt_alignment_map)-1][1]+1]
        aaQuerySequence = dnaTranslate(trimmedQuerySequence)
        stopLocation = aaQuerySequence.find('*') 
        #if missense mutations, disregard
        if (stopLocation != -1) & ((stopLocation < len(aaQuerySequence)-1) | (read.reference_end < gene.ntSeqLength())-1): disregard(name)
        else :
            aa_alignment_map = aaAlignment(nt_alignment_map)
            SNPInfo = gene.condensedInfo()
            res = 0
            for snp in SNPInfo:
                if (aa_alignment_map.get(snp[1]-1,False)) == False:
                    continue
                for mt in snp[2]: 
                    for queryIndex in tuple(aa_alignment_map[snp[1]-1]):
                        if mt == aaQuerySequence[queryIndex]:
                            res = 1
                            break
                    if res == 1: break
                if res == 1: break
            resistant(name, res)


metamarcSNPinfo = open("extracted_SNP_files/metamarcSNPinfo.fasta", "rt")
isSequence = False
name = ""
snp = ""
sequence = ""
for line in metamarcSNPinfo:
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
metamarcSNPinfo.close()
pysam.sort("-o", "SAM_files/Sorted_Filtered_out_P_BPW_50_2_R.amr.alignment.sam", "SAM_files/Filtered_out_P_BPW_50_2_R.amr.alignment.sam")
output = open("Test/test_Sorted_Filtered_out_P_BPW_50_2_R.amr.alignment.txt", "w")
samfile = pysam.AlignmentFile("SAM_files/Sorted_Filtered_out_P_BPW_50_2_R.amr.alignment.sam", "r")

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
output.close()
