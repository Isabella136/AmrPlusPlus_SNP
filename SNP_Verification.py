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
    index = 0
    shift = 0
    translatable = False
    for op in cigar:
        if op == "S":
            index += 1
        elif (op == "M") | (op == "X") | (op == "="):
            if not(translatable):
                if (aligned_pair[index][1] % 3) != 0:
                    index += 1
                else:
                    translatable = True
                    queryLength += 1
                    alignment_map.append((op, aligned_pair[index][0], aligned_pair[index][1]), shift)
                    index += 1
            else:
                queryLength += 1
                alignment_map.append((op, aligned_pair[index][0], aligned_pair[index][1]), shift)
                index += 1
        elif op == "I":
            shift += 1
            if not(translatable):
                if ((aligned_pair[index-1][1] + 1) % 3) != 0:
                    index += 1
                else:
                    translatable = True
                    queryLength += 1
                    alignment_map.append((op, aligned_pair[index][0], aligned_pair[index-1][1]+1), shift)
                    index += 1
            else:
                queryLength += 1
                alignment_map.append((op, aligned_pair[index][0], aligned_pair[index-1][1]+1), shift)
                index += 1
        elif op == "D":
            shift -= 1
            if not(translatable):
                if (aligned_pair[index][1] % 3) != 0:
                    index += 1
                else:
                    translatable = True
                    alignment_map.append((op, aligned_pair[index-1][0]+1, aligned_pair[index][1]), shift)
                    index += 1
            else:
                alignment_map.append((op, aligned_pair[index-1][0]+1, aligned_pair[index][1]), shift)
                index += 1
        else:
            continue
    while (len(alignment_map) % 3) != 0:
        alignment_map.pop()
    return alignment_map

def aaAlignment(nt_alignment_map):
    aa_alignment_map = {

    }
    ntIndex = 0
    aaIndex = 0
    lastShift = 0
    previousQuery = []
    nextQuery = []
    disregard = False
    for nt in nt_alignment_map:
        if (ntIndex % 3) == 2:
            if nt[0] == "M":
                if nt[3] == 0:
                    disregard = False
                    if lastShift >= 0:
                        previousQuery.clear()
                        previousQuery.append(int(nt[2]/3))
                        aa_alignment_map[aaIndex] = (previousQuery.copy(), disregard)
                        lastShift = nt[3]
                        aaIndex += 1
                    else:
                        previousQuery.append(int(nt[2]/3))
                        aa_alignment_map[aaIndex] = (previousQuery.copy(), disregard)
                        previousQuery.clear()
                        previousQuery.append(int(nt[2]/3))
                        lastShift = nt[3]
                        aaIndex += 1
                else :
                    disregard = ((nt[3] % 3) == 0)
                    if lastShift < 0:
                        
            elif (nt[0] == "I"):
                lastShift = nt[3]
            else:

        ntIndex += 1


def verify(read, gene):
    name = gene.getName()
    cigarOpCount = read.get_cigar_stats()[0].tolist()
    #if insertions != deletions, disregard
    if cigarOpCount[1] != cigarOpCount[2]: disregard(name) 
    else :
        cigar = extendCigar(read.cigarstring)
        aligned_pair = read.get_aligned_pairs()
        alignment_map = mapCigarToAlignment(cigar, aligned_pair)
        querySeq = read.query_sequence
        trimmedQuerySequence = querySeq[alignment_map[0][1]:alignment_map[len(alignment_map)-1][1]+1]
        aaQuerySequence = dnaTranslate(trimmedQuerySequence)
        stopLocation = aaQuerySequence.find('*') 
        #if missense mutations, disregard
        if (stopLocation != -1) & ((stopLocation < len(aaQuerySequence)-1) | (read.reference_end < gene.ntSeqLength())-1): disregard(name)
        else :
            SNPInfo = gene.condensedInfo()
            res = False
            for snp in SNPInfo:
                if ((snp[1]-1)*3 < alignment_map[0][1]) | ((snp[1]-1)*3 >= alignment_map[len(alignment_map)-1][1]):
                    continue
                for mt in snp[2]: 
                    queryIndex = findQueryIndex(read.get_aligned_pairs(),(snp[1]-1)*3)
                    if (queryIndex != -1) & (mt == aaQuerySequence[int((queryIndex - (queryIndex%3))/3)-1]):
                        res = True
                        resistant(name, 1)
                        break
                if res:
                    break


metamarcSNPinfo = open("SNPInfoExtraction/metamarcSNPinfo.fasta", "rt")
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
output = open("test_SNP_Verification_Sorted_Fluro_P_BPW_10_3_filtered.txt", "w")
samfile = pysam.AlignmentFile("SAM_files/Sorted_Fluro_P_BPW_10_3_filtered.sam", "r")

iter = samfile.fetch()
for read in iter:
    gene = geneDict.get(read.reference_name, False)
    if (gene == False):
        continue
    elif (read.cigarstring == None) :
        continue
    verify(read, gene)
samfile.close() 
for name in argInfoDict:
    output.write(snpInfoPrint(name, str(argInfoDict[name][0]), str(argInfoDict[name][1])))
output.close()
