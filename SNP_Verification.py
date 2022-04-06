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
    refLength = 0
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
                    refLength += 1
                    queryLength += 1
                    alignment_map.append((op, aligned_pair[index][0], aligned_pair[index][1], shift))
                    index += 1
            else:
                refLength += 1
                queryLength += 1
                alignment_map.append((op, aligned_pair[index][0], aligned_pair[index][1], shift))
                index += 1
        elif op == "I":
            queryLength += 1
            shift += 1
            alignment_map.append((op, aligned_pair[index][0], aligned_pair[index][1], shift))
            index += 1
        elif op == "D":
            refLength += 1
            shift -= 1
            alignment_map.append((op, aligned_pair[index][0], aligned_pair[index][1]))
            index += 1
        else:
            continue
    while ((refLength % 3) != 0) &  ((queryLength % 3) != 0):
        poped = alignment_map.pop()
        if poped[0] == "D" :
            refLength -= 1
        elif poped[0] == "I" :
            queryLength -= 1
        else:
            refLength -= 1
            queryLength -= 1
    return alignment_map

def verify(read, gene):
    name = gene.getName()
    cigarOpCount = read.get_cigar_stats()[0].tolist()
    #if insertions != deletions, disregard
    if cigarOpCount[1] != cigarOpCount[2]: disregard(name) 
    else :
        querySeq = read.query_sequence
        cigar = extendCigar(read.cigarstring)
        aligned_pair = read.get_aligned_pairs()
        alignment_map = mapCigarToAlignment(cigar, aligned_pair)
        trimmedQuerySequence = querySeq[alignment_map[0][1]:alignment_map[len(alignment_map)][1]]
        aaQuerySequence = dnaTranslate(trimmedQuerySequence)
        stopLocation = aaQuerySequence.find('*') 
        #if missense mutations, disregard
        if (stopLocation != -1) & ((stopLocation < len(aaQuerySequence)-1) | (read.reference_end < gene.ntSeqLength())-1): disregard(name)
        else :
            SNPInfo = gene.condensedInfo()
            res = False
            for snp in SNPInfo:
                if ((snp[1]-1)*3 < (read.reference_start)) | ((snp[1]-1)*3 >= read.reference_end):
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
test = open("test.txt", "w")

iter = samfile.fetch()
for read in iter:
    gene = geneDict.get(read.reference_name, False)
    if (gene == False):
        continue
    elif (read.cigarstring == None) :
        continue
    verify(read, gene)
samfile.close()
test.close()   
for name in argInfoDict:
    output.write(snpInfoPrint(name, str(argInfoDict[name][0]), str(argInfoDict[name][1])))
output.close()
