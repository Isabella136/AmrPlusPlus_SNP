from SNP_Verification.Gene import Gene
from SNP_Verification.SNP import SNP
from SNP_Verification import dnaTranslate
from SNP_Verification import reverseTranslation

geneDict = {

}
def SNPTest(SNP_Num, fullName, start, end, cigar, SEQ):
    line1 = []
    line2 = []
    line1.append("Test_" + SNP_Num.__str__())
    line2.append(line1[0])
    line1.append("99")
    line2.append("147")
    line1.append(fullName)
    line2.append(line1[2])
    line1.append(start.__str__())
    line2.append(line1[3])
    line1.append("255")
    line2.append(line1[4])
    line1.append(cigar)
    line2.append(cigar)
    line1.append("=")
    line2.append(line1[6])
    line1.append(line1[3])
    line2.append(line1[3])
    line1.append((end-(start-1)).__str__())
    line2.append("-" + line1[8])
    line1.append(SEQ)
    line2.append(SEQ)
    line1.append("*")
    line2.append(line1[10])
    toReturn = ""
    for tab in line1:
        toReturn = toReturn + tab + "\t"
    toReturn = toReturn + "\n"
    for tab in line2:
        toReturn = toReturn + tab + "\t"
    toReturn = toReturn + "\n"
    return toReturn

def sortSnpInfo(snpInfo):
    snpFirst = snpInfo[0:int(len(snpInfo)/2)]
    if len(snpFirst) > 1:
        snpFirst = sortSnpInfo(snpFirst)
    snpLast = snpInfo[int(len(snpInfo)/2):len(snpInfo)]
    if len(snpLast) > 1:
        snpLast = sortSnpInfo(snpLast)
    toReturn = []
    for i in range(0, len(snpInfo)):
        if len(snpFirst) == 0:
            toReturn.append(snpLast.pop(0))
        elif len(snpLast) == 0:
            toReturn.append(snpFirst.pop(0))
        else:
            if snpFirst[0][1] < snpLast[0][1]:
                toReturn.append(snpFirst.pop(0))
            else:
                toReturn.append(snpLast.pop(0))
    return toReturn


def multipleSNPTest(snpInfo, gene, SNP_Num):
    snpInfo = sortSnpInfo(snpInfo)
    start = snpInfo[0][1]*3-23
    end = snpInfo[-1][1]*3+21
    if (start < 1) : start = 1
    if end > gene.ntSeqLength(): end = gene.ntSeqLength()
    inBetweenH = []
    lastNt = start - 1
    lastIsM = False
    for snp in snpInfo:
        if (lastNt < snp[1]*3 - 3):
            if not(lastIsM): #wt435-, wt437-/mt
                inBetweenH.append((((snp[1]*3-3) - lastNt),'M'))
                lastIsM = True
            else: #wt435mt, wt437-/mt
                inBetweenH[-1] = (inBetweenH[-1][0] + ((snp[1]*3-3) - lastNt), 'M')
        if snp[2][0] == "-":
            if (lastIsM): #wt436mt, wt437-; wt435-/mt, wt437-
                lastIsM = False
                inBetweenH.append((3, 'D'))
            else: #wt436-,wt437-
                inBetweenH[-1] = (inBetweenH[-1][0] + 3, 'D')
        else:
            if (lastIsM): #wt436mt, wt437mt; #wt435-/mt, wt437mt:
                inBetweenH[-1] = (inBetweenH[-1][0] + 3, 'M')
            else: #wt436-,wt437mt
                lastIsM = True
                inBetweenH.append((3, 'M'))
        lastNt = snp[1]*3
    if lastIsM:
        inBetweenH[-1] = (inBetweenH[-1][0] + end - snpInfo[-1][1]*3, 'M')
    elif end > snpInfo[-1][1]*3:
        inBetweenH.append((end - snpInfo[-1][1]*3, 'M'))
    cigar = "20H"
    for next in inBetweenH:
        cigar = cigar + next[0].__str__() + next[1]
    cigar = cigar + "20H"
    sequence = gene.ntSequence()
    SEQ = sequence[start-1:(snpInfo[0][1]*3-3)]
    for i in range(0, len(snpInfo)):
        SEQ = SEQ + reverseTranslation(snpInfo[i][2][0])
        if i < (len(snpInfo) - 1):
            SEQ = SEQ + sequence[snpInfo[i][1]*3:snpInfo[i+1][1]*3-3]
        else:
            SEQ = SEQ + sequence[snpInfo[i][1]*3:end]
    return SNPTest(SNP_Num, gene.getFullName(), start, end, cigar, SEQ)

def singleSNPTest(snpInfo, gene, SNP_Num):
    start = snpInfo[1]*3-23
    end = snpInfo[1]*3+21
    if (start < 1) : start = 1
    if end > gene.ntSeqLength(): end = gene.ntSeqLength()
    cigar = "20H" + (end - (start -1)).__str__() + "M20H"
    if snpInfo[2][0] == "-":
        cigar = "20H"
        if snpInfo[1] != 1:
            cigar = cigar + ((snpInfo[1]*3-3)-(start-1)).__str__() + "M3D"
        if snpInfo[1] != gene.ntSeqLength():
            cigar = cigar + (end-snpInfo[1]*3).__str__() + "M20H"
    sequence = gene.ntSequence()
    SEQ = sequence[start-1:(snpInfo[1]*3-3)] + reverseTranslation(snpInfo[2][0]) + sequence[(snpInfo[1]*3):end]
    return SNPTest(SNP_Num, gene.getFullName(), start, end, cigar, SEQ)

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

SAM_file = open("Test/Test.sam", "w")
SNP_Num = 1
for gene in geneDict.values():
    for snp in gene.condensedRegDelInfo():
        SAM_file.write(singleSNPTest(snp, gene, SNP_Num))
        SNP_Num+=1
    for snp in gene.condensedNonInfo():
        SAM_file.write(singleSNPTest(snp, gene, SNP_Num))
        SNP_Num+=1
    for snp in gene.condensedMultInfo():
        SAM_file.write(multipleSNPTest(snp, gene, SNP_Num))
        SNP_Num+=1
SAM_file.close()
SAM_file = open("Test/Test2.sam", "w")
SNP_Num = 1
SAM_file.write(multipleSNPTest([('Y', 147, ('*',)), ('G', 146, ('*',))], geneDict.get("MEG_2866|Drugs|Mycobacterium_tuberculosis-specific_Drug|Ethionamide-resistant_mutant|ETHA|RequiresSNPConfirmation"), SNP_Num))
SAM_file.write(multipleSNPTest([('Y', 147, ('*',)), ('Y', 140, ('*',))], geneDict.get("MEG_2866|Drugs|Mycobacterium_tuberculosis-specific_Drug|Ethionamide-resistant_mutant|ETHA|RequiresSNPConfirmation"), SNP_Num))
SAM_file.write(multipleSNPTest([('Y', 147, ('*',)), ('C', 131, ('*',))], geneDict.get("MEG_2866|Drugs|Mycobacterium_tuberculosis-specific_Drug|Ethionamide-resistant_mutant|ETHA|RequiresSNPConfirmation"), SNP_Num))
SAM_file.write(multipleSNPTest([('W', 98, ('*',)), ('R', 99, ('*',))], geneDict.get("MEG_7250|Drugs|Mycobacterium_tuberculosis-specific_Drug|Para-aminosalicylic_acid_resistant_mutant|THYA|RequiresSNPConfirmation"), SNP_Num))
SAM_file.write(multipleSNPTest([('W', 98, ('*',)), ('S', 105, ('*',))], geneDict.get("MEG_7250|Drugs|Mycobacterium_tuberculosis-specific_Drug|Para-aminosalicylic_acid_resistant_mutant|THYA|RequiresSNPConfirmation"), SNP_Num))
SAM_file.write(multipleSNPTest([('W', 98, ('*',)), ('A', 115, ('*',))], geneDict.get("MEG_7250|Drugs|Mycobacterium_tuberculosis-specific_Drug|Para-aminosalicylic_acid_resistant_mutant|THYA|RequiresSNPConfirmation"), SNP_Num))
SAM_file.write(multipleSNPTest([('G', 84, ('D',)), ('L', 100, ('P',)), ('A', 115, ('T',)), ('Y', 122, ('N',))], geneDict.get("MEG_4132|Drugs|Mycobacterium_tuberculosis-specific_Drug|Isoniazid-resistant_mutant|NDH|RequiresSNPConfirmation"), SNP_Num))
SAM_file.write(multipleSNPTest([('Q', 10, ('*',)), ('D', 12, ('A',)), ('L', 19, ('P',)), ('C', 14, ('Y',)), ('G', 23, ('V',))], geneDict.get("MEG_5803|Drugs|Mycobacterium_tuberculosis-specific_Drug|Pyrazinamide-resistant_mutant|PNCA|RequiresSNPConfirmation"), SNP_Num))
SAM_file.close()