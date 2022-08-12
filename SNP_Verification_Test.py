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
from SNP_Verification_Tools import dnaTranslate
from SNP_Verification_Tools import reverseTranslation

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

def makeTest(SNP_Num, fullName, aa_seq, cigar, start, end):
    nt_seq = ""
    for aa in aa_seq:
        nt_seq = nt_seq + reverseTranslation(aa)
    indexStart = (start-1) % 3
    indexEnd = end % 3
    nt_seq = nt_seq[indexStart:]
    if indexEnd != 0:
        indexEnd = 3 - indexEnd
        nt_seq = nt_seq[:-1*indexEnd]
    return SNPTest(SNP_Num, fullName, start, end, cigar, nt_seq)

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
            if (snpFirst[0][1][0] if type(snpFirst[0][1]) == list else snpFirst[0][1]) < (snpLast[0][1][0] if type(snpLast[0][1]) == list else snpLast[0][1]):
                toReturn.append(snpFirst.pop(0))
            else:
                toReturn.append(snpLast.pop(0))
    return toReturn

def nucleicMustTest(mustInfo, gene, SNP_Num):
    start = mustInfo.getPos()-26
    if start < 1: start = 1
    sequence = gene.ntSequence()
    SEQ = sequence[start-1:(mustInfo.getPos()-1)]
    cigar = "20H"
    lastNt = mustInfo.getPos()-1
    while mustInfo != None:
        if lastNt < mustInfo.getPos()-1:
            SEQ = SEQ + sequence[lastNt:mustInfo.getPos()-1]
        SEQ = SEQ + mustInfo.getWt()
        lastNt = mustInfo.getPos()
        end = mustInfo.getPos()+ 25
        mustInfo = mustInfo.getNext()
    if end > gene.ntSeqLength(): end = gene.ntSeqLength()
    cigar = cigar + str(end - (start - 1)) + "M20H"
    SEQ = SEQ + sequence[lastNt:end]
    return SNPTest(SNP_Num, gene.getFullName(), start, end, cigar, SEQ)

def mustTest(mustInfo, gene, SNP_Num):
    if gene.rRna(): return nucleicMustTest(mustInfo, gene, SNP_Num)
    start = mustInfo.getPos()*3-23
    if start < 1: start = 1
    sequence = gene.ntSequence()
    SEQ = sequence[start-1:(mustInfo.getPos()*3-3)]
    cigar = "20H"
    lastNt = mustInfo.getPos()*3-3
    while mustInfo != None:
        if lastNt < mustInfo.getPos()*3-3:
            SEQ = SEQ + sequence[lastNt:mustInfo.getPos()*3-3]
        SEQ = SEQ + reverseTranslation(mustInfo.getWt())
        lastNt = mustInfo.getPos()*3
        end = mustInfo.getPos()*3 + 21
        mustInfo = mustInfo.getNext()
    if end > gene.ntSeqLength(): end = gene.ntSeqLength()
    cigar = cigar + str(end - (start - 1)) + "M20H"
    SEQ = SEQ + sequence[lastNt:end]
    return SNPTest(SNP_Num, gene.getFullName(), start, end, cigar, SEQ)

def multipleSNPTest(snpInfo, gene, SNP_Num):
    if gene.rRna(): return multipleNucleicSNPTest(snpInfo, gene, SNP_Num)
    snpInfo = sortSnpInfo(snpInfo)
    start = (snpInfo[0][1][0] if type(snpInfo[0][1]) == list else snpInfo[0][1])*3-23
    end = (snpInfo[-1][1][0] if type(snpInfo[-1][1]) == list else snpInfo[-1][1])*3+21
    if (start < 1) : start = 1
    if end > gene.ntSeqLength(): end = gene.ntSeqLength()
    inBetweenH = []
    lastNt = start - 1
    lastIsM = False
    lastIsI = False
    lastIsD = False
    for snp in snpInfo:
        if (lastNt < (snp[1][0] if type(snp[1]) == list else snp[1])*3 - 3):
            if not(lastIsM): #wt435-, wt437-/mt
                inBetweenH.append(((((snp[1][0] if type(snp[1]) == list else snp[1])*3-3) - lastNt),'M'))
                lastIsM = True
                lastIfI = False
                lastIsD = False
            else: #wt435mt, wt437-/mt
                inBetweenH[-1] = (inBetweenH[-1][0] + (((snp[1][0] if type(snp[1]) == list else snp[1])*3-3) - lastNt), 'M')
        if snp[2][0] == "-":
            if not(lastIsD): #wt436mt, wt437-; wt435-/mt, wt437-
                lastIsD = True
                lastIsI = False
                lastIsM = False
                inBetweenH.append((3, 'D'))
            else: #wt436-,wt437-
                inBetweenH[-1] = (inBetweenH[-1][0] + 3, 'D')
        elif snp[2][0] == "+":
            if not(lastIsI): #wt436mt, wt437-; wt435-/mt, wt437-
                lastIsM = False
                lastIsD = False
                lastIsI = True
                inBetweenH.append((3 * len(snp[0]), 'I'))
            else:
                inBetweenH[-1] = (inBetweenH[-1][0] + 3 * len(snp[0]), 'I')

        else:
            if (lastIsM): #wt436mt, wt437mt; #wt435-/mt, wt437mt:
                inBetweenH[-1] = (inBetweenH[-1][0] + 3, 'M')
            else: #wt436-,wt437mt
                lastIsM = True
                lastIsD = False
                lastIsI = False
                inBetweenH.append((3, 'M'))
        if (type(snp[1]) == list):
            lastNt = snp[1][0]*3
        else:
            lastNt = snp[1]*3
    if lastIsM:
        inBetweenH[-1] = (inBetweenH[-1][0] + end - (snpInfo[-1][1][0] if type(snpInfo[-1][1]) == list else snpInfo[-1][1])*3, 'M')
    elif end > (snpInfo[-1][1][0] if type(snpInfo[-1][1]) == list else snpInfo[-1][1])*3:
        inBetweenH.append((end - (snpInfo[-1][1][0] if type(snpInfo[-1][1]) == list else snpInfo[-1][1])*3, 'M'))
    cigar = "20H"
    for next in inBetweenH:
        cigar = cigar + str(next[0]) + next[1]
    cigar = cigar + "20H"
    sequence = gene.ntSequence()
    SEQ = sequence[start-1:((snpInfo[0][1][0] if type(snpInfo[0][1]) == list else snpInfo[0][1])*3-3)]
    for i in range(0, len(snpInfo)):
        if snpInfo[i][2][-1] != '-':
            if snpInfo[i][2][-1] == '+':
                for aa in snpInfo[i][0]:
                    SEQ = SEQ + reverseTranslation(aa)
            else:
                SEQ = SEQ + reverseTranslation(snpInfo[i][2][-1])
        if i < (len(snpInfo) - 1):
            SEQ = SEQ + sequence[(snpInfo[i][1][0] if type(snpInfo[i][1]) == list else snpInfo[i][1])*3:(snpInfo[i+1][1][0] if type(snpInfo[i+1][1]) == list else snpInfo[i+1][1])*3-3]
        else:
            SEQ = SEQ + sequence[(snpInfo[i][1][0] if type(snpInfo[i][1]) == list else snpInfo[i][1])*3:end]
    return SNPTest(SNP_Num, gene.getFullName(), start, end, cigar, SEQ)

def multipleNucleicSNPTest(snpInfo, gene, SNP_Num):
    snpInfo = sortSnpInfo(snpInfo)
    start = snpInfo[0][1]-26
    end = snpInfo[-1][1]+25
    if (start < 1) : start = 1
    if end > gene.ntSeqLength(): end = gene.ntSeqLength()
    inBetweenH = []
    lastNt = start - 1
    lastIsM = False
    for snp in snpInfo:
        if (lastNt < (snp[1][0] if type(snp[1]) == list else snp[1])-1):
            if not(lastIsM): #wt435-, wt437-/mt
                inBetweenH.append((((snp[1][0] if type(snp[1]) == list else snp[1])-1 - lastNt),'M'))
                lastIsM = True
            else: #wt435mt, wt437-/mt
                inBetweenH[-1] = (inBetweenH[-1][0] + ((snp[1][0] if type(snp[1]) == list else snp[1])-1 - lastNt), 'M')
        if snp[2][0] == "-":
            if (lastIsM): #wt436mt, wt437-; wt435-/mt, wt437-
                lastIsM = False
                inBetweenH.append((1, 'D'))
            else: #wt436-,wt437-
                inBetweenH[-1] = (inBetweenH[-1][0] + 1, 'D')
        else:
            if (lastIsM): #wt436mt, wt437mt; #wt435-/mt, wt437mt:
                inBetweenH[-1] = (inBetweenH[-1][0] + 1, 'M')
            else: #wt436-,wt437mt
                lastIsM = True
                inBetweenH.append((1, 'M'))
        if (type(snp[1]) == list):
            lastNt = snp[1][0]
        else:
            lastNt = snp[1]
    if lastIsM:
        inBetweenH[-1] = (inBetweenH[-1][0] + end - (snpInfo[-1][1][0] if type(snpInfo[-1][1]) == list else snpInfo[-1][1]), 'M')
    elif end > (snpInfo[-1][1][0] if type(snpInfo[-1][1]) == list else snpInfo[-1][1]):
        inBetweenH.append((end - (snpInfo[-1][1][0] if type(snpInfo[-1][1]) == list else snpInfo[-1][1]), 'M'))
    cigar = "20H"
    for next in inBetweenH:
        cigar = cigar + str(next[0]) + next[1]
    cigar = cigar + "20H"
    sequence = gene.ntSequence()
    SEQ = sequence[start-1:((snpInfo[0][1][0] if type(snpInfo[0][1]) == list else snpInfo[0][1])-1)]
    for i in range(0, len(snpInfo)):
        if snpInfo[i][2][-1] != '-':
            SEQ = SEQ + snpInfo[i][2][-1]
        if i < (len(snpInfo) - 1):
            SEQ = SEQ + sequence[(snpInfo[i][1][0] if type(snpInfo[i][1]) == list else snpInfo[i][1]):(snpInfo[i+1][1][0] if type(snpInfo[i+1][1]) == list else snpInfo[i+1][1])-1]
        else:
            SEQ = SEQ + sequence[(snpInfo[i][1][0] if type(snpInfo[i][1]) == list else snpInfo[i][1]):end]
    return SNPTest(SNP_Num, gene.getFullName(), start, end, cigar, SEQ)

def deletionTest(snpInfo, gene, SNP_Num):
    toReturn = ""
    for i in snpInfo[1]:
        start = i*3-23
        end = i*3+21
        if (start < 1) : start = 1
        if end > gene.ntSeqLength(): end = gene.ntSeqLength()
        cigar = "20H"
        if i != 1:
            cigar = cigar + str((i*3-3)-(start-1))+ "M3D"
        if i != gene.ntSeqLength():
            cigar = cigar + str(end-i*3) + "M20H"
        sequence = gene.ntSequence()
        SEQ = sequence[start-1:(i*3-3)] + sequence[(i*3):end]
        toReturn = toReturn + SNPTest(str(SNP_Num) +"a"+str(i), gene.getFullName(), start, end, cigar, SEQ)
        cigar = "20H"
        if i != 1:
            cigar = cigar + str((i*3-4)-(start-1)) + "M3D"
        if i != gene.ntSeqLength():
            cigar = cigar + str(end-(i*3)+1) + "M20H"
        sequence = gene.ntSequence()
        SEQ = sequence[start-1:(i*3-4)] + sequence[(i*3)-1:end]
        toReturn = toReturn + SNPTest(str(SNP_Num) +"a"+str(i), gene.getFullName(), start, end, cigar, SEQ)
        cigar = "20H"
        if i != 1:
            cigar = cigar + str((i*3-2)-(start-1)) + "M3D"
        if i != gene.ntSeqLength():
            cigar = cigar + str(end-(i*3)-1) + "M20H"
        sequence = gene.ntSequence()
        SEQ = sequence[start-1:(i*3-2)] + sequence[(i*3)+1:end]
        toReturn = toReturn + SNPTest(str(SNP_Num) +"a"+str(i), gene.getFullName(), start, end, cigar, SEQ)
    return toReturn

def insertionTest(snpInfo, gene, SNP_Num):
    toReturn = ""
    for i in snpInfo[1]:
        start = i*3-23
        end = i*3+21
        if (start < 1) : start = 1
        if end > gene.ntSeqLength(): end = gene.ntSeqLength()
        cigar = "20H"
        if i != 1:
            cigar = cigar + str((i*3-3)-(start-1))+ "M" + str(3 * len(snpInfo[0])) + "I"
        if i != gene.ntSeqLength():
            cigar = cigar + str(end-i*3) + "M20H"
        sequence = gene.ntSequence()
        inserted = ""
        for aa in snpInfo[0]:
            inserted = inserted + reverseTranslation(aa)
        SEQ = sequence[start-1:(i*3-3)] + inserted + sequence[(i*3):end]
        toReturn = toReturn + SNPTest(str(SNP_Num) +"a"+str(i), gene.getFullName(), start, end, cigar, SEQ)
    if (len(snpInfo[1]) > 1) and (len(snpInfo[0]) > 1):
        start = snpInfo[1][0]*3-20
        end = snpInfo[1][0]*3+24
        if (start < 1) : start = 1
        if end > gene.ntSeqLength(): end = gene.ntSeqLength()
        cigar = "20H"
        if snpInfo[1][0] != 1:
            cigar = cigar + str((snpInfo[1][0]*3)-(start-1))+ "M" + 3 * len(snpInfo[0]) + "I"
        if snpInfo[1][0] != gene.ntSeqLength():
            cigar = cigar + str(end-snpInfo[1][0]*3+3) + "M20H"
        sequence = gene.ntSequence()
        inserted = ""
        for aa in snpInfo[0][1:]:
            inserted = inserted + reverseTranslation(aa)
        inserted = inserted + reverseTranslation(snpInfo[0][0])
        SEQ = sequence[start-1:(snpInfo[1][0]*3)] + inserted + sequence[(snpInfo[1][0]*3+3):end]
        toReturn = toReturn + SNPTest(str(SNP_Num) +"a"+str(snpInfo[1]), gene.getFullName(), start, end, cigar, SEQ)
    return toReturn

def singleSNPTest(snpInfo, gene, SNP_Num):
    if gene.rRna(): return singleNucleicSNPTest(snpInfo, gene, SNP_Num)
    if snpInfo[2] == "-": return deletionTest(snpInfo, gene, SNP_Num)
    elif snpInfo[2] == "+": return insertionTest(snpInfo, gene, SNP_Num)
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
    SEQ = sequence[start-1:(snpInfo[1]*3-3)] + reverseTranslation(snpInfo[2][-1]) + sequence[(snpInfo[1]*3):end]
    return SNPTest(SNP_Num, gene.getFullName(), start, end, cigar, SEQ)

def singleNucleicSNPTest(snpInfo, gene, SNP_Num):
    start = snpInfo[1] - 26
    end = snpInfo[1] + 25
    if start < 1: start = 1
    if end > gene.ntSeqLength(): end = gene.ntSeqLength()
    cigar = "20H" + (end - (start -1)).__str__() + "M20H"
    if snpInfo[2][0] == "-":
        cigar = "20H"
        if snpInfo[1] != 1:
            cigar = cigar + ((snpInfo[1])-(start-1)).__str__() + "M1D"
        if snpInfo[1] != gene.ntSeqLength():
            cigar = cigar + (end-snpInfo[1]).__str__() + "M20H"
    sequence = gene.ntSequence()
    SEQ = sequence[start-1:(snpInfo[1]-1)] + snpInfo[2][-1] + sequence[(snpInfo[1]):end]
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
    for snp in gene.condensedMisInDelInfo():
        SAM_file.write(singleSNPTest(snp, gene, SNP_Num))
        SNP_Num+=1
    for snp in gene.condensedNonInfo():
        SAM_file.write(singleSNPTest(snp, gene, SNP_Num))
        SNP_Num+=1
    for snp in gene.condensedMultInfo():
        SAM_file.write(multipleSNPTest(snp[0], gene, SNP_Num))
        SNP_Num+=1
    must = gene.getFirstMustBetweenParams(1, gene.ntSeqLength())
    if must == None:
        continue
    SAM_file.write(mustTest(must[0], gene, SNP_Num))
    SNP_Num+=1
    if must[0].getNext() != None:
        SAM_file.write(mustTest(must[0].getNext(), gene, SNP_Num))
        SNP_Num+=1

SAM_file.close()
#SAM_file = open("Test/Test2.sam", "w")
#SNP_Num = 1
#SAM_file.write(multipleSNPTest([('Y', 147, ('*',)), ('G', 146, ('*',))], geneDict.get("MEG_2866|Drugs|Mycobacterium_tuberculosis-specific_Drug|Ethionamide-resistant_mutant|ETHA|RequiresSNPConfirmation"), SNP_Num))
#SNP_Num+=1
#SAM_file.write(multipleSNPTest([('Y', 147, ('*',)), ('Y', 140, ('*',))], geneDict.get("MEG_2866|Drugs|Mycobacterium_tuberculosis-specific_Drug|Ethionamide-resistant_mutant|ETHA|RequiresSNPConfirmation"), SNP_Num))
#SNP_Num+=1
#SAM_file.write(multipleSNPTest([('Y', 147, ('*',)), ('C', 131, ('*',))], geneDict.get("MEG_2866|Drugs|Mycobacterium_tuberculosis-specific_Drug|Ethionamide-resistant_mutant|ETHA|RequiresSNPConfirmation"), SNP_Num))
#SNP_Num+=1
#SAM_file.write(multipleSNPTest([('W', 98, ('*',)), ('R', 99, ('*',))], geneDict.get("MEG_7250|Drugs|Mycobacterium_tuberculosis-specific_Drug|Para-aminosalicylic_acid_resistant_mutant|THYA|RequiresSNPConfirmation"), SNP_Num))
#SNP_Num+=1
#SAM_file.write(multipleSNPTest([('W', 98, ('*',)), ('S', 105, ('*',))], geneDict.get("MEG_7250|Drugs|Mycobacterium_tuberculosis-specific_Drug|Para-aminosalicylic_acid_resistant_mutant|THYA|RequiresSNPConfirmation"), SNP_Num))
#SNP_Num+=1
#SAM_file.write(multipleSNPTest([('W', 98, ('*',)), ('A', 115, ('*',))], geneDict.get("MEG_7250|Drugs|Mycobacterium_tuberculosis-specific_Drug|Para-aminosalicylic_acid_resistant_mutant|THYA|RequiresSNPConfirmation"), SNP_Num))
#SNP_Num+=1
#SAM_file.write(multipleSNPTest([('G', 84, ('D',)), ('L', 100, ('P',)), ('A', 115, ('T',)), ('Y', 122, ('N',))], geneDict.get("MEG_4132|Drugs|Mycobacterium_tuberculosis-specific_Drug|Isoniazid-resistant_mutant|NDH|RequiresSNPConfirmation"), SNP_Num))
#SNP_Num+=1
#SAM_file.write(multipleSNPTest([('Q', 10, ('*',)), ('D', 12, ('A',)), ('L', 19, ('P',)), ('C', 14, ('Y',)), ('G', 23, ('V',))], geneDict.get("MEG_5803|Drugs|Mycobacterium_tuberculosis-specific_Drug|Pyrazinamide-resistant_mutant|PNCA|RequiresSNPConfirmation"), SNP_Num))
#SAM_file.close()

FullName = []
FullName.append("MEG_4094|Drugs|Fosfomycin|Fosfomycin_target_mutation|MURA|RequiresSNPConfirmation")
FullName.append("MEG_4130|Drugs|Mycobacterium_tuberculosis-specific_Drug|Isoniazid-resistant_mutant|NDH|RequiresSNPConfirmation")
FullName.append("MEG_5328|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation")
FullName.append("MEG_5331|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation")
FullName.append("MEG_5401|Drugs|betalactams|Penicillin_binding_protein|PBP2|RequiresSNPConfirmation")
FullName.append("MEG_413|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_RND_efflux_regulator|ACRR|RequiresSNPConfirmation")
FullName.append("MEG_6090|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation")

#SAM_file = open("Test/Insertion1.sam", "w")
#SNP_Num = 1
#SAM_file.write(makeTest(SNP_Num, FullName[0], "GLTL", "20H1M3I8M20H", 65 * 3 - 2, 67 * 3)) #Test with mt second
#SNP_Num  += 1
#SAM_file.write(makeTest(SNP_Num, FullName[1], "ALII", "20H1M3I8M20H", 18 * 3 - 2, 20 * 3)) #Test with mt first
#SNP_Num  += 1
#SAM_file.write(makeTest(SNP_Num, FullName[2], "WWSI", "20H1M3I8M20H", 83 * 3 - 2, 85 * 3)) #Test with mt in both
#SNP_Num  += 1
#SAM_file.write(makeTest(SNP_Num, FullName[3], "VINN", "20H1M3I8M20H", 103 * 3 - 2, 105 * 3)) #Test with wt first, mt second
#SNP_Num  += 1
#SAM_file.write(makeTest(SNP_Num, FullName[4], "TVIA", "20H1M3I8M20H", 316 * 3 - 2, 318 * 3)) #Test with mt first, wt second
#SAM_file.close()


#SAM_file = open("Test/Insertion3.sam", "w")
#SNP_Num = 1
#SAM_file.write(makeTest(SNP_Num, FullName[0], "DGLTL", "20H3M1I3M2I6M20H", 64 * 3 - 2, 67 * 3)) #Test with mt second
#SNP_Num  += 1
#SAM_file.write(makeTest(SNP_Num, FullName[1], "VALII", "20H3M1I3M2I6M20H", 17 * 3 - 2, 20 * 3)) #Test with mt first
#SNP_Num  += 1
#SAM_file.write(makeTest(SNP_Num, FullName[2], "DWWSI", "20H3M1I3M2I6M20H", 82 * 3 - 2, 85 * 3)) #Test with mt in both
#SNP_Num  += 1
#SAM_file.write(makeTest(SNP_Num, FullName[3], "TVINN", "20H3M1I3M2I6M20H", 102 * 3 - 2, 105 * 3)) #Test with wt first, mt second
#SNP_Num  += 1
#SAM_file.write(makeTest(SNP_Num, FullName[4], "FTVIA", "20H3M1I3M2I6M20H", 315 * 3 - 2, 318 * 3)) #Test with mt first, wt second
#SAM_file.close()


#SAM_file = open("Test/Insertion4_1.sam", "w")
#SNP_Num = 1
#SAM_file.write(makeTest(SNP_Num, FullName[0], "ADGLT", "20H3M1I6M2I3M20H", 63 * 3 - 2, 66 * 3)) #Test with mt second; second insertion chunk
#SNP_Num  += 1
#SAM_file.write(makeTest(SNP_Num, FullName[1], "RVALI", "20H3M1I6M2I3M20H", 16 * 3 - 2, 19 * 3)) #Test with mt first; second insertion chunk
#SNP_Num  += 1
#SAM_file.write(makeTest(SNP_Num, FullName[2], "GDWWS", "20H3M1I6M2I3M20H", 81 * 3 - 2, 84 * 3)) #Test with mt in both; second insertion chunk
#SNP_Num  += 1
#SAM_file.write(makeTest(SNP_Num, FullName[3], "RTVIN", "20H3M1I6M2I3M20H", 101 * 3 - 2, 104 * 3)) #Test with wt first, mt second; second insertion chunk
#SNP_Num  += 1
#SAM_file.write(makeTest(SNP_Num, FullName[4], "PFTVI", "20H3M1I6M2I3M20H", 314 * 3 - 2, 317 * 3)) #Test with mt first, wt second; second insertion chunk
#SAM_file.close()

#SAM_file = open("Test/Insertion4_2.sam", "w")
#SNP_Num = 1
#SAM_file.write(makeTest(SNP_Num, FullName[0], "DGLTL", "20H3M1I6M2I3M20H", 64 * 3 - 2, 67 * 3)) #Test with mt second; first insertion chunk
#SNP_Num  += 1
#SAM_file.write(makeTest(SNP_Num, FullName[1], "VALII", "20H3M1I6M2I3M20H", 17 * 3 - 2, 20 * 3)) #Test with mt first; first insertion chunk
#SNP_Num  += 1
#SAM_file.write(makeTest(SNP_Num, FullName[2], "DWWSI", "20H3M1I6M2I3M20H", 82 * 3 - 2, 85 * 3)) #Test with mt in both; first insertion chunk
#SNP_Num  += 1
#SAM_file.write(makeTest(SNP_Num, FullName[3], "TVINN", "20H3M1I6M2I3M20H", 102 * 3 - 2, 105 * 3)) #Test with wt first, mt second; first insertion chunk
#SNP_Num  += 1
#SAM_file.write(makeTest(SNP_Num, FullName[4], "FTVIA", "20H3M1I6M2I3M20H", 315 * 3 - 2, 318 * 3)) #Test with mt first, wt second; first insertion chunk
#SAM_file.close()


#SAM_file = open("Test/Insertion2_1.sam", "w")
#SNP_Num = 1
#SAM_file.write(makeTest(SNP_Num, FullName[5], "TTCGA", "20H3M3I9M20H", 44 * 3 - 2, 47 * 3)) #Test with mt after insertion
#SAM_file.close()

#SAM_file = open("Test/Insertion2_2.sam", "w")
#SAM_file.write(makeTest(SNP_Num, FullName[5], "TRCGA", "20H3M3I9M20H", 44 * 3 - 2, 47 * 3)) #Test with wt in insertion, mt after insertion
#SAM_file.close()

#SAM_file = open("Test/Insertion2_3.sam", "w")
#SAM_file.write(makeTest(SNP_Num, FullName[5], "TCRGA", "20H3M3I9M20H", 44 * 3 - 2, 47 * 3)) #Test with mt in insertion, wt after insertion
#SAM_file.close()


#SAM_file = open("Test/Insertion5.sam", "w")
#SAM_file.write(makeTest(SNP_Num, FullName[5], "RCGA", "20H1M3I7M20H", 45 * 3 - 1, 47 * 3)) #Test with mt in insertion
#SAM_file.close()


#SAM_file = open("Test/Deletion1.sam", "w")
#SAM_file.write(makeTest(SNP_Num, FullName[6], "HQN", "20H1M9D8M20H", 432 * 3 - 2, 437 * 3)) #Test with deletion mutations
#SNP_Num  += 1
#SAM_file.write(makeTest(SNP_Num, FullName[5], "CGA", "20H1M3D8M20H", 44 * 3 - 2, 47 * 3)) #Test with deletion making mt
#SAM_file.close()


#SAM_file = open("Test/Deletion2.sam", "w")
SNP_Num = 1
#SAM_file.write(makeTest(SNP_Num, FullName[5], "GVTC", "20H1M2D8M1D2M20H", 41 * 3 - 1, 45 * 3)) #Test with mt at last deletion chunk
#SAM_file.close()


#SAM_file = open("Test/Insertion6.sam", "w")
#SAM_file.write(makeTest(SNP_Num, FullName[5], "GVTCG", "20H1M1I6M2I4M20H", 43 * 3 - 1, 46 * 3)) #Test with mt last insertion chunk
#SAM_file.close()


#SAM_file = open("Test/Deletion3_1.sam", "w")
#SAM_file.write(makeTest(SNP_Num, FullName[5], "CGAI", "20H1M1D8M2D3M20H", 45 * 3 - 2, 49 * 3)) #Test with mt at first deletion chunk
#SAM_file.close()

#SAM_file = open("Test/Deletion3_2.sam", "w")
#SAM_file.write(makeTest(SNP_Num, FullName[5], "GVTC", "20H1M1D8M2D3M20H", 41 * 3 - 2, 45 * 3)) #Test with mt at last deletion chunk
#SAM_file.close()


#SAM_file = open("Test/InsertionDeletion1_1.sam", "w")
#SAM_file.write(makeTest(SNP_Num, FullName[5], "TCGAI", "20H3M1D9M1I2M20H", 44 * 3 - 2, 48 * 3)) #Test with mt at deletion chunk
#SAM_file.close()

#SAM_file = open("Test/InsertionDeletion1_2.sam", "w")
#SAM_file.write(makeTest(SNP_Num, FullName[5], "AGVTC", "20H3M1D9M1I2M20H", 41 * 3 - 2, 45 * 3)) #Test with mt at insertion chunk
#SAM_file.close()


#SAM_file = open("Test/InsertionDeletion2_1.sam", "w")
#SAM_file.write(makeTest(SNP_Num, FullName[5], "CGAI", "20H1M1I8M1D2M20H", 45 * 3 - 2, 48 * 3)) #Test with mt at insertion chunk
#SAM_file.close()

#SAM_file = open("Test/InsertionDeletion2_2.sam", "w")
#SAM_file.write(makeTest(SNP_Num, FullName[5], "GVTC", "20H1M1I8M1D2M20H", 42 * 3 - 2, 45 * 3)) #Test with mt at deletion chunk
#SAM_file.close()


#SAM_file = open("Test/InsertionDeletion3.sam", "w")
#SAM_file.write(makeTest(SNP_Num, FullName[5], "C", "20H1D1M1I1M20H", 45 * 3 - 2, 45 * 3)) 
#SAM_file.close()


#SAM_file = open("Test/InsertionDeletion4.sam", "w")
#SAM_file.write(makeTest(SNP_Num, FullName[5], "C", "20H1I1M1D1M20H", 45 * 3 - 2, 45 * 3)) 
#SAM_file.close()


#SAM_file = open("Test/Deletion4.sam", "w")
#SAM_file.write(makeTest(SNP_Num, FullName[6], "FNP", "20H1M12D6M20H", 433 * 3, 439 * 3)) #Test with deletion mutations
#SAM_file.close()


#SAM_file = open("Test/Deletion5.sam", "w")
#SAM_file.write(makeTest(SNP_Num, FullName[6], "FN", "20H3M12D1M20H", 433 * 3 - 2, 438 * 3 - 2)) #Test with deletion mutations
#SAM_file.close()