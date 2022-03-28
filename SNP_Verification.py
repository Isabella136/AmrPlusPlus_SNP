from re import I
from SNP_Verification.Gene import Gene
from SNP_Verification.SNP import SNP


metamarcSNPinfo = open("SNPInfoExtraction/metamarcSNPinfo.fasta", "rt")
genesList = []
isSequence = False
header = ""
sequence = ""
for line in metamarcSNPinfo:
    if isSequence:
        sequence = line
        temp = 0
        for i in range(0, 5):
            temp = header[temp+1:].find('|') + temp + 1
        genesList.append(Gene(header[1:temp], sequence[:-1], header[temp+1:]))
        isSequence = False
    else:
        header = line
        isSequence = True
