from SNP_Verification.Gene import Gene
from SNP_Verification.SNP import SNP
from SNP_Verification import dnaTranslate

geneDict = {

}

metamarcSNPinfo = open("extracted_SNP_files/SNPinfo.fasta", "rt")
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