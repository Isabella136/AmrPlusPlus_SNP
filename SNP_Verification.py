from SNP_Verification.Gene import Gene
from SNP_Verification.SNP import SNP
from SNP_Verification import dnaTranslate
import pysam

geneDict = {

}

def verify(iter, gene):
    return

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
output = open("test_SNP_Verification_Sorted_Fluro_P_CI_50_3_filtered.csv", "w")
samfile = pysam.AlignmentFile("SAM_files/Sorted_Fluro_P_CI_50_3_filtered.sam", "r")


tid = 0
while (samfile.is_valid_tid(tid)):
    gene = geneDict.get(samfile.get_reference_name(tid), False)
    if (gene == False):
        tid += 1
        continue
    output.write(verify(samfile.fetch(tid=tid), gene))
    tid += 1
samfile.close()   
    
