from Bio import pairwise2

SNP = {
    }
sortedSNP = {
    }

def makeSNPdict(snp):
    pos = 0
    wt = ""
    if snp == "L525S": snp = "L524S"
    elif snp[:4] == "D515":
        temp = snp[4:]
        snp = "D516" + temp
    elif snp == "T516I": snp = "T525I"
    elif snp == "D525Y": snp = "D516Y"
    elif snp == "D518H": snp = "N518H"
    if snp[0] == '-':
        snp = snp[1:]
        if snp[0].isnumeric():
            pos = (int) (snp[:-1])
            wt = snp[-1:]
        else:
            wt = snp[0]
            pos = (int) (snp[1:])
    else:
        wt = snp[0]
        snp = snp[1:]
        if snp[-2].isnumeric() :
            pos = (int) (snp[:-1])
        else:
            pos = (int) (snp[:-5])
    dictEntry = SNP.get(pos, False)
    if (dictEntry == False):
        SNP.update({pos:wt})
    elif (dictEntry.find(wt) == -1):
        SNP[pos] += wt

MEG_6090_sequence = open("MEG_6090_sequence.fasta", "r")
MEG_6090_sequence.readline()
seq = MEG_6090_sequence.readline()
MEG_6090_sequence.close()
MEG_6090_kargva = open("MEG_6090_kargva.txt", "r")
for line in MEG_6090_kargva:
    headerList = line.split("|")
    if (headerList[1].find(";") != -1):
        snpList = headerList[1].split(";")
        for snp in snpList:
            makeSNPdict(snp)
    else:
        makeSNPdict(headerList[1])
MEG_6090_kargva.close()
posList = []
for pos in SNP:
    posList.append(pos)
posList.sort()
for pos in posList:
    sortedSNP.update({pos:SNP[pos]})
begin = 1
SNPsequence = ""
for pos, wt in sortedSNP.items():
    for i in range (begin,pos):
        SNPsequence += "X"
    SNPsequence += wt
    begin = pos+1
for i in range(begin, len(seq) + 1):
    SNPsequence += "X"
toReturn = open("alignment.txt", "w")
for a in pairwise2.align.globalxs(SNPsequence, seq, -.5, -.1):
    toReturn.write(pairwise2.format_alignment(*a))
toReturn.close()

        
