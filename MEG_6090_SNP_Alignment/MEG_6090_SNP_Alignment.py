from Bio import Align
from Bio import pairwise2
from Bio.Align import substitution_matrices

E_coli_pos_SNP = {
    }
M_tuberculosis_MEG_6090_SNP = []
E_coli_MEG_6090_SNP = []
E_coli_to_M_tuber = {
    }
E_coli_to_M_tuber_SNP = {
    }

def makeSNPdict(snp):
    pos = 0
    wt = ""
    wt = snp[0]
    snp = snp[1:]
    if snp[-2].isnumeric() :
        pos = (int) (snp[:-1])
    else:
        pos = (int) (snp[:-4])
    dictEntry = E_coli_pos_SNP.get(pos, False)
    if (dictEntry == False):
        E_coli_pos_SNP.update({pos:wt})
    elif (dictEntry.find(wt) == -1):
        E_coli_pos_SNP[pos] += wt

def changePos(snp):
    pos = 0
    mt = ""
    wt = snp[0]
    snp = snp[1:]
    if snp[-2].isnumeric() :
        pos = (int) (snp[:-1])
        mt = snp[-1:]
    else:
        pos = (int) (snp[:-4])
        mt = snp[-4:]
    newPos = E_coli_to_M_tuber_SNP.get(pos, "null")
    toAdd = wt + str(newPos) + mt
    return toAdd
    
MEG_6090_sequence = open("MEG_6090_sequence.fasta", "r")
MEG_6090_sequence.readline()
M_tuber_seq = MEG_6090_sequence.readline()
MEG_6090_sequence.close()
MEG_6090_sequence = open("MEG_6090_E_coli_sequence.fasta", "r")
MEG_6090_sequence.readline()
E_coli_seq = MEG_6090_sequence.readline()
MEG_6090_sequence.close()
MEG_6090_SNP = open("MEG_6090_SNP_list.txt", "r")
E_coli = True
for line in MEG_6090_SNP:
    if E_coli:
        if (line.find("/") != -1):
            E_coli = False
        else:
            if line[-1] == "\n":
                line = line[:-1]
            E_coli_MEG_6090_SNP.append(line)
            if (line.find(",") != -1):
                snpList = line.split(",")
                for snp in snpList:
                    makeSNPdict(snp)
            else:
                makeSNPdict(line)
    else:
        if line[-1] == "\n":
            line = line[:-1]
        M_tuberculosis_MEG_6090_SNP.append(line)
MEG_6090_SNP.close()
aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.mode = "global"
aligner.open_gap_score = -.5
aligner.extend_gap_score = -.1
alignments = aligner.align(E_coli_seq[:-1], M_tuber_seq[:-1])
alignString = alignments[0].__str__()
toReturn = open("alignment.txt", "w")
toReturn.write(alignString)
toReturn.write(alignments[0].score.__str__())
toReturn.close()
alignLines = alignString.split("\n")
index = 0
e_pos = 0
m_pos = 0
for i in alignLines[1]:
    e = alignLines[0][index]
    m = alignLines[2][index]
    if e != '-': e_pos += 1
    if m != '-': m_pos += 1
    if i == '|':
        E_coli_to_M_tuber.update({e_pos:m_pos})
    index += 1
for e, m in E_coli_to_M_tuber.items():
    dictEntry = E_coli_pos_SNP.get(e, False)
    if (dictEntry != False):
        E_coli_to_M_tuber_SNP.update({e:m})
for line in E_coli_MEG_6090_SNP:
    toAdd = ""
    if (line.find(",") != -1):
        snpList = line.split(",")
        for snp in snpList:
            toAdd = toAdd + changePos(snp) + ","
        toAdd = toAdd[:-1]
    else:
        toAdd = changePos(line)
    if M_tuberculosis_MEG_6090_SNP.count(toAdd) == 0:
        M_tuberculosis_MEG_6090_SNP.append(toAdd)
toReturn = open("MEG_6090_SNP_list_update.txt", "w")
for line in M_tuberculosis_MEG_6090_SNP:
    toReturn.write(line)
    toReturn.write("\n")
toReturn.close()

        
