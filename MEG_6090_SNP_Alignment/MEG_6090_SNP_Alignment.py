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

from Bio import Align
from Bio import pairwise2
from Bio.Align import substitution_matrices

#0-1: MEG_6090
#2-3: MEG_3241
#4-5: MEG_3243
#6-7: MEG_3246
#8-9: MEG_5779
#even: SNP seq; odd: actual seq

All_Seq = ["", "", "", "", "", "", "", "", "", ""]
All_SNP = [[], [], [], [], [], [], [], [], [], []]

#0: MEG_6090
#1: MEG_3241
#2: MEG_3243
#3: MEG_3246
#4: MEG_5779

posSNP = [{}, {}, {}, {}, {}]
SNP_seq_to_actual = [{}, {}, {}, {}, {}]
SNP_seq_to_actual_SNP = [{}, {}, {}, {}, {}]
MEG = ["MEG_6090", "MEG_3241", "MEG_3243", "MEG_3246", "MEG_5779"]

def makeSNPdict(snp, i):
    pos = 0
    wt = ""
    wt = snp[0]
    snp = snp[1:]
    if snp[-2].isnumeric() :
        pos = (int) (snp[:-1])
    else:
        pos = (int) (snp[:-4])
    dictEntry = posSNP[i].get(pos, False)
    if (dictEntry == False):
        posSNP[i].update({pos:wt})
    elif (dictEntry.find(wt) == -1):
        posSNP[i][pos] += wt

def changePos(snp, i):
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
    newPos = SNP_seq_to_actual_SNP[i].get(pos, "null")
    toAdd = wt + str(newPos) + mt
    return toAdd
  
for i in range(0,5):
    MEG_sequence = open(MEG[i] + "_actual_sequence.fasta", "r")
    MEG_sequence.readline()
    All_Seq[2*i+1] = MEG_sequence.readline()
    MEG_sequence.close()
    MEG_sequence = open(MEG[i] + "_SNP_sequence.fasta", "r")
    MEG_sequence.readline()
    All_Seq[2*i] = MEG_sequence.readline()
    MEG_sequence.close()
    MEG_SNP = open(MEG[i] + "_SNP_list.txt", "r")
    E_coli = True
    for line in MEG_SNP:
        if E_coli:
            if (line.find("/") != -1):
                E_coli = False
            else:
                if line[-1] == "\n":
                    line = line[:-1]
                All_SNP[2*i].append(line)
                if (line.find(",") != -1):
                    snpList = line.split(",")
                    for snp in snpList:
                        makeSNPdict(snp, i)
                else:
                    makeSNPdict(line, i)
        else:
            if line[-1] == "\n":
                line = line[:-1]
            All_SNP[2*i+1].append(line)
    MEG_SNP.close()
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.mode = "global"
    aligner.open_gap_score = -.5
    aligner.extend_gap_score = -.1
    alignments = aligner.align(All_Seq[2*i][:-1], All_Seq[2*i+1][:-1])
    alignString = alignments[0].__str__()
    toReturn = open(MEG[i] + "alignment.txt", "w")
    toReturn.write(alignString)
    toReturn.write(alignments[0].score.__str__())
    toReturn.close()
    alignLines = alignString.split("\n")
    index = 0
    e_pos = 0
    m_pos = 0
    for k in alignLines[1]:
        e = alignLines[0][index]
        m = alignLines[2][index]
        if e != '-': e_pos += 1
        if m != '-': m_pos += 1
        if k == '|':
            SNP_seq_to_actual[i].update({e_pos:m_pos})
        index += 1
    for e, m in SNP_seq_to_actual[i].items():
        dictEntry = posSNP[i].get(e, False)
        if (dictEntry != False):
            SNP_seq_to_actual_SNP[i].update({e:m})
    for line in All_SNP[2*i]:
        toAdd = ""
        if (line.find(",") != -1):
            snpList = line.split(",")
            for snp in snpList:
                toAdd = toAdd + changePos(snp, i) + ","
            toAdd = toAdd[:-1]
        else:
            toAdd = changePos(line, i)
        if All_SNP[2*i+1].count(toAdd) == 0:
            All_SNP[2*i+1].append(toAdd)
    toReturn = open(MEG[i] + "_SNP_list_update.txt", "w")
    for line in All_SNP[2*i+1]:
        toReturn.write(line)
        toReturn.write("\n")
    toReturn.close()

        
