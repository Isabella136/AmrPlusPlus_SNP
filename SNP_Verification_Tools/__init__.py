aa = {
    "TTT" : 'F',    "TCT" : 'S',    "TAT" : 'Y',    "TGT" : 'C',
    "TTC" : 'F',    "TCC" : 'S',    "TAC" : 'Y',    "TGC" : 'C',
    "TTA" : 'L',    "TCA" : 'S',    "TAA" : '*',    "TGA" : '*',
    "TTG" : 'L',    "TCG" : 'S',    "TAG" : '*',    "TGG" : 'W',
    "CTT" : 'L',    "CCT" : 'P',    "CAT" : 'H',    "CGT" : 'R',
    "CTC" : 'L',    "CCC" : 'P',    "CAC" : 'H',    "CGC" : 'R',
    "CTA" : 'L',    "CCA" : 'P',    "CAA" : 'Q',    "CGA" : 'R',
    "CTG" : 'L',    "CCG" : 'P',    "CAG" : 'Q',    "CGG" : 'R',
    "ATT" : 'I',    "ACT" : 'T',    "AAT" : 'N',    "AGT" : 'S',
    "ATC" : 'I',    "ACC" : 'T',    "AAC" : 'N',    "AGC" : 'S',
    "ATA" : 'I',    "ACA" : 'T',    "AAA" : 'K',    "AGA" : 'R',
    "ATG" : 'M',    "ACG" : 'T',    "AAG" : 'K',    "AGG" : 'R',
    "GTT" : 'V',    "GCT" : 'A',    "GAT" : 'D',    "GGT" : 'G',
    "GTC" : 'V',    "GCC" : 'A',    "GAC" : 'D',    "GGC" : 'G',
    "GTA" : 'V',    "GCA" : 'A',    "GAA" : 'E',    "GGA" : 'G',
    "GTG" : 'V',    "GCG" : 'A',    "GAG" : 'E',    "GGG" : 'G'
}

def dnaTranslate(dna):
    toReturn = ""
    for i in range(0, len(dna), 3):
        if dna[i:i+3].find('N') == -1:
            toReturn += aa[dna[i:i+3]]
        else:
            toReturn += 'N'
    return toReturn
def reverseTranslation(amino_acid):
    for codon, aAcid in aa.items():
        if amino_acid == aAcid:
            return codon
    return ""