mycoplasmataceaeList = ["MEG_3185", "MEG_3240", "MEG_5328", "MEG_5330", "MEG_5331"]

mt_and_wt = True    #used in case of insertion leading to presence of both mt and wt; if true, mark as resistant; if false, mark as susceptible
detailed = False    #determines whether a more detailed output will be given

geneDict = {}


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

def resistant(name, increment, argInfoDict):
    argInfo = argInfoDict.get(name, False)
    if (argInfo == False):
        argInfoDict.update({name:(increment, 1)})
    else:
        temp = list(argInfo)
        temp[0] += increment
        temp[1] += 1
        argInfo = tuple(temp)
        argInfoDict.update({name:argInfo})

def disregard(name, argInfoDict):
    resistant(name, 0, argInfoDict)

def dnaTranslate(dna, name):
    toReturn = ""
    for i in range(0, len(dna), 3):
        if dna[i:i+3].find('N') == -1:
            if (name in mycoplasmataceaeList) and (dna[i:i+3] == "TGA"):
                toReturn += 'W'
            else:
                toReturn += aa[dna[i:i+3]]
        else:
            toReturn += 'N'
    return toReturn

def reverseTranslation(amino_acid):
    if amino_acid == '*':
        return "TAA"
    for codon, aAcid in aa.items():
        if amino_acid == aAcid:
            return codon
    return ""