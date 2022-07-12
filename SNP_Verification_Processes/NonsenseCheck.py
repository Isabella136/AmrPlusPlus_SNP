from SNP_Verification_Processes import disregard, resistant
from SNP_Verification_Tools import Gene
from SNP_Verification_Tools import SNP

def NonsenseCheck(seqOfInterest, mapOfInterest, gene, read, name, argInfoDict):
    stopLocation = seqOfInterest.find('*') 
    #if nonsense mutations
    if (stopLocation != -1) & ((stopLocation < len(seqOfInterest)-1) | (read.reference_end < gene.ntSeqLength())): 
        return verifyNonsense(stopLocation, mapOfInterest, gene, name, argInfoDict)
    return None

def verifyNonsense(stopLocation, aa_alignment_map, gene, name, argInfoDict):
    SNPInfo = gene.condensedNonInfo()
    if (len(SNPInfo) == 0):
        disregard(name, argInfoDict)
        return False
    else:
        toReturn = False
        res = 0
        for snp in SNPInfo:
            stopLocTemp = aa_alignment_map.get(snp[1]-1,False)
            if stopLocTemp != False:
                if stopLocTemp[0] == stopLocation:
                    res = 1
                    toReturn = True
                    break
        resistant(name, res, argInfoDict)
        return toReturn
