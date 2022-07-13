from SNP_Verification_Tools import resistant
from SNP_Verification_Tools import Gene
from SNP_Verification_Tools import SNP

def MissenseAndDeletionCheck(gene, name, mapOfInterest, seqOfInterest, mt_and_wt, argInfoDict):
    SNPInfo = gene.condensedRegDelInfo()
    res = 0
    for snp in SNPInfo:
        if (mapOfInterest.get(snp[1]-1,False)) == False:
            continue
        for mt in snp[2]: 
            for queryIndex in tuple(mapOfInterest[snp[1]-1]):
                if queryIndex == '-': continue
                if mt_and_wt:
                    if mt == seqOfInterest[queryIndex]:
                        res = 1
                        break
                else:
                    if mt == seqOfInterest[queryIndex]:
                        res = 1
                    elif snp[0] == seqOfInterest[queryIndex]:
                        res = -1
                        break
            if res == 1: break
            elif res == -1:
                res = 0
                break
        if res == 1: 
            resistant(name, res, argInfoDict)
            return True
    return False