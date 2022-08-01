from SNP_Verification_Tools import resistant
from SNP_Verification_Tools import Gene
from SNP_Verification_Tools import SNP
from SNP_Verification_Tools import InDel

def MisInDelCheck(gene, name, mapOfInterest, seqOfInterest, mt_and_wt, argInfoDict, meg_3180InfoDict):
    hasHsMt = False
    if gene.getName() == "MEG_3180":
        for snp in gene.condensedHyperInfo():
            if (mapOfInterest.get(snp[1]-1,False)) == False:
                break
            for mt in snp[2]:
                for queryIndex in tuple(mapOfInterest[snp[1]-1]):
                    #Regardless of whether the wild-type is also present due to InDels
                    if mt == seqOfInterest[queryIndex]:
                        hasHsMt = True
                        break
                    hasHsMt = False
                if hasHsMt == True:
                    break
            if not(hasHsMt):
                break

    mtInfoList = gene.condensedMisInDelInfo()
    res = 0
    for mtInfo in mtInfoList:
        if mtInfo[2] == "+":
            count = 0       #must be equal to len(mtInfo[1]) to be considered resistant
            for pos in mtInfo[1]:
                if (mapOfInterest.get(pos-1,False)) == False:
                    continue
                for queryIndex in tuple(mapOfInterest[pos-1]):
                    fullInsertionString = True
                    for i in range(len(mtInfo[0])):
                        if mtInfo[0][i] != seqOfInterest[queryIndex+i]:
                            fullInsertionString = False
                            break
                    if fullInsertionString:
                        count += 1
                        break
            if count == len(mtInfo[1]):
                resistant(name, 1, argInfoDict)
                return True

        elif mtInfo[2] == "-":
            for pos in mtInfo[1]:
                if (mapOfInterest.get(pos-1,False)) == False:
                    continue
                for queryIndex in tuple(mapOfInterest[pos-1]):
                    if queryIndex == "-":
                        resistant(name, 1, argInfoDict)
                        return True

        else:
            resMtCount = 0  #will be used for MEG_3180 if hypersusceptible double mutation is present
            if (mapOfInterest.get(mtInfo[1]-1,False)) == False:
                continue
            for mt in mtInfo[2]: 
                for queryIndex in tuple(mapOfInterest[mtInfo[1]-1]):
                    if queryIndex == '-': continue
                    elif mt_and_wt:
                        if mt == seqOfInterest[queryIndex]:
                            res = 1
                            break
                    else:
                        if mt == seqOfInterest[queryIndex]:
                            res = 1
                        elif mtInfo[0] == seqOfInterest[queryIndex]:
                            res = -1
                            break
                if hasHsMt:
                    if res == 1:
                        resMtCount+=1
                    res = 0
                elif res == 1: break
                elif res == -1:
                    res = 0
                    break
            if res == 1: 
                resistant(name, res, argInfoDict)
                return True

    if gene.getName() == "MEG_3180":
        if not(hasHsMt) and (resMtCount > 0):
            if "resistant" not in meg_3180InfoDict:
                meg_3180InfoDict["resistant"] = 0
            meg_3180InfoDict["resistant"] += 1    
            return True
        elif hasHsMt:
            if resMtCount > 0:
                if resMtCount not in meg_3180InfoDict:
                    meg_3180InfoDict[resMtCount] = 0
                meg_3180InfoDict[resMtCount] += 1
                return True
            else:
                if "susceptible" not in meg_3180InfoDict:
                    meg_3180InfoDict["susceptible"] = 0
                meg_3180InfoDict["susceptible"] += 1
                return True
    return False