from SNP_Verification_Tools import resistant, disregard
from SNP_Verification_Tools import Gene
from SNP_Verification_Tools import SNP
from SNP_Verification_Processes.FrameshiftCheck import addRead
from SNP_Verification import argInfoDict, meg_3180InfoDict, resistantFrameshiftInfoDict

def NonsenseCheck(read, gene, mapOfInterest, seqOfInterest, frameshiftCheckResult):
    stopLocation = seqOfInterest.find('*') 
    #if nonsense mutations
    if (stopLocation != -1) & ((stopLocation < len(seqOfInterest)-1) | (read.reference_end < gene.ntSeqLength())): 
        return verifyNonsense(read, gene, stopLocation, mapOfInterest, frameshiftCheckResult)
    return None

def verifyNonsense(read, gene, stopLocation, mapOfInterest, frameshiftCheckResult):
    SNPInfo = gene.condensedNonInfo()
    frameshiftInfo = gene.getFrameshiftInfo()
    if (len(SNPInfo) == 0) and (frameshiftInfo == None):
        if gene.getName() == "MEG_3180":
            if "susceptible" not in meg_3180InfoDict:
                meg_3180InfoDict["susceptible"] = 0
            meg_3180InfoDict["susceptible"] += 1
        else:
            disregard(gene.getName(), argInfoDict)
        return False
    else:
        toReturn = False
        res = 0
        for snp in SNPInfo:
            stopLocTemp = mapOfInterest.get(snp[1]-1,False)
            if stopLocTemp != False:
                if stopLocTemp[0] == stopLocation:
                    if gene.getName() == "MEG_3180":
                        if "resistant" not in meg_3180InfoDict:
                            meg_3180InfoDict["resistant"] = 0
                        meg_3180InfoDict["resistant"] += 1
                        return True
                    elif frameshiftInfo != None:
                        frameshiftCheckResult = frameshiftCheckResult + "Nonsense mutation at codon " + str(stopLocation + 1) + " which was previously recorded in literature"
                        addRead(gene.getName(), read.query_name, resistantFrameshiftInfoDict, frameshiftCheckResult)
                    res = 1
                    toReturn = True
                    break
        if (frameshiftInfo != None) and (res == 0):
            frameshiftCheckResult = frameshiftCheckResult + "Nonsense mutation at codon " + str(stopLocation + 1) + " which was NOT previously recorded in literature"
            addRead(gene.getName(), read.query_name, resistantFrameshiftInfoDict, frameshiftCheckResult)
        resistant(gene.getName(), res, argInfoDict)
        return toReturn
