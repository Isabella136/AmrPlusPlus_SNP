from SNP_Verification_Tools import resistant, disregard
from SNP_Verification_Tools import Gene
from SNP_Verification_Tools import SNP
from SNP_Verification_Processes.FrameshiftCheck import addRead
from SNP_Verification_Tools import argInfoDict, meg_3180InfoDict, resistantFrameshiftInfoDict, meg_6094InfoDict

def NonsenseCheck(read, gene, mapOfInterest, seqOfInterest, frameshiftCheckResult):
    stopLocation = seqOfInterest.find('*') 
    #if nonsense mutations
    if (stopLocation != -1) & ((stopLocation < len(seqOfInterest)-1) | (read.reference_end < gene.ntSeqLength())): 
        return verifyNonsense(read, gene, stopLocation, mapOfInterest, frameshiftCheckResult)
    return None

def verifyNonsense(read, gene, stopLocation, mapOfInterest, frameshiftCheckResult):
    SNPInfo = gene.condensedNonInfo()
    frameshiftInfo = gene.getFrameshiftInfo()
    if (len(SNPInfo) == 0) and ((frameshiftInfo == None) or (gene.getName() == "MEG_6094")):
        if gene.getName() == "MEG_3180":
            if "susceptible" not in meg_3180InfoDict:
                meg_3180InfoDict["susceptible"] = 0
            meg_3180InfoDict["susceptible"] += 1
        elif gene.getName() == "MEG_6094":
            values = list(mapOfInterest.values())
            referenceStopLocation = -1
            i = 0
            for tuple in values:
                if stopLocation in tuple:
                    referenceStopLocation = list(mapOfInterest.keys())[i]
                    break
                i+=1
            if referenceStopLocation == -1:
                raise NameError("mapOfInterest and seqOfInterest don't match")
            addRead(gene.getName(), read.query_name, meg_6094InfoDict, "Nonsense mutation at codon " + str(referenceStopLocation + 1) + " which was NOT previously recorded in literature")
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
                        if frameshiftCheckResult == None: frameshiftCheckResult = ""
                        frameshiftCheckResult = frameshiftCheckResult + "Nonsense mutation at codon " + str(snp[1]) + " which was previously recorded in literature"
                        addRead(gene.getName(), read.query_name, resistantFrameshiftInfoDict, frameshiftCheckResult)
                    res = 1
                    toReturn = True
                    break
        if (frameshiftInfo != None) and (res == 0):
            values = list(mapOfInterest.values())
            referenceStopLocation = -1
            i = 0
            for tuple in values:
                if stopLocation in tuple:
                    referenceStopLocation = list(mapOfInterest.keys())[i]
                    break
                i+=1
            if referenceStopLocation == -1:
                raise NameError("mapOfInterest and seqOfInterest don't match")
            if frameshiftCheckResult == None: frameshiftCheckResult = ""
            frameshiftCheckResult = frameshiftCheckResult + "Nonsense mutation at codon " + str(referenceStopLocation) + " which was NOT previously recorded in literature"
            addRead(gene.getName(), read.query_name, resistantFrameshiftInfoDict, frameshiftCheckResult)
            return True
        resistant(gene.getName(), res, argInfoDict)
        return toReturn
