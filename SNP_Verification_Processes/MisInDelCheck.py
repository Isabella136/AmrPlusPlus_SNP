from SNP_Verification_Tools.Gene import Gene
from SNP_Verification_Tools.SNP import SNP
from SNP_Verification_Tools.InDel import InDel
from SNP_Verification_Tools import mt_and_wt

def MisInDelCheck(read, gene, mapOfInterest, seqOfInterest):
    def missenseCheck(mtInfo, nonstop = False):
        res = 0
        if (mapOfInterest.get(mtInfo[1]-1,False)) != False:
            for mt in mtInfo[2]:
                firstAndLastTupleIndex = [0,len(mapOfInterest[mtInfo[1]-1]) - 1]
                tupleIndex = -1
                for queryIndex in tuple(mapOfInterest[mtInfo[1]-1]):
                    tupleIndex += 1
                    if gene.rRna():
                        if tupleIndex not in firstAndLastTupleIndex:
                            continue
                    if (queryIndex == '-') or (queryIndex == None): continue
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
                if res == 1: break
                elif res == -1:
                    res = 0
                    break
            if res == 1: 
                if nonstop:
                    gene.addDetails(read, "nonstop")
                else:
                    gene.addDetails(read, mtInfo)

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
                        gene.addDetails(read, "hypersusceptible")
                        break
                if hasHsMt == True:
                    break
            if not(hasHsMt):
                break

    mtInfoList = gene.condensedMisInDelInfo()
    for mtInfo in mtInfoList:
        if mtInfo[2] == "+":
            count = 0                                               #must be equal to len(mtInfo[1]) to be considered resistant
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
                gene.addDetails(read, mtInfo)

        elif mtInfo[2] == "-":
            for pos in mtInfo[1]:
                if (mapOfInterest.get(pos-1,False)) == False:
                    continue
                foundADeletion = False
                remainingResidueIsEqualToOriginal = (False, False)  #1/2M3D2/1M, must be both True or False to be res
                for queryIndex in tuple(mapOfInterest[pos-1]):
                    if (queryIndex == "-") or (queryIndex == None):
                        foundADeletion = True
                    else:
                        if seqOfInterest[queryIndex] == mtInfo[0]:
                            if remainingResidueIsEqualToOriginal[0]:
                                remainingResidueIsEqualToOriginal = (True,True)
                            else:
                                remainingResidueIsEqualToOriginal = (True,False)
                if foundADeletion and (remainingResidueIsEqualToOriginal[0]==remainingResidueIsEqualToOriginal[1]):
                    gene.addDetails(read, mtInfo)

        else:
            resMtCount = missenseCheck(mtInfo)

    nonstop, alongWithNonstop = gene.getNonstopInfo()
    if nonstop:
        if (read.reference_end == gene.ntSeqLength()):
            if (mapOfInterest.get(gene.ntSeqLength()/3-1, False)) != False:
                for queryIndex in tuple(mapOfInterest[gene.ntSeqLength()/3-1]):
                    if seqOfInterest[queryIndex] == "*":
                        nonstop = False
                        break
                if nonstop:
                    if alongWithNonstop == None:
                        gene.addDetails(read, "nonstop")
                    else: 
                        missenseCheck(alongWithNonstop, True)