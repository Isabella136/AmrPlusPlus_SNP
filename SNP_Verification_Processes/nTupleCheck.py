from SNP_Verification_Tools import resistant
from SNP_Verification_Tools import Gene
from SNP_Verification_Tools import SNP
from SNP_Verification_Processes.FrameshiftCheck import addRead
from SNP_Verification_Tools import argInfoDict, meg_3180InfoDict, meg_6094InfoDict, resistantFrameshiftInfoDict, mt_and_wt

def nTupleCheck(read, gene, mapOfInterest, seqOfInterest):
    SNPInfo = gene.condensedMultInfo()
    if (len(SNPInfo) == 0):
        return False
    else:
        codonNotFullyDeleted = 0                #codons that count for one deletion mutation but aren't fully deleted in read; if >1, snpMult is disregarded        
        for snpMult in SNPInfo:
            codonNotFullyDeleted = 0
            if gene.getName() == "MEG_1731":    #NFQ deletion
                lastPos = 0                     #if firstPos < beginningPos, all remaining positions in between 
                beginningPos = None             #first and beginning must be deleted
                firstPos = None
                delMut = 0
                misMut = 0
                notValid = False
                for mtInfo in snpMult:
                    foundDel = False
                    foundWt = False
                    notFound = 0
                    for pos in mtInfo[1]:
                        deletedOrWt = False
                        if (mapOfInterest.get(pos-1,False)) == False:
                            notFound += 1
                            continue
                        for queryIndex in tuple(mapOfInterest[pos-1]):
                            if queryIndex == "-":
                                if not(deletedOrWt):
                                    foundDel = True
                                    delMut += 1
                                deletedOrWt = True
                            else:
                                misMut += 1
                                if mtInfo[0] == queryIndex:
                                    foundDel = False
                                    delMut = 0
                                    deletedOrWt = True
                                    foundWt = True
                        if not(deletedOrWt):
                            notValid = True
                            break
                        elif foundDel and (delMut > 0):
                            if beginningPos == None:
                                beginningPos = pos
                                lastPos = pos
                                firstPos = pos
                            elif pos > lastPos+1:
                                notValid = True
                                break
                            elif (beginningPos > pos):
                                if (firstPos == beginningPos):
                                    firstPos = pos
                                else:
                                    notValid = True
                                    break
                            elif (firstPos != beginningPos) and (pos < firstPos):
                                notValid = True
                                break
                        if (delMut > 0) & (misMut > 0):
                            codonNotFullyDeleted += 1
                    if notValid or (codonNotFullyDeleted > 2) or (notFound > 1) or not(foundDel) or not(foundWt):
                        notValid = True
                        break
                if notValid:
                    continue
                resistant(gene.getName(), 0 if notValid else 1, argInfoDict)
                return True

            resBool = True
            wtPresent = False                   #can only be true if mt_and_wt is false
            for mtInfo in snpMult:
                if mtInfo[2] == "+":
                    count = 0                   #must be equal to len(mtInfo[1]) to be considered resistant
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
                    if count != len(mtInfo[1]):
                        resBool = False
                        break

                elif mtInfo[2] == "-":
                    delMut = 0
                    misMut = 0
                    currentResBool = False
                    if (mapOfInterest.get(mtInfo[1][0]-1,False)) == False:
                        resBool = False
                        break
                    remainingResidueIsEqualToOriginal = (False, False)  #1/2M3D2/1M, must be both True or False to be res
                    for queryIndex in tuple(mapOfInterest[mtInfo[1][0]-1]):
                        if queryIndex == "-":
                            delMut += 1
                            currentResBool = True
                        else:
                            if seqOfInterest[queryIndex] == mtInfo[0]:
                                if remainingResidueIsEqualToOriginal[0]:
                                    remainingResidueIsEqualToOriginal = (True,True)
                                else:
                                    remainingResidueIsEqualToOriginal = (True,False)
                            misMut += 1
                    if (codonNotFullyDeleted > 1) or ((delMut > 0) and (remainingResidueIsEqualToOriginal[0]!=remainingResidueIsEqualToOriginal[1])):
                        resBool = False
                        break
                    resBool = currentResBool
                    if not(resBool):
                        break
                        
                else:
                    if (mapOfInterest.get(mtInfo[1]-1,False)) == False:
                        resBool = False
                        break
                    for mt in mtInfo[2]: 
                        currentResBool = False
                        for queryIndex in tuple(mapOfInterest[mtInfo[1]-1]):
                            if queryIndex == '-': continue
                            if mt == seqOfInterest[queryIndex]:
                                currentResBool = True
                                if mt_and_wt:
                                    break
                            elif (mtInfo[0] == seqOfInterest[queryIndex]) and not(mt_and_wt):
                                wtPresent = True
                                currentResBool = False
                                resBool = False
                                break
                            else:
                                resBool = False
                        if currentResBool: 
                            resBool = currentResBool
                            break
                        if wtPresent:
                            break
                    if not(resBool): break

            if resBool: 
                if (gene.getName() == "MEG_3180"):
                    if "resistant" not in meg_3180InfoDict:
                        meg_3180InfoDict["resistant"] = 0
                    meg_3180InfoDict["resistant"] += 1   
                elif (gene.getName() == "MEG_6094"):
                    addRead(gene.getName(), read.query_name, meg_6094InfoDict, "Has a resistance-conferring n-tuple mutation")
                elif (gene.getFrameshiftInfo() != None):
                    addRead(gene.getName(), read.query_name, resistantFrameshiftInfoDict, "Has a resistance-conferring n-tuple mutation")
                else:
                    resistant(gene.getName(), 1, argInfoDict)
                return True
    return False