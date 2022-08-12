from SNP_Verification_Tools.Gene import Gene
from SNP_Verification_Tools.SNP import SNP
from SNP_Verification_Tools.InDel import InDel
from SNP_Verification_Tools import mt_and_wt

def nTupleCheck(read, gene, mapOfInterest, seqOfInterest):
    SNPInfo = gene.condensedMultInfo()
    if (len(SNPInfo) != 0):
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
                for mtInfo in snpMult[0]:
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
                gene.addDetails(read, snpMult)

            resBool = True
            wtPresent = False                   #can only be true if mt_and_wt is false
            for mtInfo in snpMult[0]:
                if mtInfo[2] == "+":
                    count = 0                                                   #must be equal to len(mtInfo[1]) to be considered resistant
                    if len(mtInfo[0] > 1):
                        if (mapOfInterest.get(mtInfo[1][0]-1,False)) == False:
                            continue
                        for queryIndex in tuple(mapOfInterest[mtInfo[1][0]-1]):
                            if mtInfo[0][i] != seqOfInterest[queryIndex]:
                                continue
                            for startingIndex in range(2 * mtInfo[1][0] - mtInfo[1][1],mtInfo[1][1]+1, mtInfo[1][1] - mtInfo[1][0]):
                                fullInsertionString = True
                                for i in range(len(mtInfo[0])):
                                    if mtInfo[0][i] != seqOfInterest[startingIndex+i]:
                                        fullInsertionString = False
                                        break
                                if fullInsertionString:
                                    count += 1
                    else:
                        if (mapOfInterest.get(mtInfo[1][0]-1,False)) == False:
                            continue
                        for queryIndex in tuple(mapOfInterest[mtInfo[1][0]-1]):
                            if mtInfo[0] == seqOfInterest[queryIndex]:          #should be first of inserted residue
                                count += 1
                                i = queryIndex - 1
                                while seqOfInterest[i] == mtInfo[0]:
                                    count += 1
                                    i -= 1
                                i = queryIndex + 1
                                while seqOfInterest[i] == mtInfo[0]:
                                    count += 1
                                    i += 1
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
                                remainingResidueIsEqualToOriginal = (True,False)
                                if (mapOfInterest.get(mtInfo[1][0],False)) != False:
                                    if queryIndex in mapOfInterest[mtInfo[1][0]]:
                                        remainingResidueIsEqualToOriginal = (True,True)
                                elif (mapOfInterest.get(mtInfo[1][0]-2,False)) != False:
                                    if queryIndex in mapOfInterest[mtInfo[1][0]-2]:
                                        remainingResidueIsEqualToOriginal = (True,True)
                            misMut += 1
                    if (delMut > 0) & (misMut > 0):
                        codonNotFullyDeleted += 1
                    if (codonNotFullyDeleted > 1) or ((delMut > 0) and (remainingResidueIsEqualToOriginal[0]!=remainingResidueIsEqualToOriginal[1])):
                        currentResBool = False
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
                gene.addDetails(read, snpMult)