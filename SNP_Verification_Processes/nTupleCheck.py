from SNP_Verification_Tools.Gene import Gene
from SNP_Verification_Tools.SNP import SNP
from SNP_Verification_Tools.InDel import InDel

def nTupleCheck(read, gene, mapOfInterest, seqOfInterest, config):
    SNPInfo = gene.condensedMultInfo()
    if (len(SNPInfo) != 0):
        codonNotFullyDeleted = 0                #codons that count for one deletion mutation but aren't fully deleted in read; if >1, snpMult is disregarded        
        for snpMult in SNPInfo:
            codonNotFullyDeleted = 0
            if gene.getName() == "MEG_1731":    #NFQ deletion
                lastPos = 0                     #if firstPos < beginningPos, all remaining positions in between 
                beginningPos = None             #first and beginning must be deleted
                firstPos = None
                notValid = False
                for mtInfo in snpMult[0]:
                    foundDel = False
                    foundWt = False
                    notFound = 0
                    for pos in mtInfo[1]:
                        delMut = 0
                        misMut = 0
                        deletedOrWt = False
                        if (mapOfInterest.get(pos-1,False)) == False:
                            notFound += 1
                            continue
                        for queryIndex in tuple(mapOfInterest[pos-1]):
                            if queryIndex == "-":
                                if not(deletedOrWt):
                                    delMut += 1
                                deletedOrWt = True
                            else:
                                misMut += 1
                                if mtInfo[0] == seqOfInterest[queryIndex]:
                                    if delMut > 0:
                                        delMut = 0
                                    deletedOrWt = True
                                    foundWt = True
                        if not(deletedOrWt):
                            notValid = True
                            break
                        if delMut > 0:
                            foundDel = True
                        if foundDel and (delMut > 0):
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
                            else:
                                lastPos += 1
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
                    if len(mtInfo[0]) > 1:
                        if (mapOfInterest.get(mtInfo[1][0]-1,False)) == False:
                            continue
                        for queryIndex in tuple(mapOfInterest[mtInfo[1][0]-1]):
                            if queryIndex == '-': continue
                            if len(mtInfo[1]) > 1:
                                if mtInfo[0][0] != seqOfInterest[queryIndex]:
                                    continue
                                for startingIndex in range(queryIndex - len(mtInfo[0]), queryIndex + len(mtInfo[0]) + 1, len(mtInfo[0])):
                                    fullInsertionString = True
                                    for i in range(len(mtInfo[0])):
                                        if (len(seqOfInterest) <= (startingIndex + i)) or (0>(startingIndex + i)):
                                            fullInsertionString = False
                                            break
                                        if mtInfo[0][i] != seqOfInterest[startingIndex+i]:
                                            fullInsertionString = False
                                            break
                                    if fullInsertionString:
                                        count += 1
                            else:
                                fullInsertionString = True
                                if seqOfInterest[queryIndex] in mtInfo[0]:
                                    for i in range(len(mtInfo[0])):
                                        if mtInfo[0][i] == seqOfInterest[queryIndex]:
                                            queryIndex -= i
                                            break
                                for i in range(len(mtInfo[0])):
                                    if (len(seqOfInterest) <= (queryIndex + i)) or (0>(queryIndex + i)):
                                        fullInsertionString = False
                                        break
                                    if mtInfo[0][i] != seqOfInterest[queryIndex+i]:
                                        fullInsertionString = False
                                        break
                                if fullInsertionString:
                                    count += 1
                    else:
                        if (mapOfInterest.get(mtInfo[1][0]-1,False)) == False:
                            continue
                        for queryIndex in tuple(mapOfInterest[mtInfo[1][0]-1]):
                            if queryIndex == '-': continue
                            if mtInfo[0][0] == seqOfInterest[queryIndex]:          #should be first of inserted residue
                                count += 1
                                i = queryIndex - 1
                                while (i >= 0) and (seqOfInterest[i] == mtInfo[0][0]):
                                    count += 1
                                    i -= 1
                                i = queryIndex + 1
                                while (i < len(seqOfInterest)) and (seqOfInterest[i] == mtInfo[0][0]):
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
                                if ((mapOfInterest.get(mtInfo[1][0],False)) != False) and (queryIndex in mapOfInterest[mtInfo[1][0]]):
                                    if gene.finalSequence()[mtInfo[1][0]] == mtInfo[0]:
                                        remainingResidueIsEqualToOriginal = (True,True)
                                elif ((mapOfInterest.get(mtInfo[1][0]-2,False)) != False) and queryIndex in mapOfInterest[mtInfo[1][0]-2]:
                                    if gene.finalSequence()[mtInfo[1][0]-2] == mtInfo[0]:
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
                                if config.getboolean('SETTINGS', 'MT_AND_WT'):
                                    break
                            elif (mtInfo[0] == seqOfInterest[queryIndex]) and not(config.getboolean('SETTINGS', 'MT_AND_WT')):
                                if ((mapOfInterest.get(mtInfo[1],False)) != False) and (queryIndex in mapOfInterest[mtInfo[1]]):
                                    if gene.finalSequence()[mtInfo[1]] != mtInfo[0]:
                                        wtPresent = True
                                        currentResBool = False
                                        resBool = False
                                        break
                                elif ((mapOfInterest.get(mtInfo[1]-2,False)) != False) and queryIndex in mapOfInterest[mtInfo[1]-2]:
                                    if gene.finalSequence()[mtInfo[1]-2] != mtInfo[0]:
                                        wtPresent = True
                                        currentResBool = False
                                        resBool = False
                                        break
                                else:
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