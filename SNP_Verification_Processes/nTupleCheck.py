from SNP_Verification_Tools import resistant
from SNP_Verification_Tools import Gene
from SNP_Verification_Tools import SNP

def nTupleCheck(mapOfInterest, gene, name, seqOfInterest, mt_and_wt, argInfoDict, meg_3180InfoDict):
    SNPInfo = gene.condensedMultInfo()
    if (len(SNPInfo) == 0):
        return False
    else:
        codonNotFullyDeleted = 0    #codons that count for one deletion mutation but aren't fully deleted in read; if >1, snpMult is disregarded        
        for snpMult in SNPInfo:
            codonNotFullyDeleted = 0
            if snpMult[0][2] == "-":        #n-tuple that consists only of mutations
                lastPos = 0                 #if firstPos < beginningPos, all remaining positions in between 
                beginningPos = None         #first and beginning must be deleted
                firstPos = None             #firstPos will stay None if there is a gap in deletions
                delMut = 0
                misMut = 0
                notValid = False
                for mtInfo in snpMult:
                    for pos in mtInfo[1]:
                        deletedOrWt = False
                        if (mapOfInterest.get(pos-1,False)) == False:
                            continue
                        for queryIndex in tuple(mapOfInterest[pos-1]):
                            if queryIndex == "-":
                                deletedOrWt = True
                                delMut += 1
                                if beginningPos == None:
                                    beginningPos = pos
                                    lastPos = pos
                                    firstPos = pos
                                elif pos > lastPos:
                                    if firstPos != beginningPos:
                                        notValid = True
                                        break
                                    if pos-lastPos > 1:
                                        firstPos = None
                                    lastPos = pos
                                elif beginningPos > pos:
                                    if firstPos == None:
                                        notValid = True
                                        break
                                    firstPos = pos
                            else:
                                misMut += 1
                                if mtInfo[0] == queryIndex:
                                    deletedOrWt = True
                        if (delMut > 0) & (misMut > 0):
                            codonNotFullyDeleted += 1
                        if not(deletedOrWt) or notValid:
                            notValid = True
                            break
                    if ((codonNotFullyDeleted > 1) and (firstPos == None)) or notValid or (codonNotFullyDeleted > 2):
                        notValid = True
                        break
                if notValid:
                    continue
                resistant(name, 1, argInfoDict)
                return True

            resBool = True
            wtPresent = False           #can only be true if mt_and_wt is false
            for mtInfo in snpMult:
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
                    if count != len(mtInfo[1]):
                        resBool = False
                        break

                elif mtInfo[2] == "-":
                    delMut = 0
                    misMut = 0
                    currentResBool = False
                    for pos in mtInfo[1]:
                        if (mapOfInterest.get(pos-1,False)) == False:
                            continue
                        for queryIndex in tuple(mapOfInterest[pos-1]):
                            if queryIndex == "-":
                                delMut = 0
                                currentResBool = True
                            else:
                                misMut = 0
                        if codonNotFullyDeleted > 1:
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
                resistant(name, 1, argInfoDict)
                return True
    return False