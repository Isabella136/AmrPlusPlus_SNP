from SNP_Verification_Tools.Gene import Gene
from SNP_Verification_Tools.SNP import SNP
from SNP_Verification_Tools.InDel import InDel

def MisInDelCheck(read, gene, mapOfInterest, seqOfInterest, config):
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
                    elif config.getboolean('SETTINGS', 'MT_AND_WT'):
                        if mt == seqOfInterest[queryIndex]:
                            res = 1
                            break
                    else:
                        if mt == seqOfInterest[queryIndex]:
                            res = 1
                        elif mtInfo[0] == seqOfInterest[queryIndex]:
                            if ((mapOfInterest.get(mtInfo[1],False)) != False) and (queryIndex in mapOfInterest[mtInfo[1]]):
                                if gene.finalSequence()[mtInfo[1]] != mtInfo[0]:
                                    res = -1
                                    break
                            elif ((mapOfInterest.get(mtInfo[1]-2,False)) != False) and queryIndex in mapOfInterest[mtInfo[1]-2]:
                                if gene.finalSequence()[mtInfo[1]-2] != mtInfo[0]:
                                    res = -1
                                    break
                            else:
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
        for snp in gene.condensedHyperInfo()[0]:
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
        if hasHsMt:
            gene.addDetails(read, "hypersusceptible")

    mtInfoList = gene.condensedMisInDelInfo()
    for mtInfo in mtInfoList:
        if mtInfo[2] == "+":
            count = 0                                                   #must be equal to len(mtInfo[1]) to be considered resistant
            if len(mtInfo[0])  > 1:
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
            if count == len(mtInfo[1]):
                gene.addDetails(read, mtInfo)

        elif mtInfo[2] == "-":
            codonNotFullyDeleted = 0
            deletionCount = 0
            for pos in mtInfo[1]:
                foundADeletion = False
                misMut = 0
                delMut = 0
                if (mapOfInterest.get(pos-1,False)) == False:
                    continue
                remainingResidueIsEqualToOriginal = (False, False)      #1/2M3D2/1M, must be both True or False to be res
                for queryIndex in tuple(mapOfInterest[pos-1]):
                    if (queryIndex == "-") or (queryIndex == None):
                        foundADeletion = True
                        delMut += 1
                    else:
                        misMut += 1
                        if seqOfInterest[queryIndex] == mtInfo[0]:
                            remainingResidueIsEqualToOriginal = (True,False)
                            if ((mapOfInterest.get(pos,False)) != False) and (queryIndex in mapOfInterest[pos]):
                                if gene.finalSequence()[pos] == mtInfo[0]:
                                    remainingResidueIsEqualToOriginal = (True,True)
                            elif ((mapOfInterest.get(pos-2,False)) != False) and queryIndex in mapOfInterest[pos-2]:
                                if gene.finalSequence()[pos-2] == mtInfo[0]:
                                    remainingResidueIsEqualToOriginal = (True,True)
                if (misMut > 0) and (delMut > 0):
                    codonNotFullyDeleted += 1
                    if ((codonNotFullyDeleted % 2) == 1) and (remainingResidueIsEqualToOriginal[0]==remainingResidueIsEqualToOriginal[1]):
                        deletionCount += 1
                elif foundADeletion:
                    deletionCount += 1
            if (deletionCount == 1):
                gene.addDetails(read, mtInfo)

        else:
            missenseCheck(mtInfo)

    nonstop, alongWithNonstop = gene.getNonstopInfo()
    if nonstop and (gene.getCurrentReadNonstopInformation() == None):                #Only MEG_3594 can fullfill this condition
        nonstop = False
        if (read.reference_end == gene.ntSeqLength()):
            if (mapOfInterest.get(gene.ntSeqLength()/3-1, False)) != False:
                for queryIndex in tuple(mapOfInterest[gene.ntSeqLength()/3-1]):
                    #only stop to Ser was found in resistant genes
                    if seqOfInterest[queryIndex] == "S":
                        nonstop = True
                    #regardless of mt_and_wt, presence of stop would make read susceptible
                    elif seqOfInterest[queryIndex] == "*":              
                        nonstop = False
                        break
                if nonstop:
                    if alongWithNonstop == None:
                        gene.addDetails(read, "nonstop")
                    else: 
                        missenseCheck(alongWithNonstop, True)