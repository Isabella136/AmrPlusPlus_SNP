from SNP_Verification_Tools import disregard, dnaTranslate
from SNP_Verification_Tools import argInfoDict, susceptibleFrameshiftInfoDict, resistantFrameshiftInfoDict, meg_6094InfoDict

def FrameshiftCheck(read, gene, rRna):
    readFrameshiftInfo = ""
    if not(rRna):
        cigarOpCount = read.get_cigar_stats()[0].tolist()
        frameshiftInfo = gene.getFrameshiftInfo()
        #if (insertions - deletions) %3 != 0, disregard unless if gene has frameshift info
        if (((cigarOpCount[1] - cigarOpCount[2]) % 3) != 0) and (frameshiftInfo == None): 
            disregard(gene.getName(), argInfoDict) 
            addRead(gene.getName(), read.query_name, susceptibleFrameshiftInfoDict, "")
            return False
        elif frameshiftInfo != None:
            if gene.getName() == "MEG_6094":
                querySequence = read.query_sequence
                aligned_pairs = read.get_aligned_pairs()[read.cigartuples[0][1] if read.cigartuples[0][0] == 4 else 0:len(querySequence)-read.cigartuples[-1][1] if read.cigartuples[-1][0] == 4 else len(querySequence)]
                shiftCount = 0
                inCodon531 = False
                hasCinsertion = False
                hasDeletionAfterCinsertion = False
                lastBeforeFull = (3 - (aligned_pairs[0][1]) % 3) % 3 - 1            #Last nt index before first full codon
                residue531To534 = ""
                valid = True
                for pair in aligned_pairs:
                    if pair[0] == None:                                             #If deletion
                        if hasCinsertion and not(hasDeletionAfterCinsertion):       #If first deletion since C insertion at codon 531
                            hasDeletionAfterCinsertion = True
                        else:
                            shiftCount -= 1
                    elif pair[1] == None:                                           #If insertion
                        if inCodon531:                                              #If in codon 531
                            if not(hasCinsertion):
                                shiftCount -= 1
                            hasCinsertion = True
                        shiftCount += 1
                    if pair[1] == 1590:
                        inCodon531 = True
                    if pair[1] == (lastBeforeFull + 3):
                        if inCodon531 or (len(residue531To534) >= 1 and len(residue531To534) < 4):
                            residue531To534 += dnaTranslate(querySequence[lastBeforeFull+1:pair[1]+1])
                        lastBeforeFull += 3
                if (shiftCount % 3) == 2:
                    if not(hasCinsertion) or hasDeletionAfterCinsertion:
                        valid = False
                elif (shiftCount % 3) == 1:
                    valid = False
                if not(valid):
                    addRead(gene.getName(), read.query_name, meg_6094InfoDict, "Insertions/deletions that lead to a frameshift for the rest of the reference in residues other than 531")
                    return False
                elif hasCinsertion:
                    if residue531To534 == "SRTR"[0:len(residue531To534)]:
                        addRead(gene.getName(), read.query_name, meg_6094InfoDict, "Has a C insertion at residue 531")
                        return True
            else:
                if (((cigarOpCount[1] - cigarOpCount[2]) % 3) != 0):
                    addRead(gene.getName(), read.query_name, resistantFrameshiftInfoDict, "Insertions/deletions that lead to a frameshift for the rest of the reference")
                    return True
                else:
                    currentPosition = read.reference_start
                    for cigarTuple in read.cigartuples:
                        if (cigarTuple[0] in range(1,3)) and (cigarTuple[1] >= 12):
                            InOrDel = "insertion" if cigarTuple[0] == 1 else "deletion"
                            readFrameshiftInfo = readFrameshiftInfo + str(cigarTuple[1]) + " bp long " + InOrDel + " at nucleotide position " + str(currentPosition + 1) + ","
                        if cigarTuple[0] not in range(4,6):
                            currentPosition += cigarTuple[1]
                    return readFrameshiftInfo
    return None

def addRead(name, queryName, frameshiftInfoDict, additionalInformation = None):
    if name not in frameshiftInfoDict:
        frameshiftInfoDict.update({name:list()})
    if (additionalInformation == None) or (additionalInformation == ''):
        frameshiftInfoDict[name].append(queryName)
    else:
        frameshiftInfoDict[name].append((queryName, additionalInformation))