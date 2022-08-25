from SNP_Verification_Tools import dnaTranslate

def FrameshiftCheck(read, gene, rRna):
    removeFromLongFrameshiftCheck = [None]
    toReturn = []
    if not(rRna):
        cigarOpCount = read.get_cigar_stats()[0].tolist()
        frameshiftInfo = gene.getFrameshiftInfo()
        #if (insertions - deletions) %3 != 0, susceptible unless if gene has frameshift info
        if (((cigarOpCount[1] - cigarOpCount[2]) % 3) != 0) and (gene.getGeneTag() == 'N') and (gene.getName() != "MEG_6142"): 
            gene.addDetails(read, 'FS till end')
            return False                                                            #Should not continue
        elif frameshiftInfo != None:
            if gene.getName() == "MEG_6142":                                        #Gene has two special cases where frameshift were found to be allowed in the literature
                pass
            if gene.getGeneTag() == 'S':                                            #MEG_6094 can have a C insertion at res 531 if followed by frameshift suppression during translation
                if not(MEG_6094Check(read, gene, removeFromLongFrameshiftCheck)): return False
            elif gene.getGeneTag() == 'F':
                if (((cigarOpCount[1] - cigarOpCount[2]) % 3) != 0):
                    gene.addDetails(read, 'FS till end')
                    extendedIndelCheck(read, gene)
                    return True
        longFrameshiftCheck(read, gene, removeFromLongFrameshiftCheck[0])
    extendedIndelCheck(read, gene)
    return True

def longFrameshiftCheck(read, gene, removeFromLongFrameshiftCheck):
    shiftcount = 0
    longFS = 0
    fsLength = 0
    for cigarTuple in read.cigartuples:
        if (cigarTuple[0] == 1):
            shiftcount += cigarTuple[1]
        elif (cigarTuple[0] == 2):
            shiftcount -= cigarTuple[1]
        elif (cigarTuple[0] == 0):
            if (shiftcount % 3) != 0:
                fsLength += cigarTuple[1]
            elif (shiftcount % 3) == 0:
                if fsLength >= 12:
                    longFS += 1
                fsLength = 0
    if (longFS > 0) and (gene.mustSuppressFrameshift()) or (removeFromLongFrameshiftCheck == True):
        longFS -= 1
    if (longFS > 0):
        gene.addDetails(read, "12+bp frameshift: " + str(longFS))
            

def extendedIndelCheck(read, gene):
    indel = 0
    for cigarTuple in read.cigartuples:
        if (cigarTuple[0] in range(1,3)) and (cigarTuple[1] >= 12):
            indel += 1
    if indel > 0:
        gene.addDetails(read, "12+bp indel: " + str(indel))

def MEG_6142Check(read, gene):
    querySequence = read.query_sequence
    startIndex = read.cigartuples[0][1] if read.cigartuples[0][0] == 4 else 0
    endIndex = len(querySequence)-read.cigartuples[-1][1] if read.cigartuples[-1][0] == 4 else len(querySequence)
    aligned_pairs = read.get_aligned_pairs()[startIndex:endIndex]       #Removes soft-clipping
    shiftCount = 0                                                      #If pos, more ins; if neg. more del
    queryIndex = -1
    stopCodon = ["TAA", "TGA", "TAG"]
    skipTill291 = False
    restartAt291 = False
    for pair in aligned_pairs:
        if skipTill291:
            if restartAt291:
                if pair[0] == None:                                     #If deletion
                    shiftCount -= 1
                elif pair[1] == None:
                    shiftCount += 1
            if pair[0] != None:
                queryIndex += 1
            if pair[1] == 276:
                if querySequence[queryIndex+1:queryIndex+7] == "AGGAGG":
                    restartAt291 = True
        else:
            if pair[0] == None:                                         #If deletion
                shiftCount -= 1
            else:
                queryIndex += 1
                if pair[1] == None:                                     #If insertion
                    shiftCount += 1
            if pair[1] == 75:
                if (querySequence[queryIndex+1:queryIndex+4] in stopCodon) and ((shiftCount%3)==2):
                    skipTill292 = True
                    shiftCount = 0


def MEG_6094Check(read, gene, removeFromLongFrameshiftCheck):
    querySequence = read.query_sequence
    startIndex = read.cigartuples[0][1] if read.cigartuples[0][0] == 4 else 0
    endIndex = len(querySequence)-read.cigartuples[-1][1] if read.cigartuples[-1][0] == 4 else len(querySequence)
    aligned_pairs = read.get_aligned_pairs()[startIndex:endIndex]       #Removes soft-clipping
    shiftCount = 0                                                      #If pos, more ins; if neg. more del
    inCodon531 = False
    hasCinsertion = False
    insertionCountAfterCinsertion = 0
    deletionCountAfterCinsertion = 0
    lastBeforeFull = (3 - (aligned_pairs[0][1]) % 3) % 3 - 1            #Last nt index before first full codon
    residue531To536 = ""
    residue531To534 = ""
    valid = True
    refIndex = aligned_pairs[0][1]
    insertionPosition = None
    temp = removeFromLongFrameshiftCheck[0]
    for pair in aligned_pairs:
        refIndex += 1
        if pair[0] == None:                                             #If deletion
            if hasCinsertion:
                deletionCountAfterCinsertion += 1
                if (temp == None) and ((refIndex - insertionPosition) != 0):
                    temp = (refIndex - insertionPosition)>= 12
            shiftCount -= 1
        elif pair[1] == None:                                           #If insertion
            refIndex -= 1
            if hasCinsertion:
                insertionCountAfterCinsertion += 1
                if (temp == None) and ((refIndex - insertionPosition) != 0):
                    temp = (refIndex - insertionPosition) >= 12
            if inCodon531 and ((shiftCount%3)==0):                      #If in codon 531 and not in frameshift
                if not(hasCinsertion):
                    shiftCount -= 1
                    hasCinsertion = True
                    insertionPosition = refIndex + 1
            shiftCount += 1
        if (insertionCountAfterCinsertion > 0) and (deletionCountAfterCinsertion > 0):
            insertionCountAfterCinsertion -= 1                          #Never will have both over 1
            deletionCountAfterCinsertion -= 1
        if pair[1] == 1590:
            inCodon531 = True
        if pair[1] == 1593:
            inCodon531 = False
        if pair[0] == (lastBeforeFull + 3):
            if inCodon531 or (len(residue531To534) >= 1 and len(residue531To534) < 4):
                residue531To534 += dnaTranslate(querySequence[lastBeforeFull+1:pair[0]+1], gene.getName())
            if inCodon531 or (len(residue531To536) >= 1 and len(residue531To536) < 6):
                residue531To536 += dnaTranslate(querySequence[lastBeforeFull+1:pair[0]+1], gene.getName())
            lastBeforeFull += 3
    removeFromLongFrameshiftCheck[0] = temp
    twoInsertionsAfter = ((insertionCountAfterCinsertion%3) == 2)
    deletionAfter = ((deletionCountAfterCinsertion%3) == 1)
    if (shiftCount % 3) == 2:
        if not(hasCinsertion):
            valid = False
    elif (shiftCount % 3) == 1:
        valid = False
    if not(valid):
        gene.addDetails(read, 'FS till end')
        return False                                                    #Should not continue
    elif hasCinsertion:                                                 #shiftCount%3 == 0 at residue 531
        if (residue531To534 == "SRTR"[0:len(residue531To534)]):         #Must have those residues due to insertion
            if (deletionAfter or twoInsertionsAfter):
                gene.addDetails(read, 'C insert followed by del/ins')   #In that case, must not have FS suppression                                                
            elif (residue531To536 == "SRTRPR"[0:len(residue531To536)]):                                                       
                gene.addDetails(read, 'Suppressible C insert')
            else:
                gene.addDetails(read, 'C insert + not SRTRPR')          #Suppression can't happen if no PR
                gene.addDetails(read, 'FS till end')
                return False
            return True 
        elif (shiftCount % 3) == 0:                                     #Because SRTR is not present, can't make prediction on resistance
            gene.addDetails(read, 'C insert + not SRTR')
            gene.addDetails(read, 'FS till end')
            return False
    return True