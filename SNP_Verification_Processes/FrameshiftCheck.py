from SNP_Verification_Tools import dnaTranslate

def FrameshiftCheck(read, gene, rRna):
    toReturn = []
    if not(rRna):
        cigarOpCount = read.get_cigar_stats()[0].tolist()
        frameshiftInfo = gene.getFrameshiftInfo()
        #if (insertions - deletions) %3 != 0, susceptible unless if gene has frameshift info
        if (((cigarOpCount[1] - cigarOpCount[2]) % 3) != 0) and (frameshiftInfo == None): 
            gene.addDetails(read, 'FS till end')
            return False                                                            #Should not continue
        elif frameshiftInfo != None:
            if gene.getGeneTag() == 'S':                                            #MEG_6094 can have a C insertion at res 531 if followed by frameshift suppression during translation
                if not(MEG_6094Check(read, gene)): return False
            else:
                if (((cigarOpCount[1] - cigarOpCount[2]) % 3) != 0):
                    gene.addDetails(read, 'FS till end')
                    extendedIndelCheck(read, gene)
                    return True
        longFrameshiftCheck(read, gene)
    extendedIndelCheck(read, gene)
    return True

def longFrameshiftCheck(read, gene):
    shiftcount = 0
    longFS = 0
    for cigarTuple in read.cigartuples:
        if (cigarTuple[0] == 1):
            shiftcount += cigarTuple[1]
        elif (cigarTuple[0] == 2):
            shiftcount -= cigarTuple[1]
        elif (cigarTuple[0] == 0):
            if ((shiftcount % 3) != 0) and (cigarTuple[1] >= 12):
                longFS += 1
    if longFS > 0:
        gene.addDetails(read, "12+fs: " + str(longFS))
            

def extendedIndelCheck(read, gene):
    indel = 0
    for cigarTuple in read.cigartuples:
        if (cigarTuple[0] in range(1,3)) and (cigarTuple[1] >= 12):
            indel += 1
    if indel > 0:
        gene.addDetails(read, "12+indel: " + str(indel))

def MEG_6094Check(read, gene):
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
    for pair in aligned_pairs:
        if pair[0] == None:                                             #If deletion
            if hasCinsertion:
                deletionCountAfterCinsertion += 1
            shiftCount -= 1
        elif pair[1] == None:                                           #If insertion
            if hasCinsertion:
                insertionCountAfterCinsertion += 1
            if inCodon531 and ((shiftCount%3)==0):                      #If in codon 531 and not in frameshift
                if not(hasCinsertion):
                    shiftCount -= 1
                    hasCinsertion = True
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
                gene.addDetails(read, 'C insert + del/ins')             #In that case, must not have FS suppression                                                
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