from SNP_Verification_Tools import dnaTranslate

def FrameshiftCheck(read, gene, rRna):
    toReturn = []
    if not(rRna):
        cigarOpCount = read.get_cigar_stats()[0].tolist()
        frameshiftInfo = gene.getFrameshiftInfo()
        #if (insertions - deletions) %3 != 0, susceptible unless if gene has frameshift info
        if (((cigarOpCount[1] - cigarOpCount[2]) % 3) != 0) and (frameshiftInfo == None): 
            gene.addDetails(read, 'FS end')
            return False                                                            #Should not continue
        elif frameshiftInfo != None:
            if gene.getGeneTag() == 'S':                                            #MEG_6094 can have a C insertion at res 531 if followed by frameshift suppression during translation
                if not(MEG_6094Check(read, gene)): return False
            else:
                if (((cigarOpCount[1] - cigarOpCount[2]) % 3) != 0):
                    gene.addDetails(read, 'FS end')
                    extendedIndelCheck(read, gene)
                    return True
    extendedIndelCheck(read, gene)
    longFrameshiftCheck(read, gene)
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
    hasInsertionAfterCinsertion = False
    hasDeletionAfterCinsertion = False
    lastBeforeFull = (3 - (aligned_pairs[0][1]) % 3) % 3 - 1            #Last nt index before first full codon
    residue531To534 = ""
    valid = True
    for pair in aligned_pairs:
        if pair[0] == None:                                             #If deletion
            if hasCinsertion and not(hasDeletionAfterCinsertion):       #If first deletion since C insertion at codon 531
                hasDeletionAfterCinsertion = True
            shiftCount -= 1
        elif pair[1] == None:                                           #If insertion
            if hasCinsertion and not(hasInsertionAfterCinsertion):
                hasInsertionAfterCinsertion = True
            if inCodon531 and ((shiftCount%3)!=0):                      #If in codon 531 and not in frameshift
                if not(hasCinsertion):
                    shiftCount -= 1
                hasCinsertion = True
            shiftCount += 1
        if hasInsertionAfterCinsertion and hasDeletionAfterCinsertion:
            hasDeletionAfterCinsertion = False
            hasInsertionAfterCinsertion = False
        if pair[1] == 1590:
            inCodon531 = True
        if pair[1] == 1593:
            inCodon531 = False
        if pair[0] == (lastBeforeFull + 3):
            if inCodon531 or (len(residue531To534) >= 0 and len(residue531To534) < 4):
                residue531To534 += dnaTranslate(querySequence[lastBeforeFull+1:pair[1]+1])
            lastBeforeFull += 3
    if (shiftCount % 3) == 2:
        if not(hasCinsertion):
            valid = False
    elif (shiftCount % 3) == 1:
        valid = False
    if not(valid):
        gene.addDetails(read, 'FS end')
        return False                                                    #Should not continue
    elif hasCinsertion:
        if ((shiftCount % 3) == 0) and hasDeletionAfterCinsertion:      #Will only be true if hasCinsertion is true and hasInsertionAfter is false
            gene.addDetails(read, 'C insert + del')
            return True                                                 #In that case, must not have FS suppression
        elif residue531To534 == "SRTR"[0:len(residue531To534)]:         #Must have those residues due to insertion
            gene.addDetails(read, 'C insert')
            return True                                                 #Read has info in gene.additionalInfo (T) but is still susceptible (F); has insert (T)
        elif (shiftCount % 3) == 0:                                     #Because SRTR is not present, can't make prediction on resistance
            gene.addDetails(read, 'C insert')
            gene.addDetails(read, 'FS end')
            return False

def addRead(name, queryName, frameshiftInfoDict, additionalInformation = None):
    if name not in frameshiftInfoDict:
        frameshiftInfoDict.update({name:list()})
    if (additionalInformation == None) or (additionalInformation == ''):
        frameshiftInfoDict[name].append(queryName)
    else:
        frameshiftInfoDict[name].append((queryName, additionalInformation))