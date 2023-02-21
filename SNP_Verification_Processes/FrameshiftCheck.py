from SNP_Verification_Tools import dnaTranslate


# Description of input for all functions:
#   read:                               Reference to current AlignmentSegment 
#   gene:                               Reference to current Gene
def FrameshiftCheck(read, gene):
    # rRNA skips straight to extendedIndelCheck
    if gene.rRna():
        extendedIndelCheck(read, gene)
        return True

    # Gets counts of CIGAR operations
    cigarOpCount = read.get_cigar_stats()[0].tolist()

    # If (insertions - deletions) %3 != 0, susceptible unless if gene has frameshift info
    if (((cigarOpCount[1] - cigarOpCount[2]) % 3) != 0) and (gene.getGeneTag() == 'N') and (gene.getName() != "MEG_6142"): 
        gene.addDetails(read, 'FS till end')
        #Should not continue
        return False                                                            
    
    # Gene has two special cases where frameshift were found to be allowed in the literature
    elif gene.getName() == "MEG_6142": 
        if not(MEG_6142Check(read, gene)): return False

    # MEG_6094 can have a C insertion at res 531 if followed by frameshift suppression during translation
    elif gene.getGeneTag() == 'S':
        if not(MEG_6094Check(read, gene)): return False

    # Frameshift-type genes want 'FS till end'
    elif gene.getGeneTag() == 'F':
        if (((cigarOpCount[1] - cigarOpCount[2]) % 3) != 0):
            gene.addDetails(read, 'FS till end')
            extendedIndelCheck(read, gene)
            return True
    longFrameshiftCheck(read, gene)
    extendedIndelCheck(read, gene)
    return True


# Descripiton of top variables:
#   shift_count:                        Incremented by 1 for insertions, decremented by 1 for deletions
#   long_frameshift_count:              Number of frameshifts that are >= 12bp long
#   frameshift_length:                  Current length of frameshift, or 0 if not in a frameshift
#   pair_index:                         Current_index of aligned_pairs
#   aligned_pairs:                      List of tuples with aligned read and reference positions
#   reset:                              For MEG_6142, will for disregard of count information
#                                       if current read has a nonsense-causing frameshift at pos 26
#
# Counts amount of frameshifts that are that are >= 12bp long; adds this info to Gene object
def longFrameshiftCheck(read, gene):
    shift_count, frameshift_length, pair_index = 0, 0, -1
    long_frameshift_count = 0
    aligned_pairs = read.get_aligned_pairs()
    reset = gene.currentReadSpecial()
    for cigar_tuple in read.cigartuples:
        # Skip clipping
        if cigar_tuple[0] in [4,5]:
            continue

        pair_index += cigar_tuple[1]
        if (aligned_pairs[pair_index][1] != None) and (aligned_pairs[pair_index][1] >= 78) and reset:
            reset = False
            frameshift_length = 0
            long_frameshift_count = 0

        if (cigar_tuple[0] == 1):    shift_count += cigar_tuple[1]
        elif (cigar_tuple[0] == 2):  shift_count -= cigar_tuple[1]
        elif (cigar_tuple[0] == 0):
            if (shift_count % 3) != 0:  frameshift_length += cigar_tuple[1]
            elif (shift_count % 3) == 0:
                long_frameshift_count += 1 if frameshift_length >= 12 else 0
                frameshift_length = 0

    # MEG_6094 can have a C insertion may lead to a long frameshift which for this analysis should not be considered
    if (long_frameshift_count > 0) and gene.removeFromLongFrameshiftCheck(): long_frameshift_count -= 1

    # Adds to additional_info list in Gene object
    if (long_frameshift_count > 0): gene.addDetails(read, "12+bp frameshift: " + str(long_frameshift_count))
            
# Counts for indels greater than or equal to 12 base pairs long
def extendedIndelCheck(read, gene):
    indel = 0
    for cigar_tuple in read.cigartuples:
        if (cigar_tuple[0] in range(1,3)) and (cigar_tuple[1] >= 12):
            indel += 1
    if indel > 0:
        gene.addDetails(read, "12+bp indel: " + str(indel))

def MEG_6142Check(read, gene):
    query_sequence = read.query_alignment_sequence
    startIndex = read.cigartuples[0][1] if read.cigartuples[0][0] == 4 else 0
    endIndex = -1 * read.cigartuples[-1][1] if read.cigartuples[-1][0] == 4 else 0
    aligned_pairs = read.get_aligned_pairs()[startIndex:]               #Removes soft-clipping
    if endIndex != 0:
        aligned_pairs = aligned_pairs[:endIndex]
    shift_count = 0                                                     #If pos, more ins; if neg. more del
    queryIndex = -1
    stopCodon = ["TAA", "TGA", "TAG"]
    deletion76 = False
    fs482 = False
    for pair in aligned_pairs:
        if pair[0] == None:                                             #If deletion
            shift_count -= 1
            if pair[1] in range(74,76):
                deletion76 = True
        else:
            queryIndex += 1
            if pair[1] == None:                                         #If insertion
                shift_count += 1
        if (pair[1] == 75) and deletion76:                              #Change if next res isn't stop
            if (len(query_sequence) < queryIndex+4) or (query_sequence[queryIndex+1:queryIndex+4] not in stopCodon) or ((shift_count%3)!=2):
                deletion76 = False
        if (pair[1] == 78) and deletion76:
            shift_count = 0
        if (pair[1] == 1442) and ((shift_count%3)==0):
            if aligned_pairs[-1][1] == 1445:
                if len(query_sequence) < queryIndex+4:                   #If deletion
                    fs482 = True
                elif (query_sequence[queryIndex+1:queryIndex+4] not in stopCodon) and (len(query_sequence) > queryIndex+4):
                    fs482 = True
            
    if fs482:
        gene.updateCurrentReadNonstopInformation(True)
        gene.addDetails(read, "nonstop")
        if deletion76:
            gene.hasSpecialCase()
        return True
    else:
        gene.updateCurrentReadNonstopInformation(False)
        if (shift_count%3) == 0:
            if deletion76:
                gene.hasSpecialCase()
            return True
    gene.addDetails(read, 'FS till end')
    return False
        


def MEG_6094Check(read, gene):
    query_sequence = read.query_alignment_sequence
    startIndex = read.cigartuples[0][1] if read.cigartuples[0][0] == 4 else 0
    endIndex = -1 * read.cigartuples[-1][1] if read.cigartuples[-1][0] == 4 else 0
    aligned_pairs = read.get_aligned_pairs()[startIndex:]               #Removes soft-clipping
    if endIndex != 0:
        aligned_pairs = aligned_pairs[:endIndex]
    shift_count = 0                                                      #If pos, more ins; if neg. more del
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
    removeFromLongFrameshiftCheck = None
    for pair in aligned_pairs:
        refIndex += 1
        if pair[0] == None:                                             #If deletion
            if hasCinsertion:
                deletionCountAfterCinsertion += 1
                if (removeFromLongFrameshiftCheck == None) and ((refIndex - insertionPosition) != 0):
                    removeFromLongFrameshiftCheck = (refIndex - insertionPosition)>= 12
            shift_count -= 1
        elif pair[1] == None:                                           #If insertion
            refIndex -= 1
            if hasCinsertion:
                insertionCountAfterCinsertion += 1
                if (removeFromLongFrameshiftCheck == None) and ((refIndex - insertionPosition) != 0):
                    removeFromLongFrameshiftCheck = (refIndex - insertionPosition) >= 12
            if inCodon531 and ((shift_count%3)==0):                      #If in codon 531 and not in frameshift
                if not(hasCinsertion):
                    shift_count -= 1
                    hasCinsertion = True
                    insertionPosition = refIndex + 1
            shift_count += 1
        if (insertionCountAfterCinsertion > 0) and (deletionCountAfterCinsertion > 0):
            insertionCountAfterCinsertion -= 1                          #Never will have both over 1
            deletionCountAfterCinsertion -= 1
        if pair[1] == 1590:
            inCodon531 = True
        if pair[1] == 1593:
            inCodon531 = False
        if pair[0] == (lastBeforeFull + 3):
            if inCodon531 or (len(residue531To534) >= 1 and len(residue531To534) < 4):
                residue531To534 += dnaTranslate(query_sequence[lastBeforeFull+1:pair[0]+1], gene.getName())
            if inCodon531 or (len(residue531To536) >= 1 and len(residue531To536) < 6):
                residue531To536 += dnaTranslate(query_sequence[lastBeforeFull+1:pair[0]+1], gene.getName())
            lastBeforeFull += 3
    gene.hasLongFrameshift(removeFromLongFrameshiftCheck)
    twoInsertionsAfter = ((insertionCountAfterCinsertion%3) == 2)
    deletionAfter = ((deletionCountAfterCinsertion%3) == 1)
    if (shift_count % 3) == 2:
        if not(hasCinsertion):
            valid = False
    elif (shift_count % 3) == 1:
        valid = False
    if not(valid):
        gene.addDetails(read, 'FS till end')
        return False                                                    #Should not continue
    elif hasCinsertion:                                                 #shift_count%3 == 0 at residue 531
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
        elif (shift_count % 3) == 0:                                     #Because SRTR is not present, can't make prediction on resistance
            gene.addDetails(read, 'C insert + not SRTR')
            gene.addDetails(read, 'FS till end')
            return False
    return True