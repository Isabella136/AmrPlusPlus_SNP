from SNP_Verification_Tools import dnaTranslate
DEBUGGING_MODE = False

# Description of input for all functions:
#   read:                               Reference to current AlignmentSegment 
#   gene:                               Reference to current Gene
def FrameshiftCheck(read, gene, config):
    global DEBUGGING_MODE
    DEBUGGING_MODE = config.getboolean('SETTINGS', 'DEBUGGING_MODE')
    # rRNA skips straight to extendedIndelCheck
    if gene.rRna():
        extendedIndelCheck(read, gene)
        return True

    # Gets counts of CIGAR operations
    cigarOpCount = read.get_cigar_stats()[0].tolist()

    # If (insertions - deletions) %3 != 0, susceptible unless if gene has frameshift info
    if (((cigarOpCount[1] - cigarOpCount[2]) % 3) != 0) and (gene.getGeneTag() in  ['N', 'H', 'I']) and (gene.getName() != "MEG_6142"): 
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


# Descripiton of local variables:
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
        if (cigar_tuple[0] in range(1,3)) and (cigar_tuple[1] >= 12): indel += 1
    if indel > 0:
        gene.addDetails(read, "12+bp indel: " + str(indel))

# Description of local variables:
#   shift_count:                        Incremented by 1 for insertions, decremented by 1 for deletions
#   query_index:                        Current index of query sequence
#   deletion76:                         True if, and only if, there is a deletion between reference index 74 and 76
#                                       that leads to a stop at codon 26
#   fs482:                              True if, and only if, there is a nonstop-causing frameshift in the last codon of the gene
#   aligned_pairs:                      List of tuples with aligned read and reference positions
#   query_sequence:                     Sequence of read portion that aligned to MEGARes
#
# Specifically checks alignements to MEG_6142 which can allow for two specific base indel without removing gene function. 
# This includes a nonstop mutation that confers resistance

def MEG_6142Check(read, gene):
    shift_count, query_index = 0 , -1
    deletion76, fs482 = False, False

    # Remove soft-clipping
    start_index = read.cigartuples[0][1] if read.cigartuples[0][0] == 4 else 0
    end_index = read.cigartuples[-1][1] if read.cigartuples[-1][0] == 4 else 0
    aligned_pairs = read.get_aligned_pairs()[start_index:len(read.get_aligned_pairs()) - end_index]
    
    query_sequence = read.query_alignment_sequence
    stop_codon = ["TAA", "TGA", "TAG"]
    
    for pair in aligned_pairs:
        query_index += 1

        # Looks for deletion
        if pair[0] == None:
            shift_count -= 1
            query_index -= 1
            if pair[1] in range(74,76):
                deletion76 = True

        # Looks for insertion
        elif pair[1] == None:
            shift_count += 1

        # Checks to see if codon 26 is now a stop
        if (pair[1] == 75) and deletion76:
            if (len(query_sequence) < query_index+4) or (query_sequence[query_index+1:query_index+4] not in stop_codon) or ((shift_count%3)!=2):
                deletion76 = False

        # Start over the frameshift check
        if (pair[1] == 78) and deletion76:
            shift_count = 0

        # Looks for nonstop-causing frameshift
        if (pair[1] == 1442) and ((shift_count%3)==0):
            if aligned_pairs[-1][1] == 1445:
                if len(query_sequence) < query_index+4:
                    fs482 = True
                elif (query_sequence[query_index+1:query_index+4] not in stop_codon) and (len(query_sequence) > query_index+4):
                    fs482 = True

    # Update output info accordingly        
    if fs482:
        gene.updateCurrentReadNonstopInformation(True)
        if DEBUGGING_MODE:
            print('Nonstop')
        gene.addDetails(read, "nonstop")
        if deletion76:
            if DEBUGGING_MODE:
                print('special deletion at 76')
            gene.hasSpecialCase()
        return True
    else:
        gene.updateCurrentReadNonstopInformation(False)
        if (shift_count%3) == 0:
            if deletion76:
                if DEBUGGING_MODE:
                    print('special deletion at 76')
                gene.hasSpecialCase()
            return True
    if DEBUGGING_MODE:
        print('FS till end')
    gene.addDetails(read, 'FS till end')

    # Return False if frameshift is unexplained
    return False
        

# Description of local variables:
#   shift_count:                        Incremented by 1 for insertions, decremented by 1 for deletions
#   insertion_count_after_C_insertion:  Number of insertions after the insertion at codon 531 that weren't negated by deletions; 
#                                       if insertion_count_after_C_insertion % 3 == 2, the insertion at codon 531 shouldn't be 
#                                       suppressed in this analysis
#   deletion_count_after_C_insertion:   Number of deltions after the insertion at codon 531 that weren't negated by insertions;
#                                       if deletion_count_after_C_insertion % 3 == 1, the insertion at codon 531 shouldn't be 
#                                       suppressed in this analysis
#   has_C_insertion:                    True if, and only if, there is a cytosine insertion in codon 531 that causes a frameshift
#   residue_531_to_534/536:             Keeps track of amino acids in query residues 531 to 534 or 536
#   aligned_pairs:                      List of tuples with aligned read and reference positions
#   query_sequence:                     Sequence of read portion that aligned to MEGARes
#   last_before_full:                   Keeps track of last query nucleotide index before next codon
#   ref_index:                          Current index of reference sequence
#   in_codon531:                        True only when analyzing codon 531
#   valid:                              False if, and only if, there is still a frameshift by the end of the query
#   insertion_position:                 If has_C_insertion is True, contains the position of the insertion; else is None
#   remove_from_long_frameshift_check:  If cytosine insertion in codon 531 causes a frameshift that isn't suppressible,
#                                       Will be True if frameshift is long and False if frameshift is short; otherwise is None
# 
# Specifically checks alignments to MEG_6094 which can have a cytosine insertion in codon 531.
# If present, this insertion can be suppressed during protein translation and confer resistance.

def MEG_6094Check(read, gene):
    shift_count, insertion_count_after_C_insertion, deletion_count_after_C_insertion = 0, 0, 0
    has_C_insertion = False
    residue_531_to_536, residue_531_to_534 = "", ""

    # Remove soft-clipping
    start_index = read.cigartuples[0][1] if read.cigartuples[0][0] == 4 else 0
    end_index = read.cigartuples[-1][1] if read.cigartuples[-1][0] == 4 else 0
    aligned_pairs = read.get_aligned_pairs()[start_index:len(read.get_aligned_pairs()) - end_index]

    query_sequence = read.query_sequence
    last_before_full = (3 - (aligned_pairs[0][1]) % 3) % 3 - 1 + start_index
    ref_index = aligned_pairs[0][1]

    in_codon531 = True if aligned_pairs[0][1] in range(1591,1594) else False

    valid = True
    insertion_position = None
    remove_from_long_frameshift_check = None

    for pair in aligned_pairs:
        ref_index += 1

        # Looks for deletion
        if pair[0] == None:
            if has_C_insertion:
                deletion_count_after_C_insertion += 1
                if (remove_from_long_frameshift_check == None) and ((ref_index - insertion_position) != 0):
                    remove_from_long_frameshift_check = (ref_index - insertion_position)>= 12
            shift_count -= 1

        # Looks for insertion
        elif pair[1] == None:
            ref_index -= 1
            if has_C_insertion:
                insertion_count_after_C_insertion += 1
                if (remove_from_long_frameshift_check == None) and ((ref_index - insertion_position) != 0):
                    remove_from_long_frameshift_check = (ref_index - insertion_position) >= 12

            # If in codon 531 and not currently in frameshift
            if in_codon531 and ((shift_count%3)==0):
                if not(has_C_insertion):
                    shift_count -= 1
                    has_C_insertion = True
                    insertion_position = ref_index + 1
            shift_count += 1

        # Never will have both over 0 at the same time
        if (insertion_count_after_C_insertion > 0) and (deletion_count_after_C_insertion > 0):
            insertion_count_after_C_insertion -= 1                         
            deletion_count_after_C_insertion -= 1

        # Sets in_codon531 as True for next iteration
        if pair[1] == 1590:
            in_codon531 = True

        # Sets in_codon531 as False for next iteration
        if pair[1] == 1593:
            in_codon531 = False

        # Updates last_before_full; keep track of amino acids for variables residue_531_to_534/536, if needed
        if pair[0] == (last_before_full + 3):
            if in_codon531 or (len(residue_531_to_534) >= 1 and len(residue_531_to_534) < 4):
                residue_531_to_534 += dnaTranslate(query_sequence[last_before_full+1:pair[0]+1], gene.getName())
            if in_codon531 or (len(residue_531_to_536) >= 1 and len(residue_531_to_536) < 6):
                residue_531_to_536 += dnaTranslate(query_sequence[last_before_full+1:pair[0]+1], gene.getName())
            last_before_full += 3

    gene.updateLongFrameshift(remove_from_long_frameshift_check)
    two_insertions_after = ((insertion_count_after_C_insertion%3) == 2)
    deletion_after = ((deletion_count_after_C_insertion%3) == 1)

    # Checks for validity
    if (shift_count % 3) == 2 and not(has_C_insertion): valid = False
    elif (shift_count % 3) == 1: valid = False

    # Alignment should not continue through SNP_Verification
    if not(valid):
        if DEBUGGING_MODE:
            print('Not valid')
        gene.addDetails(read, 'FS till end')
        return False

    # Unless if deletion_after or two_insertions_after are True, shift_count % 3 == 0
    elif has_C_insertion:
        if DEBUGGING_MODE:
            print('C insertion')      
        # Must have thoe residues SRTR due to insertion                                             
        if (residue_531_to_534 == "SRTR"):    
            # In that case, must not have FS suppression in analysis     
            if (deletion_after or two_insertions_after):
                if DEBUGGING_MODE:
                    print('C insert followed by del/ins')
                gene.addDetails(read, 'C insert followed by del/ins')                                                  
            elif (residue_531_to_536 == "SRTRPR"):  
                if DEBUGGING_MODE:
                    print('Suppressible C insert')                                                     
                gene.addDetails(read, 'Suppressible C insert')
            # Suppression won't happen during RNA translation to protein in the cell if there is no PR in residues 535 and 536
            else:
                if DEBUGGING_MODE:
                    print('C insert + not SRTRPR')
                gene.addDetails(read, 'C insert + not SRTRPR')          
                gene.addDetails(read, 'FS till end')
                return False       
            return True 
        # Because SRTR is not present, can't make prediction on resistance
        elif (shift_count % 3) == 0:     
            if DEBUGGING_MODE:
                print('C insert + not SRTR')                                    
            gene.addDetails(read, 'C insert + not SRTR')
            gene.addDetails(read, 'FS till end')
            return False
    if DEBUGGING_MODE: print('not suppressible but valid')
    return True