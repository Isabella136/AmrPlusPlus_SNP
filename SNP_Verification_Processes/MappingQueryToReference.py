from SNP_Verification_Tools import dnaTranslate
DEBUGGING_MODE = False

# Example: transforms 3M2I3M to MMMIIMMM
def extendCigar(cigar):
    extended_cigar = ""
    count = 0
    for c in cigar:
        if c.isdigit():
            count = count * 10 + int(c)
        else:
            while count > 0:
                extended_cigar = extended_cigar + c
                count -= 1
    return extended_cigar

# Description of input: 
#   read:                               Reference to current AlignmentSegment 
#   gene:                               Reference to current Gene
# Description of local variables:
#   cigar:                              Extended version of cigar string
#   aligned_pairs:                      List of tuples with aligned read and reference positions
#   query_seq:                          Query sequence
#   nt_alignment_map:                   List of tuples, where each tuple has the form 
#                                       (opcode, query pos, ref pos, ref nuc shift, and pre ref nuc shift)

def MapQueryToReference(read, gene, config):
    global DEBUGGING_MODE
    DEBUGGING_MODE = config.getboolean('SETTINGS', 'DEBUGGING_MODE')
    cigar = extendCigar(read.cigarstring)
    aligned_pair = read.get_aligned_pairs()
    query_seq = read.query_sequence

    # List of tuples, where each tuple has the form (opcode, query pos, ref pos, ref nuc shift, and pre ref nuc shift)
    nt_alignment_map = mapCigarToAlignment(cigar, aligned_pair, gene.rRna(), gene.mustSuppressFrameshift(), gene.currentReadSpecial())   

    # This can occur for MEG_6142 if query has a deletion at bp 76 and the query alignment to the reference sequence stops before bp 80
    if len(nt_alignment_map) == 0:
        return (False, False)
    
    #Makes nt_alignment_map have same format as aa_alignment_map (useful for rRNA)
    map_of_interest = transformNtAlignmentMap(nt_alignment_map)   

    # For MEG_6094, if gene.mustSuppressFrameshift() is True, calls suppressFS function
    if gene.mustSuppressFrameshift():
        soft_clipping_beginning = read.cigartuples[0][1] if read.cigartuples[0][0] == 4 else 0
        query_seq = suppressFS(1602, map_of_interest, query_seq, soft_clipping_beginning)

    # Trimms query sequence in order to not have broken codons during translation
    # for rRNA, the mapCigarToAlignment returns full map, so sequence won't be trimmed
    start = nt_alignment_map[0][1]
    i = 1
    while start == None:                                        
        start = nt_alignment_map[i][1]                          
        i += 1
    end = nt_alignment_map[len(nt_alignment_map)-1][1]
    i = 2
    while end == None:
        end = nt_alignment_map[len(nt_alignment_map)-i][1]
        i += 1
    if gene.mustSuppressFrameshift():
        end -= 1
    seq_of_interest = trimmed_query_sequence = query_seq[start:end+1]

    if not(gene.rRna()):
        seq_of_interest = dnaTranslate(trimmed_query_sequence, gene.getName())
        map_of_interest = aaAlignment(nt_alignment_map)         
    if DEBUGGING_MODE: print((seq_of_interest, map_of_interest))
    if DEBUGGING_MODE: print(gene.finalSequence()[list(map_of_interest.keys())[0]:list(map_of_interest.keys())[-1]+1])
    return (seq_of_interest, map_of_interest)


# Description of input: 
#   nucleotide_to_delete:               The reference position of the base pair that gets skiped over in translation
#   map_of_interest:                    Map reference postion to query postion(s)
#   query_sequence:                     Sequence of query that is aligned to reference
# 
# For MEG_6094, if gene.mustSuppressFrameshift() is True, the aligned_pair for reference position 1602
# was removed in map_of_interest; this function removes it from query_sequence as well

def suppressFS(nucleotide_to_delete, map_of_interest, query_sequence, soft_clipping_beginning):
    # Because aligned_pair for reference position 1602 was removed
    # Looks for query position that aligned to position 1601 and add one
    query_index_to_remove = map_of_interest[nucleotide_to_delete-1][0] + 1 + soft_clipping_beginning
    query_sequence = query_sequence[:query_index_to_remove] + query_sequence[query_index_to_remove+1:]
    return query_sequence

# Description of input:
#   cigar:                              Extended cigar string as outputted by extendCigar function
#   aligned_pair:                       List of tuples with aligned read and reference positions
#   rRNA:                               If currently aligned to rRna gene, is True; else is False
#   suppress:                           Boolean returned by gene.mustSuppressFrameshift()
#   special:                            Boolean returned by gene.currentReadSpecial()
#
# Description of local variables:
#   query_length:                       Length of trimmed query sequence
#   ref_length                          Length of reference sequence
#   index:                              Current aligned_pair index
#   shift:                              Incremented by 1 for insertions, decremented by 1 for deletions
#   prev_shift:                         Value of shift for previous iterated aligned_pair tuple
#   split_index:                        When query sequence alignment starts mid-codon, keeps track of ref_index % 3
#   suppressing_insertion:              Used if suppress is true; will be true from cytosine insertion to base pair 1602 (equivalent to ref_index 1601)
#   translatable:                       For alignments to proteins, will become true when reaching first full codon 
#   ref_index:                          Keeps track of current reference sequence
#   alignment_map:                      Reference to  list that is to be returned.  
#
# Trims map to not have broken codons if not rRNA and adapt map for unique edge cases in MEG_6094 and MEG_6142
# Returns a list of tuples, where each tuple has the form (opcode, query pos, ref pos, ref nuc shift, and pre ref nuc shift)

def mapCigarToAlignment(cigar, aligned_pair, rRNA, suppress, special):   
    query_length = ref_length = index = shift = prev_shift = 0
    split_index = -1
    suppressing_insertion = translatable = False 
    ref_index = None

    alignment_map = []

    # To facilitate aa alignment, can slightly tweak alignment in a manner such as that:
    # ...|SSS|SMI|MMM changed to ...|SSS|SMM|IMM
    # ...|SSS|SMI|III|MMM changed to ...|SSS|SMM|III|IMM
    # Requires keeping track of count of insertions that must be moved

    insertions_to_move = 0

    # Traverses through extended cigar string
    for op in cigar:
        # Not in aligned_pair, so skip
        if op in ["H", "P"]:
            continue

        # No indel
        elif op in ["M", "X", "="]:

            # Update prev_shift and ref_index
            prev_shift = shift
            ref_index = aligned_pair[index][1]

            #For MEG_6142, if has deletion at bp 76, ignores till bp 80 (reference index 79)
            if special :                                        
                if ref_index < 79:
                    index += 1
                    continue
                else:
                    special = False

            # If gene is protein and we previously haven't encountered its first fulll codon
            if not(translatable) and not(rRNA):

                # Define split_index if this is the first aligned_pair we encoutered; else increment by one
                split_index = aligned_pair[index][1] % 3 if split_index < 0 else split_index + 1

                # If we reached start of a full codon in query
                if (split_index % 3) == 0:
                    translatable = True
                    query_length += 1
                    if insertions_to_move == 0:
                        ref_length += 1
                        alignment_map.append(("M", aligned_pair[index][0], aligned_pair[index][1], shift, prev_shift))
                    else:
                        shift += 1
                        alignment_map.append(("I", aligned_pair[index][0], None, shift, prev_shift))
                        insertions_to_move -= 1

            # If gene is either rRNA or we've already encountered its first fulll codon
            else:

                # For MEG_6094, if query had suppressible cytosine insertion, mimic translation step by 
                # jumping over base pair 1602 (equivalent to ref_index 1601)
                if suppress and (aligned_pair[index][1] == 1601):
                    suppressing_insertion = False
                    index += 1
                    continue

                query_length += 1
                if insertions_to_move == 0:
                    ref_length += 1

                    # For MEG_6094, if query had suppressible cytosine insertion, maps query base pair to next reference base pair
                    alignment_map.append(("M", aligned_pair[index][0], ref_index+1 if suppressing_insertion else aligned_pair[index][1], shift, prev_shift))
                else:
                    shift += 1
                    alignment_map.append(("I", aligned_pair[index][0], None, shift, prev_shift))
                    insertions_to_move -= 1

        # Insertion
        elif op == "I":

            # Updates prev_shift; increment shift by one due to insertion, no change to ref_index because query is not aligned to reference
            prev_shift = shift
            shift += 1

            # For MEG_6142, if has deletion at bp 76, ignores till bp 80
            if special :
                if ref_index < 78:      # Because this is an insertion, having a ref_index 
                    index += 1          # of 78 singnifies that next ref_index will be 79 
                    continue            # (i.e. reference position will be 80)
                else:
                    special = False

            # If gene is protein and we previously haven't encountered its first fulll codon
            if not(translatable) and not(rRNA):

                # Define split_index if this is the first aligned_pair we encoutered; else increment by one
                if split_index < 0:
                    temp = aligned_pair[index][1]
                    i = 1                                   # Because this is an insertion, split_index can't be defined
                    while temp == None:                     # based on aligned_pair, so we must retrieve the first 
                        temp = aligned_pair[index+i][1]     # aligned pair without an insertion
                        i += 1
                    split_index = temp % 3
                else:
                    split_index += 1

                # If we reached start of a full codon in query
                if (split_index % 3) == 0:
                    translatable = True
                    query_length += 1
                    alignment_map.append((op, aligned_pair[index][0], aligned_pair[index][1], shift, prev_shift))
                else:
                    insertions_to_move += 1
                    shift -= 1

            # If gene is either rRNA or we've already encountered its first fulll codon
            else:

                # For MEG_6094, if query had suppressible cytosine insertion, mimic translation step by not counting this as an insertion
                if suppress and (ref_index + 1 > 1590) and not(suppressing_insertion):
                    suppressing_insertion = True
                    query_length += 1
                    if insertions_to_move == 0:
                        ref_length += 1
                        alignment_map.append(("M", aligned_pair[index][0], ref_index+1 + insertions_to_move, shift, prev_shift))
                        shift -= 1
                    else:
                        insertions_to_move -= 1
                        alignment_map.append((op, aligned_pair[index][0], aligned_pair[index][1], shift, prev_shift))
                else:
                    query_length += 1
                    alignment_map.append((op, aligned_pair[index][0], aligned_pair[index][1], shift, prev_shift))

        # Deletion
        elif op == "D":

            # Updates prev_shift and ref_index; decrement shift by one due to insertion
            prev_shift = shift
            shift -= 1
            ref_index = aligned_pair[index][1]

            #For MEG_6142, if has deletion at bp 76, ignores till bp 80 (reference index 79)
            if special :
                if ref_index < 79:
                    index += 1
                    continue
                else:
                    special = False

            # If gene is protein and we previously haven't encountered its first fulll codon
            if not(translatable) and not(rRNA):

                # Because this is a deletion, split_index doesn't get incremented
                # However, if split_index is currently 2, then the current reference base pair
                # is at the start of a new codon and should therefore be taken into account
                if (split_index == 2):
                    translatable = True
                    ref_length += 1
                    alignment_map.append((op, aligned_pair[index][0], aligned_pair[index][1], shift, prev_shift))

            else:
                ref_length += 1
                alignment_map.append((op, aligned_pair[index][0], aligned_pair[index][1], shift, prev_shift))

        # In aligned_pair, but irrelevant, so increment index then skip
        else:
            index += 1
            continue

        index += 1

    # For proteins, trims alignment so that only full codons are included
    if not(rRNA):
        while (query_length % 3) != 0:
            poped = alignment_map.pop()
            if poped[0] != "D" : query_length -= 1
    return alignment_map

# Description of input:
#   nt_alignment_map:           .       Reference to list that was return in the mapCigarToAlignment function
#
# Description of local variables:
#   new_nt_alignment_map:               Reference to map that is to be returned
#   nt_query_index:                     Index to portion of query nucleotide sequence that aligns to reference sequence
#   nt_ref_index:                       Index to reference sequence
# 
# Makes nt_alignment_map have same format as aa_alignment_map (useful for rRNA);
# format of map is {refNt, (correspondingQueryNt1, correspondingQueryNt2)}

def transformNtAlignmentMap(nt_alignment_map):
    new_nt_alignment_map = {}
    nt_query_index = 0 

    # Finds position of first reference base pair that aligns directly to query
    nt_ref_index = nt_alignment_map[0][2]
    for i in range (1, len(nt_alignment_map)+1):
        if i == len(nt_alignment_map):
            raise NameError("All Insertion in Alignment")
        nt_ref_index = nt_alignment_map[i][2]
        if nt_ref_index != None:
            break

    # Loops through nt_alignment_map
    for nt in nt_alignment_map:

        # Insertion
        if nt[0] == "I":

            # Maps inserted nucleotide to upcoming reference nucleotide
            if nt_ref_index not in new_nt_alignment_map:
                new_nt_alignment_map.update({nt_ref_index:tuple()})
            temp = list(new_nt_alignment_map[nt_ref_index])
            temp.append(nt_query_index)
            new_nt_alignment_map[nt_ref_index] = tuple(temp)

            # Also maps inserted nucelotide to previous reference nucleotide
            if nt_ref_index - 1 not in new_nt_alignment_map:
                new_nt_alignment_map.update({nt_ref_index-1:tuple()})
            temp = list(new_nt_alignment_map[nt_ref_index-1])
            temp.append(nt_query_index-1)
            new_nt_alignment_map[nt_ref_index-1] = tuple(temp)

            # Increment nt_query_index only
            nt_query_index += 1 

        # Deletion
        elif nt[0] == "D":

            # Maps None to current reference nucleotide
            if nt[2] not in new_nt_alignment_map:
                new_nt_alignment_map.update({nt[2]:tuple()})
            temp = list(new_nt_alignment_map[nt[2]])
            temp.append(None)
            new_nt_alignment_map[nt[2]] = tuple(temp)   

            # Increment nt_ref_index only
            nt_ref_index += 1

        # Not indel
        else:

            # Maps current query nucleotide to current reference nucleotide
            if nt[2] not in new_nt_alignment_map:
                new_nt_alignment_map.update({nt[2]:tuple()})
            temp = list(new_nt_alignment_map[nt[2]])
            temp.append(nt_query_index)
            new_nt_alignment_map[nt[2]] = tuple(temp)  

            # Increment nt_ref_index and nt_query_index
            nt_ref_index += 1
            nt_query_index += 1 

    return new_nt_alignment_map


# Description of input:
#   nt_alignment_map:           .       Reference to list that was return in the mapCigarToAlignment function
#
# Description of local variables:
#   aa_alignment_map:                   Reference to map that is to be returned
#   nt_query_index:                     Index to portion of query nucleotide sequence that aligns to reference sequence
#   aa_query_index:                     Index to reference sequence
#   insert_count:                       Number of consecutive insertions
#   delete_count:                       Number of consecutive deletions
#   nt_ref_index:                       Index to reference sequence
#   prev_aa_shift:                      When starting new reference codon, retrieve previous shift
#   first_alignment:                    References the first query codon aligned to current reference codon
#   last_alignment:                     References the last query codon aligned to current reference codon
#   inbetween:                          If in frameshift, inbetween should either be '-' or previous query codon
#   has_deletion:                       True if, and only if, current query codon has a deletion
# 
# For proteins, returns translated version of nt_alignment_map;
# format of map is {refCodon, (correspondingQueryCodon1, correspondingQueryCodon2)}

def aaAlignment(nt_alignment_map):
    aa_alignment_map = {}
    nt_query_index = aa_query_index = insert_count = delete_count = 0


    # Finds position of first reference base pair that aligns directly to query
    for i in range (0, len(nt_alignment_map)+1):
        if i == len(nt_alignment_map):
            raise NameError("All Insertion in Alignment")
        nt_ref_index = nt_alignment_map[i][2]
        if nt_ref_index != None:
            break

    prev_aa_shift = None    
    first_alignment = None
    last_alignment = 0
    has_deletion = False
    full_codon_deletion = False
    last_add_to_map_scenario = None
    first_ref_codon = True

    # Maps one query codon to one reference codon
    def addOneToOne(aa_query_index, aa_ref_index):
        aa = aa_alignment_map.get(aa_ref_index, False)
        if aa == False:
            aa_alignment_map.update({aa_ref_index:(aa_query_index,)})
        else:
            temp = list(aa)
            temp.append(aa_query_index)
            aa = tuple(temp)
            aa_alignment_map.update({aa_ref_index:aa})

    # Maps two query codons to one reference codon
    def addTwoToOne(aa_query_index1, aa_query_index2, aa_ref_index):
        addOneToOne(aa_query_index1, aa_ref_index)
        addOneToOne(aa_query_index2, aa_ref_index)

    # Maps one query codon to two reference codons
    def addOneToTwo(aa_query_index, aa_ref_index1, aa_ref_index2):
        addOneToOne(aa_query_index, aa_ref_index1)
        addOneToOne(aa_query_index, aa_ref_index2)

    def addToMapScenario(aa_ref_index, current_aa_shift):
        nonlocal last_add_to_map_scenario 

        # Scenario 1:   Full match (3M) or pseudo-full match (1D1I2M) 
        #               b/w query codon and ref codon
        # Add last query codon to current ref codon
        if ((nt_ref_index % 3 == 0)                                             # reached end of ref codon
                and (prev_aa_shift % 3 == 0)                                    # start of ref codon not in frameshift
                and (current_aa_shift % 3 == 0)                                 # end of ref codon is not in frameshift
                and (current_aa_shift - prev_aa_shift == 0)):                   # no insertions of length 3 or more (no new query codon introduced either)
            
            if DEBUGGING_MODE: print(str(nt_ref_index) + ' add to map scenario: 1')
            addOneToOne(last_alignment, aa_ref_index)
            last_add_to_map_scenario = 1

        # Scenario 2:   First full codon insertion (3I)
        # Add last query codon to previous ref codon
        elif (not has_deletion                                                  # no deletions allowed
                and (nt_ref_index % 3 == 0)                                     # between two reference codons
                and (prev_aa_shift % 3 == 0)                                    # start of ref codon not in frameshift
                and (current_aa_shift % 3 == 0)                                 # end of ref codon is not in frameshift
                and (first_alignment == last_alignment)                         # only insertions in query codon
                and (current_aa_shift - prev_aa_shift) == 3):                   # is first codon insertion
            
            if DEBUGGING_MODE: print(str(nt_ref_index) + ' add to map scenario: 2')
            addOneToOne(last_alignment, aa_ref_index-1)
            last_add_to_map_scenario = 2

        # Scenario 3:   Occurs for instances where ref codon can align 
        #               to a codon with insertion(s) followed by a codon 
        #               with deletion(s)
        # Change current ref codon to last two query codons 
        elif (has_deletion                                                      # deletions are allowed
                and not full_codon_deletion                                     # no full codon deletion
                and (nt_ref_index % 3 == 0)                                     # reached end of ref codon
                and (prev_aa_shift % 3 == 0)                                    # start of ref codon not in frameshift
                and (first_alignment <= (last_alignment - 1))                   # enough insertions present to start a new query codon
                and ((nt[3]-prev_aa_shift+delete_count)%3 == 0)):               # before first deletion, we left frameshift (insertion count is multiple of 3)
            
            if DEBUGGING_MODE: print(str(nt_ref_index) + ' add to map scenario: 3')
            if aa_ref_index in list(aa_alignment_map.keys()):
                aa_alignment_map.pop(aa_ref_index)
            addTwoToOne(last_alignment - 1, last_alignment, aa_ref_index)
            last_add_to_map_scenario = 3

        # Scenario 4:   We started or ended frameshift with insertions, 
        #               or we're right after the first full codon insertion
        # Add last two query codons to current ref codon
        elif ((nt_ref_index % 3 == 0)                                           # reached end of ref codon
                and (not has_deletion                                           # no deletions allowed if last scenario is 11
                    or last_add_to_map_scenario != 11)
                and ((prev_aa_shift % 3 == 0)                                   # start or end of ref codon not in frameshift
                    or (current_aa_shift % 3 == 0))
                and (first_alignment == (last_alignment - 1))):                 # we are or were in frameshift 

              
            if DEBUGGING_MODE: print(str(nt_ref_index) + ' add to map scenario: 4')
            addTwoToOne(first_alignment, last_alignment, aa_ref_index)
            last_add_to_map_scenario = 4

        # Scenario 5:   We're right after at least two full codon insertions
        # Change current ref codon to last two query codons 
        elif (not has_deletion                                                  # no deletions allowed
                and (nt_ref_index % 3 == 0)                                     # between two reference codons
                and (current_aa_shift % 3 == 0)                                 # end of ref codon not in frameshift
                and (first_alignment < (last_alignment - 1))                    # two or more full codon insertions
                and (aa_ref_index in list(aa_alignment_map.keys()))):           # must update map at ref codon
            
            if DEBUGGING_MODE: print(str(nt_ref_index) + ' add to map scenario: 5')
            aa_alignment_map.pop(aa_ref_index)
            addTwoToOne(last_alignment - 1, last_alignment, aa_ref_index)
            last_add_to_map_scenario = 5

        # Scenario 6:   We had a multiple of three (but not three) 
        #               insertions added inside ref codon
        # Add first and last query codons to current ref codon
        # If one query codon in the middle, also add it to current ref codon
        elif (not has_deletion                                                  # no deletions allowed
                and (nt_ref_index % 3 == 0)                                     # reached end of ref codon
                and (prev_aa_shift % 3 == 0)                                    # start of ref codon not in frameshift
                and (current_aa_shift % 3 == 0)                                 # end of ref codon is not in frameshift
                and (first_alignment < (last_alignment - 1))                    # multiple of three (but not three) insertions
                and (aa_ref_index not in list(aa_alignment_map.keys()))):       # ref codon not added to map previously (i.e insertions not b/w two ref codons)
            
            if DEBUGGING_MODE: print(str(nt_ref_index) + ' add to map scenario: 6')
            if (first_alignment+1) != (last_alignment-1): 
                addTwoToOne(first_alignment, last_alignment, aa_ref_index)
            else:
                addOneToOne(first_alignment, aa_ref_index)
                addTwoToOne(last_alignment - 1, last_alignment, aa_ref_index)
            last_add_to_map_scenario = 6

        # Scenario 7:   We had more than three insertions added 
        #               inside ref codon that ended frameshift
        # Add last two query codons to current ref codon
        elif (not has_deletion                                                  # no deletions allowed
                and (nt_ref_index % 3 == 0)                                     # reached end of ref codon
                and (prev_aa_shift % 3 != 0)                                    # start of ref codon in frameshift
                and (current_aa_shift % 3 == 0)                                 # end of ref codon is not in frameshift
                and (first_alignment < (last_alignment - 1))):                  # more than three insertions
            
            if DEBUGGING_MODE: print(str(nt_ref_index) + ' add to map scenario: 7')
            addTwoToOne(last_alignment - 1, last_alignment, aa_ref_index)
            last_add_to_map_scenario = 7


        # Scenario 8:   We had more than three insertions added 
        #               inside ref codon that started frameshift
        # Add first two query codons to current ref codon
        elif (not has_deletion                                                  # no deletions allowed
                and (nt_ref_index % 3 == 0)                                     # reached end of ref codon
                and (prev_aa_shift % 3 == 0)                                    # start of ref codon not in frameshift
                and (current_aa_shift % 3 != 0)                                 # end of ref codon currently in frameshift
                and (first_alignment < (last_alignment - 1))):                  # more than three insertions
            
            if DEBUGGING_MODE: print(str(nt_ref_index) + ' add to map scenario: 8')
            addTwoToOne(first_alignment, first_alignment+1, aa_ref_index)
            last_add_to_map_scenario = 8

        # Scenario 9: 
        # If (13 was not last called
        #     full_codon_deletion, 
        #     ref_index%3 == 0, 
        #     prev_aa_shift%3 == 0, 
        #     nt[3]%3 == 0, 
        #     prev_aa_shift - nt[3] == 3) : add '-' to aa_ref_index
        elif (full_codon_deletion
                and (nt_ref_index % 3 == 0)
                and (prev_aa_shift % 3 == 0)
                and (current_aa_shift % 3 == 0)
                and (last_add_to_map_scenario != 13)
                and ((prev_aa_shift - current_aa_shift) == 3)):
              
            if DEBUGGING_MODE: print(str(nt_ref_index) + ' add to map scenario: 9')
            addOneToOne('-', aa_ref_index)
            last_add_to_map_scenario = 9
        #10
        # If (first==last, 
        #     ref_index%3 == 0, 
        #     prev_aa_shift%3 == 0 or first_ref_codon,
        #     nt[3]%3 != 0,
        #     has_deletion or first_ref_codon) : add last to aa_ref_index
        elif ((nt_ref_index % 3 == 0)
                and (current_aa_shift % 3 != 0)
                and (has_deletion or first_ref_codon)
                and (first_alignment == last_alignment)
                and ((prev_aa_shift % 3 == 0) or first_ref_codon)):
            
            if DEBUGGING_MODE: print(str(nt_ref_index) + ' add to map scenario: 10')
            addOneToOne(last_alignment, aa_ref_index)
            last_add_to_map_scenario = 10
        #11
        # If (3 or 10 was last called,
        #     first==last-1 , 
        #     ref_index%3 == 0, 
        #     prev_aa_shift%3 != 0,
        #     nt[3]%3 != 0) : add first to aa_ref_index
        elif ((nt_ref_index % 3 == 0)
                and (prev_aa_shift % 3 != 0)
                and (current_aa_shift % 3 != 0)
                and (last_add_to_map_scenario in [3, 10])
                and (first_alignment == (last_alignment-1))):
                
            if DEBUGGING_MODE: print(str(nt_ref_index) + ' add to map scenario: 11')
            addOneToOne(first_alignment, aa_ref_index)
            last_add_to_map_scenario = 11
        #12
        # If (3 or 10 was last called,
        #     first==last, 
        #     ref_index%3 == 0, 
        #     nt[3]%3 == 0),
        #     has_deletion) : add '-' to aa_ref_index and aa_ref_index-1, add last to aa_ref_index
        elif (has_deletion
                and (nt_ref_index % 3 == 0)
                and (current_aa_shift % 3 == 0)
                and (first_alignment == last_alignment)
                and (last_add_to_map_scenario in [3, 10])):
            
            if DEBUGGING_MODE: print(str(nt_ref_index) + ' add to map scenario: 12')
            addOneToTwo('-', aa_ref_index-1, aa_ref_index)
            addOneToOne(last_alignment, aa_ref_index)
            last_add_to_map_scenario = 12
        #13
        # If (ref_index%3 != 0, 
        #     full_codon_deletion) : add '-' to aa_ref_index, aa_ref_index-1 unless if last_add_to_map_scenario == 8, 13 or None;
        #                            in that scenario, only add '-' to aa_ref_index
        elif (full_codon_deletion 
                and ((nt_ref_index % 3 != 0) or (current_aa_shift % 3 != 0))):

            if DEBUGGING_MODE: print(str(nt_ref_index) + ' add to map scenario: 13')
            if last_add_to_map_scenario not in [8, 11, 13, None]:
                addOneToTwo('-', aa_ref_index-1, aa_ref_index)
            else:
                addOneToOne('-', aa_ref_index)
            last_add_to_map_scenario = 13
        #14
        # If (first==last, 
        #     ref_index%3 == 0, 
        #     prev_aa_shift%3 != 0
        #     nt[3]%3 == 0,
        #     full_codon_deletion) : add last to aa_ref_index
        elif (full_codon_deletion
                and (nt_ref_index % 3 == 0)
                and (prev_aa_shift % 3 != 0)
                and (current_aa_shift % 3 == 0)
                and (first_alignment == last_alignment)):
            
            if DEBUGGING_MODE: print(str(nt_ref_index) + ' add to map scenario: 14')
            addOneToOne(last_alignment, aa_ref_index)
            last_add_to_map_scenario = 14
        #15
        # If (11 was last called,
        #     first==last, 
        #     ref_index%3 == 0, 
        #     nt[3]%3 == 0,
        #     has_deletion)  : add '-' to aa_ref_index-1, add last to aa_ref_index and aa_ref_index-1
        elif (has_deletion 
                and (nt_ref_index % 3 == 0) 
                and (current_aa_shift % 3 == 0) 
                and (last_add_to_map_scenario == 11) 
                and (first_alignment == last_alignment)):
            if DEBUGGING_MODE: print(str(nt_ref_index) + ' add to map scenario: 15')
            addOneToTwo(last_alignment, aa_ref_index-1, aa_ref_index)
            addOneToOne('-', aa_ref_index-1)
            last_add_to_map_scenario = 15
        #16
        # If (first==last, 
        #     ref_index%3 == 0, 
        #     nt[3]%3 == 0,
        #     has_deletion
        #     not full_codon_deletion)  : add last to aa_ref_index and aa_ref_index-1
        elif (has_deletion
                and not full_codon_deletion
                and (nt_ref_index % 3 == 0) 
                and (current_aa_shift % 3 == 0) 
                and (first_alignment == last_alignment)):
            if DEBUGGING_MODE: print(str(nt_ref_index) + ' add to map scenario: 16')
            addOneToTwo(last_alignment, aa_ref_index-1, aa_ref_index)
            last_add_to_map_scenario = 16

        else: 
            if DEBUGGING_MODE: print(str(nt_ref_index) + ' add to map scenario: else')
            last_add_to_map_scenario = None

    for count, nt in enumerate(nt_alignment_map):
        # If previous alignment pair was last query base pair in codon
        if ((nt_query_index % 3) == 0) and (delete_count == 0): 
            has_deletion = False

        # In the improbale case that insertion and deletion are together, combine them into one:
        combine = False
        if (nt[0] == "I") :
            if (count != 1) and (nt_alignment_map[count-1][0] == "D"): continue
            if (count != len(nt_alignment_map)-1) and (nt_alignment_map[count+1][0] == "D"): combine = True
        if (nt[0] == "D") :
            if (count != 1) and (nt_alignment_map[count-1][0] == "I"): continue
            if (count != len(nt_alignment_map)-1) and (nt_alignment_map[count+1][0] == "I"): combine = True

        # If we reached start of new reference codon
        if (nt_ref_index % 3) == 0:

            if prev_aa_shift == None:
                # Redefines prev_aa_shift because it is a new reference codon
                prev_aa_shift = nt[4]

            # Insertion
            if (nt[0] == "I") and not(combine):
                # Increment insert_count and nt_query_index by one, update delete_count to zero
                insert_count += 1 
                nt_query_index += 1
                delete_count = 0

                # If we didn't have an insertion beforehand
                if first_alignment == None: first_alignment = aa_query_index

                # If shift in alignment is now a multiple of three
                if (nt[3] % 3) == 0:  

                    # If there has been three consecutive insertions
                    if insert_count % 3 == 0:
                        addToMapScenario(int(nt_ref_index//3), nt[3])

                    # We are leaving a frameshift and we're heading into a new codon
                    elif (count == len(nt_alignment_map)-1) or nt_alignment_map[count+1][0] != "I":
                        if DEBUGGING_MODE: print(str(nt_ref_index) + ' add to map scenario: 4 special')
                        addOneToOne(last_alignment, int(nt_ref_index//3))
                        first_alignment = None
                        prev_aa_shift = None
                
            # Deletion
            elif (nt[0] == "D") and not(combine):
                # Increment delete_count and nt_ref_index by one, update insert_count to zero, set has_deletion to true
                nt_ref_index += 1
                delete_count += 1
                insert_count = 0
                has_deletion = True

                # For scenario such as the codon MDDDMM, maps a deletion to previous and current reference codons
                # For scenarios such as DMMDDDDDM, maps deletion to current reference codon since deletions were mapped to last
                #       two codons in line 849
                if (count == len(nt_alignment_map)-1) or nt_alignment_map[count+1][0] != "D":
                    if (delete_count >= 3) and (prev_aa_shift % 3 != 0) and (nt[3] %3 == 0):
                        full_codon_deletion = True
                        addToMapScenario(int(nt_ref_index//3), nt[3])

            # Not indel
            else:

                # Increment nt_query_index and nt_ref_index by one, update insert_count and delete_count to zero
                nt_query_index += 1
                nt_ref_index += 1
                insert_count = 0
                delete_count = 0

                # If we didn't have an insertion beforehand
                if first_alignment == None: first_alignment = aa_query_index

                # This is the last query codon and it has two insertions
                if count == (len(nt_alignment_map) - 1) and (prev_aa_shift%3)==0:
                    if DEBUGGING_MODE: print(str(nt_ref_index) + ' add to map scenario: end special')
                    addOneToOne(last_alignment, int(nt_ref_index//3))
                # This is the last query codon and the previous codon started a frameshift through a deletion 
                elif count == (len(nt_alignment_map) - 1) and last_add_to_map_scenario==8:
                    if DEBUGGING_MODE: print(str(nt_ref_index) + ' add to map scenario: end special')
                    addOneToOne(last_alignment, int(nt_ref_index//3))


        # If we reached the middle of a reference codon        
        elif (nt_ref_index % 3) == 1:

            # Insertion
            if (nt[0] == "I") and not(combine):

                # Increment nt_query_index and insert_count by one, update delete_count to zero
                nt_query_index += 1
                insert_count += 1 
                delete_count = 0

                # This is the last query codon and it has two insertions
                if count == (len(nt_alignment_map) - 1):
                    if DEBUGGING_MODE: print(str(nt_ref_index) + ' add to map scenario: end special')
                    addOneToOne(last_alignment, int(nt_ref_index//3))

            # Deletion
            elif (nt[0] == "D") and not(combine):

                # Increment delete_count and nt_ref_index by one, update insert_count to zero, set has_deletion to true
                nt_ref_index += 1
                delete_count += 1
                insert_count = 0
                has_deletion = True

                # For scenario such as the codon MMDDDM, maps a deletion to previous and current reference codons
                # For scenarios such as DMMDDDDDM, maps deletion to current reference codon since deletions were mapped to last
                #       two codons in line 849
                if (count == len(nt_alignment_map)-1) or nt_alignment_map[count+1][0] != "D":
                    if (delete_count >= 3) and (prev_aa_shift % 3 != 0) and (nt[3] %3 == 0):
                        full_codon_deletion = True
                        addToMapScenario(int(nt_ref_index//3), nt[3])

            # Not indel
            else:

                # Increment nt_query_index and nt_ref_index by one, update insert_count and delete_count to zero
                nt_query_index += 1
                nt_ref_index += 1
                insert_count = 0
                delete_count = 0

                # Would happen for scenarios where first alignment pair in nt_alignment_map is an insertion 
                # and second alignemnt pair is a match to the second base pair in the reference codon
                if prev_aa_shift == None: prev_aa_shift = nt[4]

                # If the first base in the reference codon was deleted
                if first_alignment == None: first_alignment = aa_query_index

                # This is the last query codon and it has one insertion
                if count == (len(nt_alignment_map) - 1) and (prev_aa_shift%3)==0:
                    if DEBUGGING_MODE: print(str(nt_ref_index) + ' add to map scenario: end special')

                    addOneToOne(last_alignment, int(nt_ref_index//3))
                # This is the last query codon and the previous codon started a frameshift through two deletions 
                elif count == (len(nt_alignment_map) - 1) and last_add_to_map_scenario==8:
                    if DEBUGGING_MODE: print(str(nt_ref_index) + ' add to map scenario: end special')

                    addOneToOne(last_alignment, int(nt_ref_index//3))

        # If we reached the end of a reference codon   
        else:

            # Insertion
            if (nt[0] == "I") and not(combine):

                # Increment nt_query_index and insert_count by one, update delete_count to zero
                nt_query_index += 1
                insert_count += 1 
                delete_count = 0

                # This is the last query codon and it has one insertion
                if count == (len(nt_alignment_map) - 1):
                    if DEBUGGING_MODE: print(str(nt_ref_index) + ' add to map scenario: end special')

                    addOneToOne(last_alignment, int(nt_ref_index//3))

            # Deletion
            elif (nt[0] == "D") and not(combine):

                # Increment delete_count and nt_ref_index by one, update insert_count to zero, set has_deletion to true
                nt_ref_index += 1
                delete_count += 1
                insert_count = 0
                has_deletion = True

                # For reference: in cases where we have an extended cigar like SMDMDD where 
                # the first M is considered as the second base of a codon,
                # we consider the current reference codon to be fully deleted
                if first_alignment == None: 
                    full_codon_deletion = True

                elif nt_query_index % 3 == 0:
                    last_alignment -= 1
                
                # Run addToMapScenario if one of these applies:
                if ((nt[3] % 3 == 0)                                        # we have left frameshift                                
                        or first_ref_codon                                  # first reference codon aligned to query
                        or full_codon_deletion                              # entire reference codon has been deleted
                        or (prev_aa_shift % 3 == 0)                         # we weren't in a frameshift at the start
                        or (last_add_to_map_scenario in [3, 10, 11])):      # previous reference codon went through scenario 3, 10, or 11

                    addToMapScenario(int(nt_ref_index//3)-1, nt[3])

                prev_aa_shift = None
                first_alignment = None
                first_ref_codon = False
                full_codon_deletion = False

            # Not indel
            else: #nt[0] == "M"

                # Increment nt_query_index and nt_ref_index by one, update insert_count and delete_count to zero
                nt_query_index += 1
                nt_ref_index += 1
                insert_count = 0
                delete_count = 0

                # Would happen for scenarios where first two alignment pairs in nt_alignment_map are insertions 
                # and third alignemnt pair is a match to the last base pair in the reference codon
                if prev_aa_shift == None: prev_aa_shift = nt[4]

                # If the first two bases in the reference codon was deleted
                if first_alignment == None: first_alignment = aa_query_index

                #addToMapScenario is always called at the end of a reference codon if last nt is a match
                addToMapScenario(int(nt_ref_index//3)-1, nt[3])

                prev_aa_shift = None
                first_alignment = None
                full_codon_deletion = False
                first_ref_codon = False

        aa_query_index = int(nt_query_index / 3)
        last_alignment = aa_query_index

    return aa_alignment_map    