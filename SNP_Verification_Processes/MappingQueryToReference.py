from SNP_Verification_Tools import dnaTranslate


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

def MapQueryToReference(read, gene):
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
         query_seq = suppressFS(1602, mapOfInterest, query_seq)

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

    # ?????????? what did I meant by extend past
    if not(gene.rRna()):
        seq_of_interest = dnaTranslate(trimmed_query_sequence, gene.getName())
        mapOfInterest = aaAlignment(nt_alignment_map)           #for 'F' type of genes, seq_of_interest may extend past map_of_interest
                                                                #due to frameshift that extend past the end of query sequence
    return (seq_of_interest, map_of_interest)


# Description of input: 
#   nucleotide_to_delete:               The reference position of the base pair that gets skiped over in translation
#   map_of_interest:                    Map reference postion to query postion(s)
#   query_sequence:                     Sequence of query that is aligned to reference
# 
# For MEG_6094, if gene.mustSuppressFrameshift() is True, the aligned_pair for reference position 1602
# was removed in map_of_interest; this function removes it from query_sequence as well

def suppressFS(nucleotide_to_delete, map_of_interest, query_sequence):
    # Because aligned_pair for reference position 1602 was removed
    # Looks for query position that aligned to position 1601 and add one
    query_index_to_remove = map_of_interest[nucleotide_to_delete-1][0]+1
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
                    ref_length += 1
                    alignment_map.append(("M", aligned_pair[index][0], aligned_pair[index][1], shift, prev_shift))

            # If gene is either rRNA or we've already encountered its first fulll codon
            else:

                # For MEG_6094, if query had suppressible cytosine insertion, mimic translation step by 
                # jumping over base pair 1602 (equivalent to ref_index 1601)
                if suppress and (aligned_pair[index][1] == 1601):
                    suppressing_insertion = False
                    index += 1
                    continue

                query_length += 1
                ref_length += 1

                # For MEG_6094, if query had suppressible cytosine insertion, maps query base pair to next reference base pair
                alignment_map.append(("M", aligned_pair[index][0], ref_index+1 if suppressing_insertion else aligned_pair[index][1], shift, prev_shift))

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

            # If gene is either rRNA or we've already encountered its first fulll codon
            else:

                # For MEG_6094, if query had suppressible cytosine insertion, mimic translation step by not counting this as an insertion
                if suppress and (ref_index + 1 > 1590) and not(suppressing_insertion):
                    suppressing_insertion = True
                    query_length += 1
                    ref_length += 1
                    alignment_map.append(("M", aligned_pair[index][0], ref_index+1, shift, prev_shift))
                    index += 1
                    shift -= 1

                else:
                    query_length += 1
                    alignment_map.append((op, aligned_pair[index][0], aligned_pair[index][1], shift, prev_shift))
                    index += 1

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
    nt_ref_index = nt_alignment_map[0][2]
    for i in range (1, len(nt_alignment_map)+1):
        if i == len(nt_alignment_map):
            raise NameError("All Insertion in Alignment")
        nt_ref_index = nt_alignment_map[i][2]
        if nt_ref_index != None:
            break

    prev_aa_shift = None    
    first_alignment = None
    last_alignment = 0
    inbetween = None
    has_deletion = False


    # Maps query codon in inbetween to previous reference codon then calls addToMap function
    def addToMapDeletion(index, toAdd):
        aa = aa_alignment_map.get(index-1, False)
        if aa == False:
            aa_alignment_map.update({index-1:(inbetween, )})
        else:
            temp = list(aa)
            if inbetween not in temp:
                temp.append(inbetween)
            aa = tuple(temp)
            aa_alignment_map.update({index-1:aa})
        addToMap(index, toAdd)

    # Adds both current query codon and inbetween query codon to same reference codon
    def addToMapInbetween(index, toAdd):
        aa = aa_alignment_map.get(index, False)
        if aa == False:
            if (inbetween != None):
                aa_alignment_map.update({index:(inbetween, toAdd, )})
            else:
                aa_alignment_map.update({index:(toAdd,)})
        else:
            temp = list(aa)
            if (inbetween != None):
                temp.append(inbetween)
            temp.append(toAdd)
            aa = tuple(temp)
            aa_alignment_map.update({index:aa})

    # Maps query codon to reference codon
    def addToMap(index, toAdd):
        aa = aa_alignment_map.get(index, False)
        if aa == False:
            aa_alignment_map.update({index:(toAdd,)})
        else:
            temp = list(aa)
            temp.append(toAdd)
            aa = tuple(temp)
            aa_alignment_map.update({index:aa})

    for count, nt in enumerate(nt_alignment_map):
        # If previous alignment pair was last query base pair in codon
        if ((nt_query_index % 3) == 0) and (delete_count == 0): has_deletion = False

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
                    if insert_count == 3:
                        # That didn't happend right after another inserted codon
                        if inbetween == None:
                            # Inserted codon mapped to previous reference codon
                            addToMap(int(nt_ref_index/3)-1, last_alignment)

                        # Defines inbetween as current query codon so that it gets also mapped to upcoming reference codon
                        inbetween = last_alignment

                    # If shift % 3 == 0 and insert_count is not 3, then frameshift is present and if of form #1I...1M2I or 2I...2M1I 
                    # Inbetween is previous query codon, and both the previous and current query codons are mapped to the previous reference codon
                    else:
                        addToMapInbetween(int(nt_ref_index/3)-1, last_alignment)
                        inbetween = None

                    # Resets insert_count
                    insert_count = 0
                
            # Deletion
            elif (nt[0] == "D") and not(combine):
                # Increment delete_count and nt_ref_index by one, update insert_count to zero, set has_deletion to true
                nt_ref_index += 1
                delete_count += 1
                insert_count = 0
                has_deletion = True

                # For scenario such as the codon MDDDMM, maps a deletion to previous reference codon
                # Since the current query codon has already been mapped in line 
                if delete_count == 3:
                    addToMap(int(nt_ref_index/3)-1, '-')
                    inbetween = '-'

            # Not indel
            else:

                # Increment nt_query_index and nt_ref_index by one, update insert_count and delete_count to zero
                nt_query_index += 1
                nt_ref_index += 1
                insert_count = 0
                delete_count = 0

                # We started frameshift with either two insertions or one deletion, and continued with match eversince
                if (nt_query_index % 3) == 0:

                    # Case: we were not in a frameshift, but we started the query codon with two insertions
                    if (prev_aa_shift % 3) == 0: inbetween = None

                    # Case: we are in a frameshift that started in the previous reference codon
                    elif inbetween == None:

                        # Case: we are in the second codon in IIM|MMM and must map to first codon in ABC|A..
                        if not(has_deletion):
                            addToMap(int(nt_ref_index/3)-1, last_alignment)

                        # Case: we are in the first codon in MDMM|MM... and must map to last codon in ABC|ABC
                        if (nt[3] < 0) or (combine and nt_alignment_map[count][3] < 0):
                            addToMap(int(nt_ref_index/3), last_alignment)

                        inbetween = last_alignment 

                    # Case: we are in a frameshift that started before the previous reference codon
                    else: inbetween = last_alignment

                # If we didn't have an insertion beforehand
                if first_alignment == None: first_alignment = aa_query_index


        # If we reached the middle of a reference codon        
        elif (nt_ref_index % 3) == 1:

            # Insertion
            if (nt[0] == "I") and not(combine):

                # Increment nt_query_index and insert_count by one, update delete_count to zero
                nt_query_index += 1
                insert_count += 1 
                delete_count = 0

                # Case: we are at our second insertion and reached the end of the query codon
                if (nt_query_index % 3) == 0:
                    
                    # Case: we didn't start the query codon in a frameshift, but we encountered two insertions since then
                    if (prev_aa_shift % 3) == 0:
                        inbetween = None

                    # Case: we are in the second codon in IMM|MMI and must map to the first codon in ABC|A...
                    elif inbetween == None:
                        addToMap(int(nt_ref_index/3)-1, last_alignment)
                        inbetween = last_alignment

                    # Case: the first insertion was before the prevous query codon
                    else:
                        inbetween = last_alignment

            # Deletion
            elif (nt[0] == "D") and not(combine):

                # Increment delete_count and nt_ref_index by one, update insert_count to zero, set has_deletion to true
                nt_ref_index += 1
                delete_count += 1
                insert_count = 0
                has_deletion = True

                # For scenario such as the codon MMDDDM, maps a deletion to previous reference codon
                # Since the current query codon has already been mapped in line 
                if delete_count == 3:
                    addToMap(int(nt_ref_index/3)-1, '-')
                    inbetween = '-'

            # Not indel
            else:

                # Increment nt_query_index and nt_ref_index by one, update insert_count and delete_count to zero
                nt_query_index += 1
                nt_ref_index += 1
                insert_count = 0
                delete_count = 0

                # If the first base in the reference codon was deleted
                if first_alignment == None: first_alignment = aa_query_index

                # Would happen for scenarios where first alignment pair in nt_alignment_map is an insertion 
                # and second alignemnt pair is a match to the second base pair in the reference codon
                if prev_aa_shift == None: prev_aa_shift = nt[4]

                # We started frameshift with either one insertion or two deletions, and continued with match ever since
                if (nt_query_index % 3) == 0:

                    # Case: we weren't in a frameshift, but we encountered an insertion during this reference codon
                    if (prev_aa_shift % 3) == 0: inbetween = None

                    # Case: we are in a frameshift that started in the previous reference codon
                    elif inbetween == None:

                        # Case: we are in the second codon in IMM|MMM and must map to the first codon in ABC|AB...
                        if not(has_deletion):
                            addToMap(int(nt_ref_index/3)-1, last_alignment)

                        # Case: we are in the first codon in MDDMM|M... and must map to the last codon in ABC|ABC
                        if (nt[3] < 0) or (combine and nt_alignment_map[count][3] < 0):
                            addToMap(int(nt_ref_index/3), last_alignment)

                        inbetween = last_alignment

                    # Case: we are in a frameshift that started before the previous reference codon
                    else:
                        inbetween = last_alignment

        # If we reached the end of a reference codon   
        else:

            # Insertion
            if (nt[0] == "I") and not(combine):

                # Increment nt_query_index and insert_count by one, update delete_count to zero
                nt_query_index += 1
                insert_count += 1 
                delete_count = 0

            # Deletion
            elif (nt[0] == "D") and not(combine):

                # Increment delete_count and nt_ref_index by one, update insert_count to zero, set has_deletion to true
                nt_ref_index += 1
                delete_count += 1
                insert_count = 0
                has_deletion = True

                # Case 1: the reference codon is not entirely deleted and we are in first codon in 
                #         MDDMM|M... or MMDM|MM.. and must map to first in ABC|ABC
                if (first_alignment != None) and (prev_aa_shift % 3) == 0:
                    addToMap(int(nt_ref_index-1/3), last_alignment)

                # Case 2: the reference codon has been deleted
                elif first_alignment == None:
                    addToMap(int(nt_ref_index-1/3), '-')
                    # If we are in a frameshift
                    if (prev_aa_shift % 3) != 0:  
                        inbetween = '-'

                prev_aa_shift = None
                first_alignment = None

            # Not indel
            else: #nt[0] == "M"

                # Increment nt_query_index and nt_ref_index by one, update insert_count and delete_count to zero
                nt_query_index += 1
                nt_ref_index += 1
                insert_count = 0
                delete_count = 0

                # If the first base in the reference codon was deleted
                if first_alignment == None: first_alignment = aa_query_index

                # Would happen for scenarios where first two alignment pairs in nt_alignment_map are insertions 
                # and third alignemnt pair is a match to the last base pair in the reference codon
                if prev_aa_shift == None: prev_aa_shift = nt[4]


                # Case by case scenario of what could be true for the reference codon:
                
                # Case 1: no indel or there is a deletion for each insertion (MMM or DMIM)
                if ((first_alignment == last_alignment) and
                    ((prev_aa_shift % 3) == 0) and
                    ((nt[4] % 3) == 0) and 
                    (inbetween == None)):
                        addToMap(int(nt_ref_index-1/3), last_alignment)
                # Case 2: there is has a deletion that started a frameshift (MDM..., DDM...)
                elif ((first_alignment == last_alignment) and
                      ((prev_aa_shift % 3) == 0) and
                      ((nt[4] % 3) != 0) and 
                      (inbetween == None)):
                        addToMap(int(nt_ref_index-1/3), last_alignment)
                # Case 3: query codon has a three+ consecutive base pair deletion that span previous and 
                #         current reference codons and we are no longer in frameshift (MDDDMM, MIM|MMDDDDM)
                elif ((first_alignment == last_alignment) and
                      ((prev_aa_shift % 3) != 0) and
                      ((nt[4] % 3) == 0) and 
                      (inbetween == '-')):
                        addToMapInbetween(int(nt_ref_index-1/3), last_alignment)
                        inbetween = None
                # Case 4: we just ended a frameshift either through the deletion of the last base pair in the previous reference codon  
                #         or through the deletion of at least one of the first two base pairs in the current reference codon 
                #         (MIM|MMM|MDMM, MDMM|DDMMM, MDMDDM, MIM|MMM|MMDM)
                elif ((first_alignment == last_alignment) and
                      ((nt[4] % 3) == 0) and 
                      (inbetween != '-')):
                        addToMapDeletion(int(nt_ref_index-1/3), last_alignment)
                        inbetween = None
                # Case 5: we didn't have a frameshift at the start and we had x insertions in a single reference codon, 
                #         where x is a multiple of 3 (IMM|IIM, MII|IMM, III|MMM)
                elif ((first_alignment != last_alignment) and
                      ((prev_aa_shift % 3) == 0) and
                      ((nt[4] % 3) == 0) and 
                      (inbetween == None)):
                        aa_alignment_map.update({int(nt_ref_index-1/3):(first_alignment, last_alignment)})
                # Case 6: we were already in a frameshift, and we ended it through insertions (IMM|MMM|IIM)
                elif ((first_alignment != last_alignment) and
                      ((prev_aa_shift % 3) != 0) and
                      ((nt[4] % 3) == 0) and 
                      (inbetween != None)):
                        addToMapInbetween(int(nt_ref_index-1/3), last_alignment)
                        inbetween = None  
                # Case 7: we started a frameshift through an insertion (MII|MM... or MIM|M...)
                elif ((first_alignment != last_alignment) and
                      ((prev_aa_shift % 3) == 0)):
                        addToMapInbetween(int(nt_ref_index-1/3), first_alignment)
                        inbetween = None  
                # Case 8: we started the frameshift through a three consecutive base pair deletion that span previous and 
                #         current reference codons which was followed by a three consecutive base pair insertion (MDDDMI|IIM)
                elif ((first_alignment != last_alignment) and 
                      ((prev_aa_shift % 3) == 0) and
                      ((nt[4] % 3) == 0) and
                      (inbetween == '-')):
                        aa_alignment_map.update({int(nt_ref_index-1/3):(inbetween, first_alignment, last_alignment)})     
                        inbetween = None
                else:
                    raise(RuntimeError("Not valid case in aaAlignment - contact developer to fix this"))
                prev_aa_shift = None
                first_alignment = None

        aa_query_index = int(nt_query_index / 3)
        last_alignment = aa_query_index

    return aa_alignment_map    