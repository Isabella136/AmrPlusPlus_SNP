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
# .
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
#   first_alignment:
#   third_alignment:
#   inbetween:                          If in frameshift, inbetween should either be '-' or previous query codon
#   has_deletion:                       True if, and only if, current reference codon has a deletion
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
    third_alignment = 0
    inbetween = None
    has_deletion = False

    def addToMapDeletion(index, toAdd):                         #maps queryCodon in inbetween to previous referenceCodon
        aa = aa_alignment_map.get(index-1, False)               #currently, toAdd and inbetween are the same
        if aa == False:
            aa_alignment_map.update({index-1:(inbetween, )})
        else:
            temp = list(aa)
            temp.append(inbetween)
            aa = tuple(temp)
            aa_alignment_map.update({index-1:aa})
        addToMap(index, toAdd)

    def addToMapInbetween(index, toAdd):                        #adds both queryCodon and inbetweenQueryCodon to same referenceCodon
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

    def addToMap(index, toAdd):                                 #maps queryCodon to referenceCodon
        aa = aa_alignment_map.get(index, False)
        if aa == False:
            aa_alignment_map.update({index:(toAdd,)})
        else:
            temp = list(aa)
            temp.append(toAdd)
            aa = tuple(temp)
            aa_alignment_map.update({index:aa})

    for nt in nt_alignment_map:
        # If previous alignment pair was last query base pair in codon
        if ((nt_query_index % 3) == 0) and (delete_count == 0): has_deletion = False

        # If we reached start of new reference codon
        if (nt_ref_index % 3) == 0:

            # Redefines prev_aa_shift because it is a new reference codon
            prev_aa_shift = nt[4]

            # Insertion
            if nt[0] == "I":

                # Increment insert_count and nt_query_index by one, update delete_count to zero
                insert_count += 1
                nt_query_index += 1
                delete_count = 0

                # If shift in alignment is now a multiple of three
                if (nt[3] % 3) == 0:  

                    # If there has been three consecutive insertions
                    if insert_count == 3:
                        # That didn't happend right after another inserted codon
                        if inbetween == None:
                            # Inserted codon mapped to previous reference codon
                            addToMap(int(nt_ref_index/3)-1, third_alignment)

                        # Defines inbetween as current query codon so that it gets also mapped to upcoming reference codon
                        inbetween = third_alignment

                    # If shift % 3 == 0 and insert_count is not 3, then frameshift is present and if of form #1I...1M2I or 2I...2M1I 
                    # Inbetween is previous query codon, and both the previous and current query codons are mapped to the previous reference codon
                    else:
                        addToMapInbetween(int(nt_ref_index/3)-1, third_alignment)
                        inbetween = None

                    # Resets insert_count
                    insert_count = 0
                
            # Deletion
            elif nt[0] == "D":

                # Increment delete_count and nt_ref_index by one, update insert_count to zero, set has_deletion to true
                nt_ref_index += 1
                insert_count = 0
                delete_count += 1
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

                # We started frameshift with two insertions, and continued with match eversince
                if (nt_query_index % 3) == 0:

                    # Case: we were not in a frameshift, but we started the query codon with two insertions
                    if (prev_aa_shift % 3) == 0: #??????? This used to be prev_aa_shift == 0
                        inbetween = None
                    elif inbetween == None:
                        if not(has_deletion):                                        #second codon in IIM|MMM
                            addToMap(int(nt_ref_index/3)-1, third_alignment)           #maps to first in ABC|A...
                        if nt[3] < 0:
                            addToMap(int(nt_ref_index/3), third_alignment)             #first codon in MDMM|MM...
                        inbetween = third_alignment                                  #maps to last in ABC|ABC
                    else:
                        inbetween = third_alignment                                  #last codon in IIM|MMM|MMM
                first_alignment = aa_query_index
        elif (nt_ref_index % 3) == 1:
            if nt[0] == "I":
                delete_count = 0
                nt_query_index += 1
                if (nt_query_index % 3) == 0: #1I2M2M1I
                    if prev_aa_shift == 0:                                            #codon MII
                        inbetween = None
                    elif inbetween == None:                                         #second codon in IMM|MMI
                        addToMap(int(nt_ref_index/3)-1, third_alignment)               #maps to first codon in ABC|A...
                        inbetween = third_alignment
                    else:                                                           #last codon in IMM|MMM|MMI
                        inbetween = third_alignment
            elif nt[0] == "D":
                nt_ref_index += 1
                delete_count += 1
                if delete_count == 3:                                                #for codon MMDDDM, '-' mapped to first in ABC|ABC
                    addToMap(int(nt_ref_index/3)-1, '-')                              #codon was previously mapped in else/elif D
                    inbetween = '-'
                has_deletion = True
            else: #nt[0] == "M"
                delete_count = 0
                nt_query_index += 1
                nt_ref_index += 1
                if first_alignment == None:
                    first_alignment = aa_query_index

                # Would happen for scenarios where first alignment pair in nt_alignment_map is an insertion 
                # and second alignemnt pair is a match to the second base pair in the reference codon
                if prev_aa_shift == None:
                    prev_aa_shift = nt[4]
                if (nt_query_index % 3) == 0: #1I2M3M
                    if (prev_aa_shift == 0):
                        inbetween = None
                    elif inbetween == None:
                        if not(has_deletion):                                        #second codon in IMM|MMM
                            addToMap(int(nt_ref_index/3)-1, third_alignment)           #maps to first in ABC|AB...
                        if nt[3] < 0:
                            addToMap(int(nt_ref_index/3), third_alignment)             #first codon in MDDMM|M...
                        inbetween = third_alignment                                  #maps to last in ABC|ABC
                    else:
                        inbetween = third_alignment
        else: #(nt_ref_index % 3) == 2
            if nt[0] == "I":
                delete_count = 0
                nt_query_index += 1
            elif nt[0] == "D":
                has_deletion = True
                delete_count += 1
                nt_ref_index += 1
                if first_alignment != None:
                    if (prev_aa_shift % 3) == 0: #1M2D, 2M1D                          #first codon in MDDMM|M... or just MDDDMM
                        addToMap(int(nt[2]/3), third_alignment)                      #maps to first in ABC|ABC
                else: #3D
                    addToMap(int(nt[2]/3), '-')                                     #'-' maps to ABC
                    if (prev_aa_shift % 3) != 0:  
                        inbetween = '-'
                prev_aa_shift = None
                first_alignment = None
            else: #nt[0] == "M"
                delete_count = 0
                nt_query_index += 1
                nt_ref_index += 1
                if first_alignment == None:
                    first_alignment = aa_query_index

                # Would happen for scenarios where first two alignment pairs in nt_alignment_map are insertions 
                # and third alignemnt pair is a match to the last base pair in the reference codon
                if prev_aa_shift == None:
                    prev_aa_shift = nt[4]
                if first_alignment == third_alignment:
                    if ((prev_aa_shift % 3) == 0) | ((nt[4] % 3) == 0): #3M, 1D1I2M, 1D/2I...2D1M, 2D/1I...1D2M, 2D3M, 1D3M
                        if inbetween != None:
                            if has_deletion:
                                if inbetween != '-':                                #last codon in MIM|MMM|MDMM or in MDMM|DDMMM
                                    inbetween = third_alignment                      #maps to last two in ABC|ABC|ABC
                                    addToMapDeletion(int(nt[2]/3), third_alignment)
                                else:                                               #codon such as MDDDMM and '-'
                                    addToMapInbetween(int(nt[2]/3), third_alignment) #maps to second in ABC|ABC

                            else:                                                   #???????????????
                                addToMapInbetween(int(nt[2]/3), third_alignment)     #???????????????
                            inbetween = None
                        else:                                                       #simple codon MMM or DMIM or MDMM
                            addToMap(int(nt[2]/3), third_alignment)                  #maps to ABC
                elif (prev_aa_shift % 3) == 0:
                    if (nt[4] % 3) == 0: #3Mand3xI
                        if inbetween != None: #inbetween == "-"                                                     #both codon in MMI|DDDIIM and '-'
                            aa_alignment_map.update({int(nt[2]/3):(inbetween, first_alignment, third_alignment)})     #maps to last codon in ABC|ABC 
                            inbetween = None
                        else:                                                                                       #codon such as IMM|IIM
                            aa_alignment_map.update({int(nt[2]/3):(first_alignment, third_alignment)})                #maps to ABC
                    else: #1/2Iand3M
                        addToMapInbetween(int(nt[2]/3), first_alignment)                                             #codon such as MII|MM...
                        inbetween = None                                                                            #maps to ABC
                elif (nt[4] % 3) == 0: #2D/1I...1M1D1M
                    if inbetween != None:
                        if has_deletion:
                            if inbetween != '-':                                    #last codon in MIM|MMM|MMDM
                                inbetween = third_alignment                          #maps to last two in ABC|ABC|ABC
                                addToMapDeletion(int(nt[2]/3), third_alignment)
                            else:                                                   #last codon in MIM|MMDDDDM and '-'
                                addToMapInbetween(int(nt[2]/3), third_alignment)     #maps to last codon in ABC|ABC
                        else:
                            addToMapInbetween(int(nt[2]/3), third_alignment)         #codon MDDDMM and '-'
                        inbetween = None                                            #maps to last codon in ABC|ABC

                    else:                                                           #???????????????
                        addToMap(int(nt[2]/3), third_alignment)                      #???????????????
                prev_aa_shift = None
                first_alignment = None

        aa_query_index = int(nt_query_index / 3)
        third_alignment = aa_query_index

    return aa_alignment_map    

