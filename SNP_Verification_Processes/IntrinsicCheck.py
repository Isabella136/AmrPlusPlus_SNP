from SNP_Verification_Tools.Gene import Gene
from SNP_Verification_Tools.SNP import SNP
from SNP_Verification_Tools.InDel import InDel

def IntrinsicCheck(read, gene, map_of_interest, seq_of_interest, config):
    # Find beginning and ending positions of aligned portion reference sequence
    begin = list(map_of_interest.keys())[0]+1
    end = list(map_of_interest.keys())[-1]+1

    # Returns None if no 'must' amino/nucleic acids are in the aligned portion of the reference sequence
    # Else, returns tuple of first 'must' encountered in the aligned portion of the reference sequence and a boolean
    # indicating whether it coincides with the first 'must' of the overall reference sequence
    first_info = gene.getFirstMustBetweenParams(begin, end)

    # If no 'must' amino/nucleic acids are in the aligned portion of the reference sequence, stops here
    if first_info == None:
        gene.addDetails(read, "NA")
        return True
    
    # current_must:                     updated at the end of each iteration of while loop
    # all:                              indicates whether first 'must' encountered in the aligned 
    #                                   portion is first 'must' of the overall reference sequence
    # aligned_portion_has_all_musts:    is true if all 'must' amino/nucleic acids of the aligned 
    #                                   portion encountered so far are present in query
    current_must, all = first_info
    aligned_portion_has_all_musts = True

    while((current_must != None) and (current_must.getPos() <= end)):
        # query_has_must:               becomes True if 'must' amino/nucleic acid is found in query
        # first_and_last_tuple_index:   see line 38 for comments
        query_has_must = False
        first_and_last_tuple_index = [0,len(map_of_interest[current_must.getPos()-1]) - 1]

        for tuple_index, query_index in enumerate(map_of_interest[current_must.getPos()-1]):

            # Because rRNA can have multiple query nucleotides mapped to one reference nucleotide in the
            # case of 2+ bp insertion, only the first and last query nucleotides are taken into account
            if gene.rRna() and tuple_index not in first_and_last_tuple_index: continue

            # In scenarios where a wild-type and a deleted amino acid can't simultaneously exist 
            # in resistant gene, query is not considered to have current 'must' amino/nucleic acid
            if query_index in [None,'-'] and not(config.getboolean('SETTINGS', 'MT_AND_WT')):
                query_has_must = False
                break

            # 'must' amino/nucleic acid is found in query
            elif current_must.getWt() == seq_of_interest[query_index]:
                query_has_must = True
                # In scenarios where a mutant and a wild-type amino/nucleic acid can simultaneously exist
                # in resistant gene, no need to further check rest of amino/nucleic acid mapped to current 'must'
                if config.getboolean('SETTINGS', 'MT_AND_WT'): break

            # In scenarios where a mutant and a wild-type amino/nucleic acid can't simultaneously exist 
            # in resistant gene, query is not considered to have current 'must' amino/nucleic acid
            elif not(config.getboolean('SETTINGS', 'MT_AND_WT')):
                query_has_must = False
                break

        # Can't be considered for intrinsic resistance if mutation is present at a 'must' amino/nucleic acid
        if not(query_has_must): 
            aligned_portion_has_all_musts = False
            break
        
        # Go to next 'must' amino/nucleic acid
        current_must = current_must.getNext()

    # all:  now indicates whether first 'must' encountered in the aligned portion is first 'must' of the overall reference 
    #       sequence and last 'must' encountered in the aligned portion is last 'must' of the overall reference sequence
    all = all and (current_must == None)

    if not(aligned_portion_has_all_musts):
        gene.addDetails(read, "Mutant")
    elif all:
        gene.addDetails(read, "All")
    else:
        gene.addDetails(read, "Some")
    return True


