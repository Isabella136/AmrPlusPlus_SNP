from SNP_Verification_Tools.Gene import Gene
from SNP_Verification_Tools.SNP import SNP
from SNP_Verification_Tools.InDel import InDel

def NonsenseCheck(read, gene, map_of_interest, seq_of_interest):
    stop_location = seq_of_interest.find('*') 
    # If nonsense mutation is present and it isn't located at the end of the reference sequence
    if (stop_location != -1) and ((stop_location < len(seq_of_interest)-1) or (read.reference_end < gene.ntSeqLength())): 
        return verifyNonsense(read, gene, stop_location, map_of_interest)
    return True

def verifyNonsense(read, gene, stop_location, map_of_interest):

    def findNonsenseReferenceLocation():
        # Get list of all query codons that aligned to reference
        aligned_query_codon = list(map_of_interest.values())
        reference_stop_location = -1
        last_query_index = aligned_query_codon[0][0]-1
        for i, query_tuple in enumerate(aligned_query_codon):
            for query_index in query_tuple:
                # Skip deleted codons
                if query_index == '-':
                    continue

                # query_index was found
                if stop_location == query_index:
                    reference_stop_location = list(map_of_interest.keys())[i]
                    break

                # the nonsense is within a frameshift
                elif stop_location < query_index:
                    reference_stop_location = list(map_of_interest.keys())[i-1] + stop_location - last_query_index
                    break
                last_query_index = query_index
            # the matching reference codon index was found
            if reference_stop_location != -1:
                break

        # the nonsense is within a frameshift that extends past th alignment
        if reference_stop_location == -1:
            reference_stop_location = list(map_of_interest.keys())[i-1] + stop_location - last_query_index
        gene.addDetails(read, "Newly found nonsense: " + str(reference_stop_location+1))

    # Retrieve information about nonsense mutation variants of gene
    SNPInfo = gene.condensedNonInfo()

    # If no information is present and gene is not of type frameshift, find location in reference and return false
    if (len(SNPInfo) == 0) and (gene.getGeneTag() != 'F'):
        findNonsenseReferenceLocation()
        return False
    else:
        for snp in SNPInfo:
            reference_index = map_of_interest.get(snp[1]-1,False)
            if reference_index != False:
                # If information on variant matches current alignement, return true
                if reference_index[0] == stop_location:
                    gene.addDetails(read, snp)
                    return True
        
        # If no matches were found, find location in reference and return true if gene is of type frameshift, false otherwise
        findNonsenseReferenceLocation()
        return (gene.getGeneTag() == 'F')