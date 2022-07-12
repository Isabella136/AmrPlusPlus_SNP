from SNP_Verification_Processes import disregard

def DifferenceCheck(name, read, rRna, argInfoDict):
    if not(rRna):
        cigarOpCount = read.get_cigar_stats()[0].tolist()
        #if (insertions - deletions) %3 != 0, disregard
        if ((cigarOpCount[1] - cigarOpCount[2]) % 3) != 0: 
            disregard(name, argInfoDict) 
            return False
    return True
