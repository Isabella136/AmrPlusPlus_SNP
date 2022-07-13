from SNP_Verification_Tools import disregard

def DifferenceCheck(name, read, rRna, argInfoDict, frameshiftInfoDict):
    if not(rRna):
        cigarOpCount = read.get_cigar_stats()[0].tolist()
        #if (insertions - deletions) %3 != 0, disregard
        if ((cigarOpCount[1] - cigarOpCount[2]) % 3) != 0: 
            disregard(name, argInfoDict) 
            addRead(name, read.query_name, frameshiftInfoDict)
            return False
    return True

def addRead(name, queryName, frameshiftInfoDict):
    if name not in frameshiftInfoDict:
        frameshiftInfoDict.update({name:list()})
    frameshiftInfoDict[name].append(queryName)