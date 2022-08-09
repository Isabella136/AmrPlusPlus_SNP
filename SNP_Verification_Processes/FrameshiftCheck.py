from SNP_Verification_Tools import dnaTranslate
from SNP_Verification_Tools import detailed

def FrameshiftCheck(read, gene, rRna):
    if not(rRna):
        cigarOpCount = read.get_cigar_stats()[0].tolist()
        frameshiftInfo = gene.getFrameshiftInfo()
        #if (insertions - deletions) %3 != 0, susceptible unless if gene has frameshift info
        if (((cigarOpCount[1] - cigarOpCount[2]) % 3) != 0) and (frameshiftInfo == None): 
            gene.addToOutputInfo(9 if gene.getGeneTag() == 'I' else 10)
            if detailed:
                gene.addDetails(read.query_name, 'FS end', True)
            return False
        elif frameshiftInfo != None:
            if gene.getGeneTag() == 'S':                                            #MEG_6094 can have a C insertion at res 531 if followed by frameshift suppression during translation
                querySequence = read.query_sequence
                aligned_pairs = read.get_aligned_pairs()[read.cigartuples[0][1] if read.cigartuples[0][0] == 4 else 0:len(querySequence)-read.cigartuples[-1][1] if read.cigartuples[-1][0] == 4 else len(querySequence)]
                shiftCount = 0
                inCodon531 = False
                hasCinsertion = False
                hasDeletionAfterCinsertion = False
                lastBeforeFull = (3 - (aligned_pairs[0][1]) % 3) % 3 - 1            #Last nt index before first full codon
                residue531To534 = ""
                valid = True
                for pair in aligned_pairs:
                    if pair[0] == None:                                             #If deletion
                        if hasCinsertion and not(hasDeletionAfterCinsertion):       #If first deletion since C insertion at codon 531
                            hasDeletionAfterCinsertion = True
                        else:
                            shiftCount -= 1
                    elif pair[1] == None:                                           #If insertion
                        if inCodon531 and ((shiftCount%3)!=0):                      #If in codon 531 and not in frameshift
                            if not(hasCinsertion):
                                shiftCount -= 1
                            hasCinsertion = True
                        shiftCount += 1
                    if pair[1] == 1590:
                        inCodon531 = True
                    if pair[1] == 1593:
                        inCodon531 = False
                    if pair[1] == (lastBeforeFull + 3):
                        if inCodon531 or (len(residue531To534) >= 0 and len(residue531To534) < 4):
                            residue531To534 += dnaTranslate(querySequence[lastBeforeFull+1:pair[1]+1])
                        lastBeforeFull += 3
                if (shiftCount % 3) == 2:
                    if not(hasCinsertion) or hasDeletionAfterCinsertion:
                        valid = False
                elif (shiftCount % 3) == 1:
                    valid = False
                if not(valid):
                    gene.addToOutputInfo(10)
                    if detailed:
                        gene.addDetails(read.query_name, 'FS end', True)
                    return False
                elif hasCinsertion:
                    if residue531To534 == "SRTR"[0:len(residue531To534)]:           #Must have those residues due to insertion
                        gene.addToOutputInfo(11)
                        if detailed:
                            gene.addDetails(read.query_name, 'C insert', True)
                        return "TFT"                                                #Read has info in gene.additionalInfo (T) but is still susceptible (F); has insert (T)
                    elif (shiftCount % 3) == 0:
                        gene.addToOutputInfo(10)
                        if detailed:
                            gene.addDetails(read.query_name, 'FS end', True)
                        return False

            else:
                if (((cigarOpCount[1] - cigarOpCount[2]) % 3) != 0):
                    gene.addToOutputInfo(2)
                    gene.addToOutputInfo(10)
                    if detailed:
                        gene.addDetails(read.query_name, 'FS end', True)
                    return True
                else:
                    indel = 0
                    for cigarTuple in read.cigartuples:
                        if (cigarTuple[0] in range(1,3)) and (cigarTuple[1] >= 12):
                            indel += 1
                    if indel > 0:
                        gene.addToOutputInfo(7)
                        if detailed:
                            gene.addDetails(read.query_name, "12+: " + str(indel))
                        return "TF"                                                 #Read has info in gene.additionalInfo (T) but is still susceptible (F)
    return "FF"

def addRead(name, queryName, frameshiftInfoDict, additionalInformation = None):
    if name not in frameshiftInfoDict:
        frameshiftInfoDict.update({name:list()})
    if (additionalInformation == None) or (additionalInformation == ''):
        frameshiftInfoDict[name].append(queryName)
    else:
        frameshiftInfoDict[name].append((queryName, additionalInformation))