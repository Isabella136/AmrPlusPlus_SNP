from SNP_Verification_Processes import resistant, disregard
from SNP_Verification_Tools import Gene
from SNP_Verification_Tools import SNP

def nTupleCheck(aa_alignment_map, gene, name, aaQuerySequence, mt_and_wt, argInfoDict):
    SNPInfo = gene.condensedMultInfo()
    if (len(SNPInfo) == 0):
        return False
    else:
        resBool = True
        for snpMult in SNPInfo:
            wtPresent = False       # Can only be true if mt_and_wt is false
            for snp in snpMult:
                #codons that count for one deletion mutation but aren't fully deleted in read
                #if >1, snp is disregarded
                codonNotFullyDeleted = 0 
                if (aa_alignment_map.get(snp[1]-1,False)) == False:
                    resBool = False
                    break
                for mt in snp[2]: 
                    delMut = 0
                    misMut = 0
                    actualResBool = False
                    for queryIndex in tuple(aa_alignment_map[snp[1]-1]):
                        if queryIndex == '-':
                            if mt == queryIndex:
                                delMut += 1
                                actualResBool = True
                                continue
                        elif mt == aaQuerySequence[queryIndex]:
                            actualResBool = True
                            if mt_and_wt:
                                break
                        elif (snp[0] == aaQuerySequence[queryIndex]) and not(mt_and_wt):
                            wtPresent = True
                            actualResBool = False
                            resBool = False
                            break
                        else:
                            misMut += 1
                            resBool = False
                    if (delMut > 0) & (misMut > 0):
                        codonNotFullyDeleted += 1
                    if actualResBool: 
                        resBool = actualResBool
                        break
                    if wtPresent:
                        break
                if not(resBool): break
            if resBool: 
                if codonNotFullyDeleted < 2:
                    resistant(name, 1, argInfoDict)
                    return True
                else:
                    resBool = False
    return False