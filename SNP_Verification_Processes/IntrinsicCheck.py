from SNP_Verification_Tools.Gene import Gene
from SNP_Verification_Tools.SNP import SNP
from SNP_Verification_Tools.InDel import InDel

def IntrinsicCheck(read, gene, mapOfInterest, seqOfInterest, config):
    if gene.getGeneTag() != 'I':
        return True
    begin = list(mapOfInterest.keys())[0]+1
    end = list(mapOfInterest.keys())[-1]+1
    firstInfo = gene.getFirstMustBetweenParams(begin, end)
    missingList = []
    if firstInfo != None:
        must, all = firstInfo
        hasAllMust = True
        while(must.getPos() < end):
            hasMust = False
            firstAndLastTupleIndex = [0,len(mapOfInterest[must.getPos()-1]) - 1]
            tupleIndex = -1
            for queryIndex in mapOfInterest[must.getPos()-1]:
                tupleIndex += 1
                if gene.rRna():
                    if tupleIndex not in firstAndLastTupleIndex:
                        continue
                if (queryIndex == None) or (queryIndex == '-'):
                    if not(config.getboolean('SETTINGS', 'MT_AND_WT')):
                        hasMust = False
                        break
                elif must.getWt() == seqOfInterest[queryIndex]:
                    hasMust = True
                    if config.getboolean('SETTINGS', 'MT_AND_WT'):
                        break
                else:
                    if not(config.getboolean('SETTINGS', 'MT_AND_WT')):
                        hasMust = False
                        break
            if not(hasMust): 
                hasAllMust = False
                missingList.append(must.condensedInfo())
            must = must.getNext()
            if must == None:
                all = all & True
                break
        if must != None:
            all = False
        if hasAllMust:
            missingList = None
        messageType = "Some"
        if missingList != None:
            messageType = "Mutant"
        elif all:
            messageType = "All"
        gene.addDetails(read, messageType)
        return True
    gene.addDetails(read, "NA")
    return True
    

