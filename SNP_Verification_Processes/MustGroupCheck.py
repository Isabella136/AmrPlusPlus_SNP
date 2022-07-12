from SNP_Verification_Processes import resistant, disregard
from SNP_Verification_Tools import Gene
from SNP_Verification_Tools import SNP

mustArg = ["MEG_1628", "MEG_3979", "MEG_3983", "MEG_4279", "MEG_4280", "MEG_4281", "MEG_4282", "MEG_6092", "MEG_6093"]

def MustGroupCheck(mapOfInterest, seqOfInterest, gene, name, read, mustGroupInfoDict, argInfoDict):
    begin = list(mapOfInterest.keys())[0]+1
    end = list(mapOfInterest.keys())[-1]+1
    mustList = gene.getFirstMustBetweenParams(begin, end)
    missingList = []
    if mustList != None:
        for must, all in mustList:
            listInMissingList = []
            hasAllMust = True
            while(must.getPos() < end):
                hasMust = False
                for queryIndex in mapOfInterest[must.getPos()-1]:
                    if queryIndex == None:
                        continue
                    if must.getWt() == seqOfInterest[queryIndex]:
                        hasMust = True
                        continue
                if not(hasMust): 
                    hasAllMust = False
                    listInMissingList.append(must.condensedInfo())
                must = must.getNext()
                if must == None:
                    all = all & True
                    break
            if must != None:
                all = False
            if hasAllMust:
                missingList = None
                break
            elif name == "MEG_6093":
                missingList.append(listInMissingList)
            else:
                missingList = listInMissingList
        messageType = None
        if name == "MEG_6093": 
            messageType = "MEG_6093"
        if all and (missingList == None): messageType = "All"
        mustResistant(name, read.query_name, missingList, messageType, gene.aaOrNu(), mustGroupInfoDict)
        if missingList == None:
            resistant(name, 1, argInfoDict)
            return True
    elif name in mustArg:
        mustResistant(name, read.query_name, None, "NA", gene.aaOrNu(), mustGroupInfoDict)
    return False

def mustResistant(name, queryName, missing, messageType, aaOrNu, mustGroupInfoDict):
    if name not in mustGroupInfoDict:
        mustGroupInfoDict.update({name:list()})
    if missing == None:
        if (messageType == None) or (messageType == "MEG_6093"):
            mustGroupInfoDict[name].append("Based on the sequence given, " + queryName + " contains the " + aaOrNu + " required for intrinsic resistance")
        elif messageType == "All":
            mustGroupInfoDict[name].append(queryName + " contains ALL " + aaOrNu + " required for intrinsic resistance")
        else:
            mustGroupInfoDict[name].append("The sequence given for " + queryName + " does not include the position where the " + aaOrNu + " required for intrinsic resistance are located")
    elif messageType == "MEG_6093":
        mustGroupInfoDict[name].append("The sequence given for " + queryName + " does not have the following amino acids required for resistance: " + str(missing[0]) + "; nor does it have this alternative group of amino acids that can also induce resistance: " + str(missing[1]))
    else:
        mustGroupInfoDict[name].append("The sequence given for " + queryName + " does not have the following " + aaOrNu + " required for intrinsic resistance:" + str(missing))

