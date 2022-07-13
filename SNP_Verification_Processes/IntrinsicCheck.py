from SNP_Verification_Tools import resistant
from SNP_Verification_Tools import Gene
from SNP_Verification_Tools import SNP

intrinsicArg = ["MEG_1628", "MEG_3979", "MEG_3983", "MEG_4279", "MEG_4280", "MEG_4281", "MEG_4282", "MEG_6092"]

def IntrinsicCheck(mapOfInterest, seqOfInterest, gene, name, read, intrinsicInfoDict, argInfoDict):
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
            else:
                missingList = listInMissingList
        messageType = None
        if all: messageType = "All"
        intrinsicResistant(name, read.query_name, missingList, messageType, intrinsicInfoDict)
        if missingList == None:
            resistant(name, 1, argInfoDict)
            return True
    elif name in intrinsicArg:
        intrinsicResistant(name, read.query_name, None, "NA", intrinsicInfoDict)
    return False

def intrinsicResistant(name, queryName, missing, messageType, intrinsicInfoDict):
    if name not in intrinsicInfoDict:
        intrinsicInfoDict.update({name:{}})
    if missing == None:
        if (messageType == None):
            if "Some" not in intrinsicInfoDict[name]:
                intrinsicInfoDict[name].update({"Some":list()})
            intrinsicInfoDict[name]["Some"].append(queryName)
        elif messageType == "All":
            if "All" not in intrinsicInfoDict[name]:
                intrinsicInfoDict[name].update({"All":list()})
            intrinsicInfoDict[name]["All"].append(queryName)
        else:
            if "NA" not in intrinsicInfoDict[name]:
                intrinsicInfoDict[name].update({"NA":list()})
            intrinsicInfoDict[name]["NA"].append(queryName)
    else:
        if "Mutant" not in intrinsicInfoDict[name]:
            intrinsicInfoDict[name].update({"Mutant":list()})
        intrinsicInfoDict[name]["Mutant"].append(queryName)

