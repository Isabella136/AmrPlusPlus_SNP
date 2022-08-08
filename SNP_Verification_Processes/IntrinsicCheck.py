from SNP_Verification_Tools import resistant
from SNP_Verification_Tools import Gene
from SNP_Verification_Tools import SNP
from SNP_Verification_Tools import intrinsicInfoDict, mt_and_wt

intrinsicArg = ["MEG_1628", "MEG_2694", "MEG_3983", "MEG_4092", "MEG_4093",
                "MEG_4096", "MEG_4279", "MEG_4280", "MEG_4281", "MEG_4282", 
                "MEG_6092", "MEG_7196", "MEG_7306"]

def IntrinsicCheck(read, gene, mapOfInterest, seqOfInterest):
    if gene.getName() == "MEG_4282":
        test = True
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
                    if not(mt_and_wt):
                        hasMust = False
                        break
                elif must.getWt() == seqOfInterest[queryIndex]:
                    hasMust = True
                    if mt_and_wt:
                        break
                else:
                    if not(mt_and_wt):
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
        messageType = None
        if missingList != None:
            messageType = "Mutant"
        elif all: messageType = "All"
        if (gene.getName()== "MEG_3979") and not(all) :
            return messageType
        intrinsicResistant(gene.getName(), read.query_name, messageType)
        return True
    elif gene.getName()== "MEG_3979" :
        return "NA"
    elif gene.getName() in intrinsicArg:
        intrinsicResistant(gene.getName(), read.query_name, "NA")
        return True
    return False

def intrinsicResistant(name, queryName, messageType):
    if name not in intrinsicInfoDict:
        intrinsicInfoDict.update({name:{}})
    if (messageType == None):
        if "Some" not in intrinsicInfoDict[name]:
            intrinsicInfoDict[name].update({"Some":list()})
        intrinsicInfoDict[name]["Some"].append(queryName)
    elif messageType == "All":
        if "All" not in intrinsicInfoDict[name]:
            intrinsicInfoDict[name].update({"All":list()})
        intrinsicInfoDict[name]["All"].append(queryName)
    elif messageType == "NA":
        if "NA" not in intrinsicInfoDict[name]:
            intrinsicInfoDict[name].update({"NA":list()})
        intrinsicInfoDict[name]["NA"].append(queryName)
    elif messageType == "Mutant":
        if "Mutant" not in intrinsicInfoDict[name]:
            intrinsicInfoDict[name].update({"Mutant":list()})
        intrinsicInfoDict[name]["Mutant"].append(queryName)
    else:   #if messageType == "Aquired"
        if "Aquired" not in intrinsicInfoDict[name]:
            intrinsicInfoDict[name].update({"Aquired":list()})
        intrinsicInfoDict[name]["Aquired"].append(queryName)
    

