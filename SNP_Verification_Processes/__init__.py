from SNP_Verification_Tools.Gene import Gene
from SNP_Verification_Tools.SNP import SNP
from SNP_Verification_Processes.DifferenceCheck import DifferenceCheck
from SNP_Verification_Processes.MappingQueryToReference import MapQueryToReference
from SNP_Verification_Processes.NonsenseCheck import NonsenseCheck
from SNP_Verification_Processes.MustGroupCheck import MustGroupCheck
from SNP_Verification_Processes.MissenseAndDeletionCheck import MissenseAndDeletionCheck
from SNP_Verification_Processes.nTupleCheck import nTupleCheck
import pysam

def resistant(name, increment, argInfoDict):
    argInfo = argInfoDict.get(name, False)
    if (argInfo == False):
        argInfoDict.update({name:(increment, 1)})
    else:
        temp = list(argInfo)
        temp[0] += increment
        temp[1] += 1
        argInfo = tuple(temp)
        argInfoDict.update({name:argInfo})

def disregard(name, argInfoDict):
    resistant(name, 0, argInfoDict)

def verify(read, gene, argInfoDict, mustGroupInfoDict, mt_and_wt):
    name = gene.getName()
    rRna = gene.rRna()
    if not(DifferenceCheck(name, read, rRna, argInfoDict)): 
        return False
    seqOfInterest, mapOfInterest = MapQueryToReference(rRna, read)
    if not(rRna):
        nonsense = NonsenseCheck(seqOfInterest, mapOfInterest, gene, read, name, argInfoDict)
        if nonsense != None: return nonsense
    if (MustGroupCheck(mapOfInterest, seqOfInterest, gene, name, read, mustGroupInfoDict, argInfoDict)):
        return True
    elif (MissenseAndDeletionCheck(gene, name, mapOfInterest, seqOfInterest, mt_and_wt, argInfoDict)):
        return True
    elif (nTupleCheck(mapOfInterest, gene, name, seqOfInterest, mt_and_wt, argInfoDict)):
        return True
    else: disregard(name, argInfoDict)
    return False