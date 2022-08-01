from SNP_Verification_Tools.Gene import Gene
from SNP_Verification_Tools.SNP import SNP
from SNP_Verification_Tools import disregard
from SNP_Verification_Processes.DifferenceCheck import DifferenceCheck
from SNP_Verification_Processes.MappingQueryToReference import MapQueryToReference
from SNP_Verification_Processes.NonsenseCheck import NonsenseCheck
from SNP_Verification_Processes.IntrinsicCheck import IntrinsicCheck
from SNP_Verification_Processes.MisInDelCheck import MisInDelCheck
from SNP_Verification_Processes.nTupleCheck import nTupleCheck
import pysam

def verify(read, gene, argInfoDict, intrinsicInfoDict, frameshiftInfoDict, meg_3180InfoDict, mt_and_wt):
    name = gene.getName()
    rRna = gene.rRna()
    if not(DifferenceCheck(name, read, rRna, argInfoDict, frameshiftInfoDict)): 
        return False
    seqOfInterest, mapOfInterest = MapQueryToReference(rRna, read)
    if not(rRna):
        nonsense = NonsenseCheck(seqOfInterest, mapOfInterest, gene, read, name, argInfoDict)
        if nonsense != None: return nonsense
    if (IntrinsicCheck(mapOfInterest, seqOfInterest, gene, name, read, intrinsicInfoDict, argInfoDict)):
        return True
    elif (MisInDelCheck(gene, name, mapOfInterest, seqOfInterest, mt_and_wt, argInfoDict, meg_3180InfoDict)):
        return True
    elif (nTupleCheck(mapOfInterest, gene, name, seqOfInterest, mt_and_wt, argInfoDict, meg_3180InfoDict)):
        return True
    else: disregard(name, argInfoDict)
    return False