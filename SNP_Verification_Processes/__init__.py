from SNP_Verification_Tools.Gene import Gene
from SNP_Verification_Tools.SNP import SNP
from SNP_Verification_Tools import disregard
from SNP_Verification_Processes.FrameshiftCheck import FrameshiftCheck
from SNP_Verification_Processes.MappingQueryToReference import MapQueryToReference
from SNP_Verification_Processes.NonsenseCheck import NonsenseCheck
from SNP_Verification_Processes.IntrinsicCheck import IntrinsicCheck
from SNP_Verification_Processes.MisInDelCheck import MisInDelCheck
from SNP_Verification_Processes.nTupleCheck import nTupleCheck
import pysam

def verify(read, gene, argInfoDict, intrinsicInfoDict, frameshiftInfoDict, meg_3180InfoDict, mt_and_wt):
    name = gene.getName()
    rRna = gene.rRna()
    if not(FrameshiftCheck(name, read, rRna, argInfoDict, frameshiftInfoDict)): 
        return False
    seqOfInterest, mapOfInterest = MapQueryToReference(rRna, read, name)
    if not(rRna):
        nonsense = NonsenseCheck(read, gene, mapOfInterest, seqOfInterest, argInfoDict, meg_3180InfoDict)
        if nonsense != None: return nonsense
    if (IntrinsicCheck(mapOfInterest, seqOfInterest, gene, name, read, intrinsicInfoDict, argInfoDict)):
        return True
    elif (MisInDelCheck(read, gene, mapOfInterest, seqOfInterest, mt_and_wt, argInfoDict, meg_3180InfoDict)):
        return True
    elif (nTupleCheck(gene, mapOfInterest, seqOfInterest, mt_and_wt, argInfoDict, meg_3180InfoDict)):
        return True
    else: 
        if name == "MEG_3180":
            if "susceptible" not in meg_3180InfoDict:
                meg_3180InfoDict["susceptible"] = 0
            meg_3180InfoDict["susceptible"] += 1
        disregard(name, argInfoDict)
    return False