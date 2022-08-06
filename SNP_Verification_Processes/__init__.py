from SNP_Verification_Tools.Gene import Gene
from SNP_Verification_Tools.SNP import SNP
from SNP_Verification_Tools import disregard
from SNP_Verification_Processes.FrameshiftCheck import FrameshiftCheck, addRead
from SNP_Verification_Processes.MappingQueryToReference import MapQueryToReference
from SNP_Verification_Processes.NonsenseCheck import NonsenseCheck
from SNP_Verification_Processes.IntrinsicCheck import IntrinsicCheck
from SNP_Verification_Processes.MisInDelCheck import MisInDelCheck
from SNP_Verification_Processes.nTupleCheck import nTupleCheck
from SNP_Verification import argInfoDict, resistantFrameshiftInfoDict, meg_6094InfoDict, meg_3180InfoDict
import pysam

def verify(read, gene):
    name = gene.getName()
    rRna = gene.rRna()
    frameshiftCheckResult = FrameshiftCheck(read, gene, rRna)
    if type(frameshiftCheckResult) == bool: 
        return frameshiftCheckResult
    seqOfInterest, mapOfInterest = MapQueryToReference(rRna, read, name)
    if not(rRna):
        nonsense = NonsenseCheck(read, gene, mapOfInterest, seqOfInterest)
        if nonsense != None: return nonsense
    if (IntrinsicCheck(read, gene, mapOfInterest, seqOfInterest)):
        return True
    elif (MisInDelCheck(read, gene, mapOfInterest, seqOfInterest)):
        return True
    elif (nTupleCheck(gene, mapOfInterest, seqOfInterest)):
        return True
    elif frameshiftCheckResult != "":
        addRead(read, name, read.query_name, resistantFrameshiftInfoDict, frameshiftCheckResult)
        return True
    else: 
        if name == "MEG_3180":
            if "susceptible" not in meg_3180InfoDict:
                meg_3180InfoDict["susceptible"] = 0
            meg_3180InfoDict["susceptible"] += 1
        elif name == "MEG_6094":
            addRead(name, read.query_name, meg_6094InfoDict, "No resistance-conferring mutations")
        elif frameshiftCheckResult == "":
            addRead(name, read.query_name, resistantFrameshiftInfoDict, "No resistance-conferring mutations")
        disregard(name, argInfoDict)
    return False