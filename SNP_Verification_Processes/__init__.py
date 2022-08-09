from SNP_Verification_Tools.Gene import Gene
from SNP_Verification_Tools.SNP import SNP
from SNP_Verification_Tools import disregard
from SNP_Verification_Processes.FrameshiftCheck import FrameshiftCheck, addRead
from SNP_Verification_Processes.MappingQueryToReference import MapQueryToReference
from SNP_Verification_Processes.NonsenseCheck import NonsenseCheck
from SNP_Verification_Processes.IntrinsicCheck import IntrinsicCheck, intrinsicResistant
from SNP_Verification_Processes.MisInDelCheck import MisInDelCheck
from SNP_Verification_Processes.nTupleCheck import nTupleCheck
from SNP_Verification_Tools import argInfoDict, resistantFrameshiftInfoDict, meg_6094InfoDict, meg_3180InfoDict
import pysam

def verify(read, gene):
    gene.addToOutputInfo(0)
    name = gene.getName()
    rRna = gene.rRna()
    checkResult = FrameshiftCheck(read, gene, rRna)
    if type(checkResult) == bool: 
        return checkResult
    seqOfInterest, mapOfInterest = MapQueryToReference(rRna, read, name)
    if not(rRna):
        nonsense = NonsenseCheck(read, gene, mapOfInterest, seqOfInterest, checkResult)
        if nonsense != None: return nonsense
    intrinsic = IntrinsicCheck(read, gene, mapOfInterest, seqOfInterest)
    if (intrinsic == True):
        return True
    elif (MisInDelCheck(read, gene, mapOfInterest, seqOfInterest)):
        return True
    elif (nTupleCheck(read, gene, mapOfInterest, seqOfInterest)):
        return True
    elif (checkResult != "") and (checkResult != None):
        addRead(name, read.query_name, resistantFrameshiftInfoDict, checkResult)
        return True
    else: 
        if name == "MEG_3180":
            if "susceptible" not in meg_3180InfoDict:
                meg_3180InfoDict["susceptible"] = 0
            meg_3180InfoDict["susceptible"] += 1
        elif name == "MEG_6094":
            addRead(name, read.query_name, meg_6094InfoDict, "No resistance-conferring mutations")
        elif name == "MEG_3979":
            intrinsicResistant(name, read.query_name, intrinsic)
        elif frameshiftCheckResult == "":
            addRead(name, read.query_name, resistantFrameshiftInfoDict, "No resistance-conferring mutations")
        elif gene.listOfMusts != None:
            intrinsicResistant(name, read.query_name, "Mutant")
        else:
            disregard(name, argInfoDict)
    return False