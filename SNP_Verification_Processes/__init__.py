from SNP_Verification_Tools.Gene import Gene
from SNP_Verification_Processes.FrameshiftCheck import FrameshiftCheck
from SNP_Verification_Processes.MappingQueryToReference import MapQueryToReference
from SNP_Verification_Processes.NonsenseCheck import NonsenseCheck
from SNP_Verification_Processes.IntrinsicCheck import IntrinsicCheck
from SNP_Verification_Processes.MisInDelCheck import MisInDelCheck
from SNP_Verification_Processes.nTupleCheck import nTupleCheck
import pysam

def verify(read, gene):
    rRna = gene.rRna()
    gene.addToOutputInfo(0)                                                     #Counts read; other data will be counted in FinalCount method

                                                                                #Checks for frameshifts (if not rRNA), but also for extended indels;
    checkResult = FrameshiftCheck(read, gene, rRna)                             #for S-tagged, also determines if needs suppression
    if not(checkResult):                                                        #If not F-tagged and has frameshifts till the end of query sequence
        FinalCount(gene)                                
    
    seqOfInterest, mapOfInterest = MapQueryToReference(rRna, read, gene)        #If S-tagged and needs suppression, seq and map of interest, removes index 1602
    
    if not(rRna):                                                               #rRNA stays as nucleotide sequence; nonsense mutations don't matter   
        checkResult = NonsenseCheck(read, gene, mapOfInterest, seqOfInterest)   #Checks for nonsense previously found in literature and for new nonsense
        if not(checkResult):                                                    #If not F-tagged and has new nonsense, can't determine resistance
            FinalCount(gene)
    
    IntrinsicCheck(read, gene, mapOfInterest, seqOfInterest)                    #For I-tagged specifically
    MisInDelCheck(read, gene, mapOfInterest, seqOfInterest)                     #Counts all resistance-conferring mutations
    nTupleCheck(read, gene, mapOfInterest, seqOfInterest)                       #Counts all resistance-conferring mutations
    


def FinalCount(gene):
    def nTypeCount()
    def fTypeCount()
    def hTypeCount()
    def sTypeCount()
    def iTypeCount()