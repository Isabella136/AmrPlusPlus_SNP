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
        FinalCount(gene, read)
        return None                                
    
    seqOfInterest, mapOfInterest = MapQueryToReference(rRna, read, gene)        #If S-tagged and needs suppression, seq and map of interest, removes index 1602
    if seqOfInterest == False:
        FinalCount(gene, read)
        return None

    if not(rRna):                                                               #rRNA stays as nucleotide sequence; nonsense mutations don't matter   
        checkResult = NonsenseCheck(read, gene, mapOfInterest, seqOfInterest)   #Checks for nonsense previously found in literature and for new nonsense
        if not(checkResult):                                                    #If not F-tagged and has new nonsense, can't determine resistance
            FinalCount(gene, read)
            return None      
    
    IntrinsicCheck(read, gene, mapOfInterest, seqOfInterest)                    #For I-tagged specifically
    MisInDelCheck(read, gene, mapOfInterest, seqOfInterest)                     #Counts all resistance-conferring mutations
    nTupleCheck(read, gene, mapOfInterest, seqOfInterest)                       #Counts all resistance-conferring mutations
    FinalCount(gene, read)    


def FinalCount(gene, read):
    gene.redefineLastTupleInfo(read)
    additionalInfo = gene.getLastTupleInfo()
    def nTypeCount(): 
        insertion = False
        deletion  = False
        missense  = False
        nTuple    = False
        resistant = False
        eight     = False
        nine      = False
        for info in additionalInfo[1:]:
            if None == info:
                break
            elif 'nonstop' == info:
                gene.addToOutputInfo(7)
                resistant = True
            elif 'Mult:' in info:
                if not(nTuple):
                    gene.addToOutputInfo(6)
                    nTuple = True
                resistant = True
            elif 'Mis:' in info:
                if not(missense):
                    gene.addToOutputInfo(2)
                    missense = True
                resistant = True
            elif 'Ins:' in info:
                if not(insertion):
                    gene.addToOutputInfo(3)
                    insertion = True
                resistant = True
            elif 'Del:' in info:
                if not(deletion):
                    gene.addToOutputInfo(4)
                    deletion = True
                resistant = True
            elif 'Nonsense:' in info:
                gene.addToOutputInfo(5)
                resistant = True
            elif '12+bp indel:' in info:
                count = int(info.split(":")[1][1:]) - gene.longIndel
                eight = count > 0
            elif '12+bp frameshift:' in info:
                nine = True
            elif "Newly found nonsense" in info:
                gene.addToOutputInfo(10)
            elif 'FS till end' == info:
                gene.addToOutputInfo(11)
        if resistant:
            gene.addDetails(read.query_name, "res")
            gene.addToOutputInfo(1)
            if eight:
                gene.addToOutputInfo(8)
            if nine:
                gene.addToOutputInfo(9)
        else:
            gene.addDetails(read.query_name, "sus")

    def fTypeCount():
        insertion = False
        deletion  = False
        missense  = False
        nTuple    = False
        resistant = False
        eight     = False
        nine      = False
        for info in additionalInfo[1:]:
            if None == info:
                break
            elif 'nonstop' == info:
                gene.addToOutputInfo(7)
                resistant = True
            elif 'Mult:' in info:
                if not(nTuple):
                    gene.addToOutputInfo(6)
                    nTuple = True
                resistant = True
            elif 'Mis:' in info:
                if not(missense):
                    gene.addToOutputInfo(2)
                    missense = True
                resistant = True
            elif 'Ins:' in info:
                if not(insertion):
                    gene.addToOutputInfo(3)
                    insertion = True
                resistant = True
            elif 'Del:' in info:
                if not(deletion):
                    gene.addToOutputInfo(4)
                    deletion = True
                resistant = True
            elif 'Nonsense:' in info:
                gene.addToOutputInfo(5)
                resistant = True
            elif '12+indel:' in info:
                count = int(info.split(":")[1][1:]) - gene.longIndel
                eight = count > 0
            elif '12+fs:' in info:
                nine = True
            elif "Newly found nonsense" in info:
                gene.addToOutputInfo(10)
                resistant = True
            elif 'FS till end' == info:
                gene.addToOutputInfo(11)
                resistant = True
        if resistant:
            gene.addToOutputInfo(1)
            gene.addDetails(read.query_name, "res")
        else:
            gene.addDetails(read.query_name, "sus")
            if eight:
                gene.addToOutputInfo(8)
            if nine:
                gene.addToOutputInfo(9)

    def hTypeCount(): 
        insertion = False
        deletion  = False
        missense  = False
        nTuple    = False
        nonsense  = False
        nonstop   = False
        resistant = False
        hyper     = False
        eight     = False
        nine      = False
        for info in additionalInfo[1:]:
            if None == info:
                break
            elif 'nonstop' == info:
                nonstop = True
                resistant = True
            elif 'Mult:' in info:
                nTuple = True
                resistant = True
            elif 'Hypersusceptible' in info:
                hyper = True
            elif 'Mis:' in info:
                missense = True
                resistant = True
            elif 'Ins:' in info:
                insertion = True
                resistant = True
            elif 'Del:' in info:
                deletion = True
                resistant = True
            elif 'Nonsense:' in info:
                nonsense = True
                resistant = True
            elif '12+indel:' in info:
                count = int(info.split(":")[1][1:]) - gene.longIndel
                eight = count > 0
            elif '12+fs:' in info:
                nine = True
            elif "Newly found nonsense" in info:
                gene.addToOutputInfo(10)
            elif 'FS till end' == info:
                gene.addToOutputInfo(11) 
        if resistant :
            gene.addDetails(read.query_name, "res")
            if not(hyper):
                if missense:
                    gene.addToOutputInfo(2)
                elif insertion:
                    gene.addToOutputInfo(3)
                elif deletion:
                    gene.addToOutputInfo(4)
                elif nonsense:
                    gene.addToOutputInfo(5)
                elif nTuple:
                    gene.addToOutputInfo(6)
                elif nonstop:
                    gene.addToOutputInfo(7)
                if eight:
                    gene.addToOutputInfo(8)
                if nine:
                    gene.addToOutputInfo(9)
                gene.addToOutputInfo(1)
            else:
                gene.addToOutputInfo(12)
        else:
            gene.addDetails(read.query_name, "sus")

    def sTypeCount():
        insertion = False
        deletion  = False
        missense  = False
        nTuple    = False
        resistant = False
        eight     = False
        nine      = False
        supr      = False
        followed  = False
        for info in additionalInfo[1:]:
            if None == info:
                break
            elif 'nonstop' == info:
                gene.addToOutputInfo(7)
                resistant = True
            elif 'Mult:' in info:
                if not(nTuple):
                    gene.addToOutputInfo(6)
                    nTuple = True
                resistant = True
            elif 'Mis:' in info:
                if not(missense):
                    gene.addToOutputInfo(2)
                    missense = True
                resistant = True
            elif 'Ins:' in info:
                if not(insertion):
                    gene.addToOutputInfo(3)
                    insertion = True
                resistant = True
            elif 'Del:' in info:
                if not(deletion):
                    gene.addToOutputInfo(4)
                    deletion = True
                resistant = True
            elif 'Nonsense:' in info:
                gene.addToOutputInfo(5)
                resistant = True
            elif '12+indel:' in info:
                count = int(info.split(":")[1][1:]) - gene.longIndel
                eight = count > 0
            elif '12+fs:' in info:
                nine = True
            elif "Newly found nonsense" in info:
                gene.addToOutputInfo(10)
                resistant = False
            elif 'FS till end' == info:
                gene.addToOutputInfo(11)
                resistant = False
            elif 'Suppressible C insert' == info:
                supr = True
                resistant = True
            elif 'C insert followed by del/ins' == info:
                followed = True
                resistant = True
        if resistant:
            gene.addDetails(read.query_name, "res")
            gene.addToOutputInfo(1)
            if eight:
                gene.addToOutputInfo(8)
            if nine:
                gene.addToOutputInfo(9)
            if supr:
                gene.addToOutputInfo(12)
            if followed:
                gene.addToOutputInfo(13)
        else:
            gene.addDetails(read.query_name, "sus")

    def iTypeCount(): 
        acquired = False
        some     = False
        none     = False
        mutant   = False
        six      = False
        seven    = False
        _all     = False
        for info in additionalInfo[1:]:
            if None == info:
                break
            elif 'All' == info:
                gene.addToOutputInfo(1)
                _all = True
                gene.addDetails(read.query_name, "res")
                if six:
                    gene.addToOutputInfo(6)
                if seven:
                    gene.addToOutputInfo(7)
                break
            elif 'Some' == info:
                some = True
            elif 'NA' == info:
                none = True
            elif 'Mutant' == info:
                mutant = True
            elif 'Mis:' in info:
                gene.addToOutputInfo(5)
                gene.addDetails(read.query_name, "res")
                acquired = True
                if six:
                    gene.addToOutputInfo(6)
                if seven:
                    gene.addToOutputInfo(7)
                break
            elif '12+indel:' in info:
                six = True
            elif '12+fs:' in info:
                seven = True
            elif "Newly found nonsense" in info:
                gene.addToOutputInfo(8)
            elif 'FS till end' == info:
                gene.addToOutputInfo(9)
        if not(acquired):
            if some:
                gene.addDetails(read.query_name, "res")
                gene.addToOutputInfo(2)
                if six:
                    gene.addToOutputInfo(6)
                if seven:
                    gene.addToOutputInfo(7)
            elif none:
                gene.addDetails(read.query_name, "sus")
                gene.addToOutputInfo(3)
            elif mutant:
                gene.addDetails(read.query_name, "sus")
                gene.addToOutputInfo(4)
            elif not _all:
                gene.addDetails(read.query_name, "sus")

    tag = gene.getGeneTag()
    if tag == 'N':
        nTypeCount()
    elif tag == 'F':
        fTypeCount()
    elif tag == 'H':
        hTypeCount()
    elif tag == 'S':
        sTypeCount()
    else:
        iTypeCount()