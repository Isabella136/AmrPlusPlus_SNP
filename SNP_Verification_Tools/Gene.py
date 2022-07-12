from typing import Tuple
from . import SNP
from . import dnaTranslate
class Gene:
    def __init__(this, name, sequence, snps):
        this.name = name[:name.find('|')]
        this.fullName = name + "|RequiresSNPConfirmation"
        this.sequence = sequence.upper()
        if "Nuc" in snps:
            this.translated = None
        else:
            this.translated = dnaTranslate(this.sequence)
        this.listOfMisSNPs = []
        this.listOfDelSNPs = []
        this.listOfNonSNPs = []
        this.listOfMultSNPs = []
        this.listsOfMusts = []
        while(snps.find('|') != -1):
            temp = snps[:snps.find('|')]
            if temp[:4] == "Mult":
                snpToAdd = SNP.SNP_Mult(this.translated, temp[5:], this.name)
                if(snpToAdd.isSnpValid()):
                    this.listOfMultSNPs.append(snpToAdd)
            elif temp[:3] == "Mis":
                snpToAdd = SNP.SNP_Mis(this.translated, temp[4:], this.name)
                if(snpToAdd.isSnpValid()):
                    this.listOfMisSNPs.append(snpToAdd)
            elif temp[:3] == "Del":
                snpToAdd = SNP.SNP_Del(this.translated, temp[4:], this.name)
                if(snpToAdd.isSnpValid()):
                    this.listOfDelSNPs.append(snpToAdd)
            elif temp[:7] == "NucMult":
                snpToAdd = SNP.SNP_Mult(this.sequence, temp[8:], this.name)
                if(snpToAdd.isSnpValid()):
                    this.listOfMultSNPs.append(snpToAdd)
            elif temp[:6] == "NucDel":
                snpToAdd = SNP.SNP_Del(this.sequence, temp[7:], this.name)
                if(snpToAdd.isSnpValid()):
                    this.listOfDelSNPs.append(snpToAdd)
            elif temp[:3] == "Nuc":
                snpToAdd = SNP.SNP_Mis(this.sequence, temp[4:], this.name)
                if(snpToAdd.isSnpValid()):
                    this.listOfMisSNPs.append(snpToAdd)
            elif temp[:4] == "Must":
                wtToAdd = SNP.MustList(temp[5:], this.name)
                this.listsOfMusts.append(wtToAdd)
            else: #temp[:3] == "Non" 
                snpToAdd = SNP.SNP_Non(this.translated, temp[4:], this.name)
                if(snpToAdd.isSnpValid()):
                    this.listOfNonSNPs.append(snpToAdd)
            snps = snps[snps.find('|')+1:]
        if snps[:4] == "Mult":
            snpToAdd = SNP.SNP_Mult(this.translated, snps[5:], this.name)
            if(snpToAdd.isSnpValid()):
                this.listOfMultSNPs.append(snpToAdd)
        elif snps[:3] == "Mis":
            snpToAdd = SNP.SNP_Mis(this.translated, snps[4:], this.name)
            if(snpToAdd.isSnpValid()):
                this.listOfMisSNPs.append(snpToAdd)
        elif snps[:3] == "Del":
            snpToAdd = SNP.SNP_Del(this.translated, snps[4:], this.name)
            if(snpToAdd.isSnpValid()):
                this.listOfDelSNPs.append(snpToAdd)
        elif snps[:7] == "NucMult":
            snpToAdd = SNP.SNP_Mult(this.sequence, snps[8:], this.name)
            if(snpToAdd.isSnpValid()):
                this.listOfMultSNPs.append(snpToAdd)
        elif snps[:6] == "NucDel":
            snpToAdd = SNP.SNP_Del(this.sequence, snps[7:], this.name)
            if(snpToAdd.isSnpValid()):
                this.listOfDelSNPs.append(snpToAdd)
        elif snps[:3] == "Nuc":
            snpToAdd = SNP.SNP_Mis(this.sequence, snps[4:], this.name)
            if(snpToAdd.isSnpValid()):
                this.listOfMisSNPs.append(snpToAdd)
        elif snps[:4] == "Must":
            wtToAdd = SNP.MustList(snps[5:], this.name)
            this.listsOfMusts.append(wtToAdd)
        else: #snps[:3] == "Non" 
            snpToAdd = SNP.SNP_Non(this.translated, snps[4:], this.name)
            if(snpToAdd.isSnpValid()):
                this.listOfNonSNPs.append(snpToAdd)
    def aaSequence(this):
        return this.translated
    def ntSequence(this):
        return this.sequence
    def condensedMultInfo(this):
        condensedInfoList = []
        for snp in this.listOfMultSNPs :
            condensedInfoList.append(snp.condensedInfo())
        return condensedInfoList
    def condensedRegDelInfo(this):
        condensedInfoList = []
        for snp in this.listOfMisSNPs :
            condensedInfoList.append(snp.condensedInfo())
        for snp in this.listOfDelSNPs :
            condensedInfoList.append(snp.condensedInfo())
        return condensedInfoList
    def condensedNonInfo(this):
        condensedInfoList = []
        for snp in this.listOfNonSNPs :
            condensedInfoList.append(snp.condensedInfo())
        return condensedInfoList
    def getFirstMustBetweenParams(this, begin, end):
        if len(this.listsOfMusts) == 0:
            return None
        toReturn = []
        for list in this.listsOfMusts:
            returnedMust = list.getFirstMustBetweenParams(begin, end)
            if returnedMust != None:
                toReturn.append(returnedMust)
        if len(toReturn) == 0:
            return None
        return toReturn
    def getName(this):
        return this.name
    def getFullName(this):
        return this.fullName
    def ntSeqLength(this):
        return len(this.sequence)
    def aaOrNu(this):
        return this.listsOfMusts[0].returnAaOrNu()
    def rRna(this):
        return (("16S" in this.fullName) or ("23S" in this.fullName)) and (this.name != "MEG_6144") 