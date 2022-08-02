from . import SNP
from . import InDel
from . import dnaTranslate
class Gene:
    def __init__(this, name, sequence, infoString):
        this.name = name[:name.find('|')]
        this.fullName = name + "|RequiresSNPConfirmation"
        this.sequence = sequence.upper()
        if "Nuc" in infoString:
            this.translated = None
        else:
            this.translated = dnaTranslate(this.sequence, name)

        this.listOfHyperSNPs = []
        this.listOfMisSNPs = []
        this.listOfNonsenseSNPs = []
        this.listOfMultSNPs = []
        this.frameshiftInfo = None
        this.isThereANonstop = (False, None)
        this.listOfDel = []
        this.listOfIns = []
        this.listOfMusts = None

        infoList = infoString.split('|')
        for info in infoList:
            if "Nonstop" in info:
                if info[:4] == "Mult":
                    snpToAdd = SNP.SNP_Mis(this.translated, info[9:info.find(";")], this.name)
                    if (snpToAdd.isSnpValid()):
                        this.isThereANonstop = (True, snpToAdd)
                else:
                    this.isThereANonstop = (True, None)
            elif info[:4] == "Mult":
                snpToAdd = SNP.SNP_Mult(this.translated, info[5:], this.name)
                if(snpToAdd.isSnpValid()):
                    this.listOfMultSNPs.append(snpToAdd)
            elif info[:3] == "Mis":
                snpToAdd = SNP.SNP_Mis(this.translated, info[4:], this.name)
                if(snpToAdd.isSnpValid()):
                    this.listOfMisSNPs.append(snpToAdd)
            elif info[:3] == "Del":
                indelToAdd = InDel.Deletion(this.translated, info[4:], this.name)
                if(indelToAdd.isSnpValid()):
                    this.listOfDel.append(indelToAdd)
            elif info[:3] == "Ins":
                indelToAdd = InDel.Deletion(this.translated, info[4:], this.name)
                if(indelToAdd.isSnpValid()):
                    this.listOfIns.append(indelToAdd)
            elif info[:5] == "Hyper":
                snpToAdd = SNP.SNP_Mult(this.translated, info[6:], this.name)
                if(snpToAdd.isSnpValid()):
                    this.listOfHyperSNPs.append(snpToAdd)
            elif info[:3] == "FS-":
                this.frameshiftInfo = info[3:]
            elif info[:7] == "NucMult":
                snpToAdd = SNP.SNP_Mult(this.sequence, info[8:], this.name)
                if(snpToAdd.isSnpValid()):
                    this.listOfMultSNPs.append(snpToAdd)
            elif info[:6] == "NucDel":
                indelToAdd = InDel.Deletion(this.sequence, info[7:], this.name)
                if(indelToAdd.isSnpValid()):
                    this.listOfDel.append(indelToAdd)
            elif info[:3] == "Nuc":
                snpToAdd = SNP.SNP_Mis(this.sequence, info[4:], this.name)
                if(snpToAdd.isSnpValid()):
                    this.listOfMisSNPs.append(snpToAdd)
            elif info[:4] == "Must":
                this.listOfMusts = SNP.MustList(info[5:], this.name)
            else: #temp[:3] == "Nonsense" 
                snpToAdd = SNP.SNP_Non(this.translated, info[4:], this.name)
                if(snpToAdd.isSnpValid()):
                    this.listOfNonsenseSNPs.append(snpToAdd)

    def aaSequence(this):
        return this.translated
    def ntSequence(this):
        return this.sequence

    def condensedMultInfo(this):
        condensedInfoList = []
        for snp in this.listOfMultSNPs :
            condensedInfoList.append(snp.condensedInfo())
        return condensedInfoList
    def condensedHyperInfo(this):
        return this.listOfHyperSNPs[0].condensedInfo()
    def condensedMisInDelInfo(this):
        condensedInfoList = []
        for snp in this.listOfMisSNPs :
            condensedInfoList.append(snp.condensedInfo())
        for snp in this.listOfIns:
            condensedInfoList.append(snp.condensedInfo())
        for snp in this.listOfDel :
            condensedInfoList.append(snp.condensedInfo())
        return condensedInfoList
    def condensedNonInfo(this):
        condensedInfoList = []
        for snp in this.listOfNonsenseSNPs :
            condensedInfoList.append(snp.condensedInfo())
        return condensedInfoList
    def getNonstopInfo(this):
        return this.isThereANonstop
    def getFirstMustBetweenParams(this, begin, end):
        if this.listOfMusts == None:
            return None
        return this.listOfMusts.getFirstMustBetweenParams(begin, end)

    def getName(this):
        return this.name
    def getFullName(this):
        return this.fullName
    def ntSeqLength(this):
        return len(this.sequence)
    def aaOrNu(this):
        return this.listOfMusts.returnAaOrNu()
    def rRna(this):
        return (("16S" in this.fullName) or ("23S" in this.fullName)) and (this.name != "MEG_6144") 