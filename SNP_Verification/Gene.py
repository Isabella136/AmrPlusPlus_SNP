from . import SNP
from . import dnaTranslate
class Gene:
    def __init__(this, name, sequence, snps):
        this.name = name[:name.find('|')]
        this.sequence = sequence.upper()
        this.translated = dnaTranslate(this.sequence)
        this.listOfRegSNPs = []
        this.listOfDelSNPs = []
        this.listOfNonSNPs = []
        this.listOfMultSNPs = []
        while(snps.find('|') != -1):
            temp = snps[:snps.find('|')]
            if temp[:4] == "Mult":
                snpToAdd = SNP.SNP_Mult(this.translated, temp[5:], this.name)
                if(snpToAdd.isSnpValid()):
                    this.listOfMultSNPs.append(snpToAdd)
                else:
                    print(name + ": " + temp)
            elif temp[:3] == "Reg":
                snpToAdd = SNP.SNP_Reg(this.translated, temp[4:], this.name)
                if(snpToAdd.isSnpValid()):
                    this.listOfRegSNPs.append(snpToAdd)
                else:
                    print(name + ": " + temp)
            elif temp[:3] == "Del":
                snpToAdd = SNP.SNP_Del(this.translated, temp[4:], this.name)
                if(snpToAdd.isSnpValid()):
                    this.listOfDelSNPs.append(snpToAdd)
                else:
                    print(name + ": " + temp)
            else: #temp[:3] == "Non" 
                snpToAdd = SNP.SNP_Non(this.translated, temp[4:], this.name)
                if(snpToAdd.isSnpValid()):
                    this.listOfNonSNPs.append(snpToAdd)
                else:
                    print(name + ": " + temp)
            snps = snps[snps.find('|')+1:]
        if snps[:4] == "Mult":
            snpToAdd = SNP.SNP_Mult(this.translated, snps[5:], this.name)
            if(snpToAdd.isSnpValid()):
                this.listOfMultSNPs.append(snpToAdd)
            else:
                print(name + ": " + snps)
        elif snps[:3] == "Reg":
            snpToAdd = SNP.SNP_Reg(this.translated, snps[4:], this.name)
            if(snpToAdd.isSnpValid()):
                this.listOfRegSNPs.append(snpToAdd)
            else:
                print(name + ": " + snps)
        elif snps[:3] == "Del":
            snpToAdd = SNP.SNP_Del(this.translated, snps[4:], this.name)
            if(snpToAdd.isSnpValid()):
                this.listOfDelSNPs.append(snpToAdd)
            else:
                print(name + ": " + snps)
        else: #snps[:3] == "Non" 
            snpToAdd = SNP.SNP_Non(this.translated, snps[4:], this.name)
            if(snpToAdd.isSnpValid()):
                this.listOfNonSNPs.append(snpToAdd)
            else:
                print(name + ": " + snps)
    def aaSequence(this):
        return this.translated
    def condensedMultInfo(this):
        condensedInfoList = []
        for snp in this.listOfMultSNPs :
            condensedInfoList.append(snp.condensedInfo())
        return condensedInfoList
    def condensedRegDelInfo(this):
        condensedInfoList = []
        for snp in this.listOfRegSNPs :
            condensedInfoList.append(snp.condensedInfo())
        for snp in this.listOfDelSNPs :
            condensedInfoList.append(snp.condensedInfo())
        return condensedInfoList
    def condensedNonInfo(this):
        condensedInfoList = []
        for snp in this.listOfNonSNPs :
            condensedInfoList.append(snp.condensedInfo())
        return condensedInfoList
    def getName(this):
        return this.name
    def ntSeqLength(this):
        return len(this.sequence)
