from . import SNP
from . import dnaTranslate
class Gene:
    def __init__(this, name, sequence, snps):
        this.name = name + "|RequiresSNPConfirmation"
        this.sequence = sequence.upper()
        this.translated = dnaTranslate(this.sequence)
        this.listOfSNPs = []
        while(snps.find('|') != -1):
            temp = snps[:snps.find('|')]
            snpToAdd = SNP.SNP(this.translated, temp)
            if(snpToAdd.isSnpValid()):
                this.listOfSNPs.append(snpToAdd)
            snps = snps[snps.find('|')+1:]
        snpToAdd = SNP.SNP(this.translated, snps)
        if(snpToAdd.isSnpValid()):
            this.listOfSNPs.append(snpToAdd)
    def aaSequence(this):
        return this.translated
    def condensedInfo(this):
        condensedInfoList = []
        for snp in this.listOfSNPs :
            condensedInfoList.append(snp.condensedInfo())
        return condensedInfoList

