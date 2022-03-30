from . import SNP
from . import dnaTranslate
class Gene:
    def __init__(this, name, sequence, snps):
        this.name = name
        this.sequence = sequence.upper()
        this.translated = dnaTranslate(this.sequence)
        this.listOfSNPs = []
        while(snps.find('|') != -1):
            temp = snps[:snps.find('|')]
            this.listOfSNPs.append(SNP.SNP(this.translated, temp, name))
            snps = snps[snps.find('|')+1:]
        this.listOfSNPs.append(SNP.SNP(this.translated, snps, name))
    def aaSequence(this):
        return this.translated

