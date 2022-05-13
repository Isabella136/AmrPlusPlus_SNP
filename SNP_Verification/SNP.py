from multiprocessing.dummy import Value
from . import Gene
class SNP:
    def __init__(this, snpStringList):
        snpString = snpStringList[0]
        this.wtOG = snpString[:1]
        snpString = snpString[1:]
        i = 0
        for x in snpString:
            if not(x.isdigit()):
                break
            i += 1
        this.posOG = int(snpString[:i])
        snpString = snpString[i:]
        snpStringList[0] = snpString
    def establishContext(this, snpString):
        i = 1
        goToNext = True
        this.leftContext = []
        temp = []
        for x in snpString:
            if x == '_':
                snpString = snpString[i:]
                break
            elif x == '[':
                goToNext = False
                temp = []
            elif x == ']':
                this.leftContext.append(temp)
                goToNext = True
            else:
                if goToNext:
                    temp = [x]
                    this.leftContext.append(temp)
                else:
                    temp.append(x)
            i += 1
        this.rightContext = []
        for x in snpString:
            if x == '[':
                goToNext = False
                temp = []
            elif x == ']':
                this.rightContext.append(temp)
                goToNext = True
            else:
                if goToNext:
                    temp = [x]
                    this.rightContext.append(temp)
                else:
                    temp.append(x)
    def findACT(this, sequence):
        this.wtACT = ""
        this.posACT = -1
        if (len(sequence)-1 < this.posOG) or (sequence[this.posOG - 1] != this.wtOG):
            begin = this.posOG - 35
            end = this.posOG + 26
            if end > len(sequence)-1:
                end = len(sequence) - 1 - len(this.rightContext) - len(this.leftContext)
                begin = end - 61
            elif begin < 0:
                begin = 0
                end = 61
            i = begin
            temp = sequence[begin:end+10]
            for x in sequence[begin:end]:
                if this.checkLeft(0, i, sequence, 0):
                    this.changeACT(sequence, i)
                    break
                i += 1
        else:
            this.wtACT = this.wtOG
            this.posACT = this.posOG
    def checkLeft(this, index, currentPos, sequence, errorMargin):
        for aa in this.leftContext[index]:
            if aa == sequence[currentPos]:
                if index == len(this.leftContext) - 1:
                    return this.checkRight(0, currentPos + 2, sequence, errorMargin)
                else:
                    return this.checkLeft(index + 1, currentPos + 1, sequence, errorMargin)
        if errorMargin < 3:
            if index == len(this.leftContext) - 1:
                return this.checkRight(0, currentPos + 2, sequence, errorMargin + 1)
            else:
                return this.checkLeft(index + 1, currentPos + 1, sequence, errorMargin + 1)
        else:
            return False
    def checkRight(this, index, currentPos, sequence, errorMargin):
        for aa in this.rightContext[index]:
            if aa == sequence[currentPos]:
                if index == len(this.rightContext) - 1:
                    return True
                else:
                    return this.checkRight(index + 1, currentPos + 1, sequence, errorMargin)
        if errorMargin < 3:
            if index == len(this.rightContext) - 1:
                return True
            else:
                return this.checkRight(index + 1, currentPos + 1, sequence, errorMargin + 1)
        else:
            return False
    def isSnpValid(this):
        return this.posACT > -1
class SNP_Reg(SNP):
    def __init__(this, sequence, snpString):
        snpStringList = [snpString]
        SNP.__init__(this, snpStringList)
        snpString = snpStringList[0]
        this.mtList = []
        i = 1
        for x in snpString:
            if x == '_':
                snpString = snpString[i:]
                break
            this.mtList.append(x)
            i += 1
        SNP.establishContext(this, snpString)
        SNP.findACT(this, sequence)
    def condensedInfo(this):
        wt = this.wtOG
        if (this.wtOG != this.wtACT):
            wt += this.wtACT
        return (wt, this.posACT, this.mtList)
    def changeACT(this, sequence, i):
        try: #if sequence contains mt
            this.mtList.index(sequence[i+5])
            this.wtACT = this.wtOG
            this.posACT = i+5+1
        except ValueError: 
            this.wtACT = sequence[i+5]
            this.posACT = i+5+1
class SNP_Del(SNP):
    def __init__(this, sequence, snpString):
        snpStringList = [snpString]
        SNP.__init__(this, snpStringList)
        snpString = snpStringList[0]
        SNP.establishContext(this, snpString)
        SNP.findACT(this, sequence)
    def condensedInfo(this):
        wt = this.wtOG
        if (this.wtOG != this.wtACT):
            wt += this.wtACT
        return (wt, this.posACT, "-")
    def changeACT(this, sequence, i):
        this.wtACT = sequence[i+5]
        this.posACT = i+5+1
class SNP_Non(SNP):
    def __init__(this, sequence, snpString):
        snpStringList = [snpString]
        SNP.__init__(this, snpStringList)
        snpString = snpStringList[0]
        SNP.establishContext(this, snpString)
        SNP.findACT(this, sequence)
    def condensedInfo(this):
        wt = this.wtOG
        if (this.wtOG != this.wtACT):
            wt += this.wtACT
        return (wt, this.posACT, "*")
    def changeACT(this, sequence, i):
        this.wtACT = sequence[i+5]
        this.posACT = i+5+1
class SNP_Mult:
    def __init__(this, sequence, snpString):
        this.listOfSNPs = []
        while(snpString.find(';') != -1):
            temp = snpString[:snpString.find(';')]
            if temp[:3] == "Reg":
                snpToAdd = SNP_Reg(sequence, temp[4:])
                if(snpToAdd.isSnpValid()):
                    this.listOfSNPs.append(snpToAdd)
            else: # temp[:3] == "Del"
                snpToAdd = SNP_Del(sequence, temp[4:])
                if(snpToAdd.isSnpValid()):
                    this.listOfSNPs.append(snpToAdd)
            snpString = snpString[snpString.find(';')+1:]
        if snpString[:3] == "Reg":
            snpToAdd = SNP_Reg(sequence, snpString[4:])
            if(snpToAdd.isSnpValid()):
                this.listOfSNPs.append(snpToAdd)
        else: # snpString[:3] == "Del"
            snpToAdd = SNP_Del(sequence, snpString[4:])
            if(snpToAdd.isSnpValid()):
                this.listOfSNPs.append(snpToAdd)
    def condensedInfo(this):
        toReturn = []
        for snp in this.listOfSNPs:
            toReturn.append(snp.condensedInfo)
        return toReturn
    def isSnpValid(this):
        for snp in this.listOfSNPs:
            if not(snp.isSnpValid()):
                return False
        return True