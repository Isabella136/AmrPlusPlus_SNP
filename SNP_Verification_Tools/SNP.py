from . import InDel
class SNP:
    def __init__(this, snpStringList, name):
        this.name = name
        snpString = snpStringList[0]
        this.snp = snpString
        this.wtOG = snpString[:1]
        snpString = snpString[1:]
        i = 0
        for x in snpString:
            if (not(x.isdigit())):
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
    def findACT(this, sequence, rRNA):
        this.wtACT = ""
        this.posACT = -1
        if (this.name == "MEG_4057") and (this.posOG == 347):
            this.wtACT = this.wtOG
            this.posACT = this.posOG
        else:
            begin = this.posOG - 35
            end = this.posOG + 26
            if end > (len(sequence) - 1 - len(this.rightContext) - len(this.leftContext)):
                end = len(sequence) - 1 - len(this.rightContext) - len(this.leftContext)
                begin = end - 61
            elif begin < 0:
                begin = 0
                end = 61
            for i in range(begin, end+1):
                if this.checkLeft(0, i, sequence, 0, rRNA):
                    this.changeACT(sequence, i)
                    break
    def checkLeft(this, index, currentPos, sequence, errorMargin, rRNA):
        if (len(this.leftContext) != 0):
            for aa in this.leftContext[index]:
                if aa == sequence[currentPos]:
                    if index == len(this.leftContext) - 1:
                        return this.checkRight(0, currentPos + 2, sequence, errorMargin, rRNA)
                    else:
                        return this.checkLeft(index + 1, currentPos + 1, sequence, errorMargin, rRNA)
        else:
            return this.checkRight(0, currentPos + 2, sequence, errorMargin, rRNA)
        if (errorMargin < 3) and not(rRNA):
            if index == len(this.leftContext) - 1:
                return this.checkRight(0, currentPos + 2, sequence, errorMargin + 1, rRNA)
            else:
                return this.checkLeft(index + 1, currentPos + 1, sequence, errorMargin + 1, rRNA)
        else:
            return False
    def checkRight(this, index, currentPos, sequence, errorMargin, rRNA):
        if (len(this.rightContext) != 0):
            for aa in this.rightContext[index]:
                if aa == sequence[currentPos]:
                    if index == len(this.rightContext) - 1:
                        return True
                    else:
                        return this.checkRight(index + 1, currentPos + 1, sequence, errorMargin, rRNA)
        else:
            return True
        if (errorMargin < 3) and not(rRNA):
            if index == len(this.rightContext) - 1:
                return True
            else:
                return this.checkRight(index + 1, currentPos + 1, sequence, errorMargin + 1, rRNA)
        else:
            return False
    def isValid(this):
        return this.posACT > -1

class SNP_Mis(SNP):
    def __init__(this, sequence, snpString, name, rRNA = False):
        snpStringList = [snpString]
        SNP.__init__(this, snpStringList, name)
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
        SNP.findACT(this, sequence, rRNA)
    def condensedInfo(this):
        wt = this.wtOG
        if (this.wtOG != this.wtACT):
            wt += this.wtACT
        return (wt, this.posACT, this.mtList, "Mis:" + this.snp)
    def changeACT(this, sequence, i):
        try: #if sequence contains mt
            this.mtList.index(sequence[i+len(this.leftContext)])
            this.wtACT = this.wtOG
            this.posACT = i+len(this.leftContext)+1
        except ValueError: 
            this.wtACT = sequence[i+len(this.leftContext)]
            this.posACT = i+len(this.leftContext)+1

class SNP_Non(SNP):
    def __init__(this, sequence, snpString, name):
        snpStringList = [snpString]
        SNP.__init__(this, snpStringList, name)
        snpString = snpStringList[0]
        snpString = snpString[2:]
        SNP.establishContext(this, snpString)
        SNP.findACT(this, sequence, False)
    def condensedInfo(this):
        wt = this.wtOG
        if (this.wtOG != this.wtACT):
            wt += this.wtACT
        return (wt, this.posACT, "*", "Nonsense:" + this.snp)
    def changeACT(this, sequence, i):
        this.wtACT = sequence[i+len(this.leftContext)]
        this.posACT = i+len(this.leftContext)+1

class Must(SNP):
    def __init__(this, wtString, name):
        tempList = [wtString]
        SNP.__init__(this, tempList, name)
        wtString = tempList[0]
        wtString = wtString[1:]
        SNP.establishContext(this, wtString)
        this.changeACT()
        this.next = None
    def condensedInfo(this):
        return (this.wtOG, this.posOG)
    def changeACT(this):
        this.wtACT = this.wtOG
        this.posACT = this.posOG
    def getPos(this):
        return this.posOG
    def defineNext(this,nextMust):
        this.next = nextMust
    def getNext(this):
        return this.next
    def getWt(this):
        return this.wtOG

class MustList:
    def __init__(this, wtString, name):
        this.listOfMust = {}
        this.aaOrNu = "amino acids"
        if wtString[:8] == "Must:Nuc":
            this.aaOrNu = "nucleic acids"
        infoList = wtString.split(";")
        for info in infoList:
            if info[:3] == "Nuc":
                temp = info[:info.find(';')][4:]
            else:
                temp = info[:info.find(';')][6:]
            wtToAdd = Must(temp, name)
            if len(this.listOfMust) > 0:
                this.listOfMust[list(this.listOfMust.keys())[-1]].defineNext(wtToAdd)
            else:
                this.firstPos = wtToAdd.getPos()
            this.listOfMust.update({wtToAdd.getPos():wtToAdd})
        this.lastPos = list(this.listOfMust.keys())[-1]
    def getFirstMustBetweenParams(this, begin, end):
        if (end < this.firstPos) or (begin > this.lastPos):
            return None
        if (end >= this.firstPos) and (begin <= this.firstPos): return (this.listOfMust[this.firstPos], True)
        for pos in this.listOfMust.keys():
            if (begin <= pos) and (end >= pos):
                return (this.listOfMust[pos], False)
        return None
    def returnAaOrNu(this):
        return this.aaOrNu

class SNP_Mult:
    def __init__(this, sequence, snpString, name):
        this.mult = snpString
        this.listOfMts = []
        snpsAndDels = snpString.split(";")
        for info in snpsAndDels:
            if (info[:3] == "Mis") or (info[:3] == "Nuc"):
                this.listOfMts.append(SNP_Mis(sequence, info[4:], name, True if info[:3]=='Nuc' else False))
            elif info[:3] == "Ins":
                this.listOfMts.append(InDel.Insertion(sequence, info[4:], name))
            else: #info[:3] = "Del"
                this.listOfMts.append(InDel.Deletion(sequence, info[4:], name))
    def condensedInfo(this):
        toReturn = []
        for snp in this.listOfMts:
            toReturn.append(snp.condensedInfo())
        return (toReturn, "Mult:" + this.mult)
    def isValid(this):
        for snp in this.listOfMts:
            if not(snp.isValid()):
                return False
        return True