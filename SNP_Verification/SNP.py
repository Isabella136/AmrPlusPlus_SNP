from . import Gene
class SNP:
    def __init__(this, sequence, snpString, name):
        this.wtOG = snpString[:1]
        snpString = snpString[1:]
        i = 0
        for x in snpString:
            if x.isalpha():
                break
            i += 1
        this.posOG = int(snpString[:i])
        snpString = snpString[i:]
        this.mtList = []
        i = 1
        for x in snpString:
            if x == '_':
                snpString = snpString[i:]
                break
            this.mtList.append(x)
            i += 1
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
        this.wtACT = ""
        this.posACT = -1
        if (len(sequence) < this.posOG) or (sequence[this.posOG - 1] != this.wtOG):
            i = 0
            for x in sequence[:0 - 1 - len(this.rightContext) - len(this.leftContext)]:
                if this.checkLeft(0, i, sequence, 0):
                    this.wtACT = sequence[i+5]
                    this.posACT = i+5
                    break
                i += 1
            if this.posACT == -1:
                print(name)
                print(this.posOG)
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