class InDel:
    def __init__(this, mtString, name, insertion):
        this.indel = mtString
        this.name = name
        i = 1
        this.inserted = None
        this.deleted = None
        if insertion:
            this.inserted = [mtString[0]]
            while not(mtString[i].isdigit()):
                this.inserted.append(mtString[i])
                i+=1
        else:
            this.deleted = mtString[0]
        posString = mtString[i:mtString.find("_")]
        mtString = mtString[mtString.find("_")+1:]
        this.posList = posString.split("/")
        this.posList = [int(i) for i in this.posList]
        InDel.establishContext(this, mtString)
    def establishContext(this, mtString) :
        contextList = mtString.split("_")
        goToNext = True
        this.leftContext = []
        for char in contextList[0]:
            if char == '[':
                goToNext = False
                this.leftContext.append([])
            elif char == ']':
                goToNext = True
            else:
                if goToNext:
                    this.leftContext.append([char])
                else:
                    this.leftContext[-1].append(char)
        this.rightContext = []
        for char in contextList[1]:
            if char == '[':
                goToNext = False
                this.rightContext.append([])
            elif char == ']':
                goToNext = True
            else:
                if goToNext:
                    this.rightContext.append([char])
                else:
                    this.rightContext[-1].append(char)
    def checkLeft(this, index, currentPos, sequence, errorMargin, rRNA = False):
        if (len(this.leftContext) != 0):
            for aa in this.leftContext[index]:
                if aa == sequence[currentPos]:
                    if index == len(this.leftContext) - 1:
                        nextPos = currentPos + (this.posList[-1]-this.posList[0]) + 1
                        nextPos = nextPos if this.inserted != None else nextPos+1
                        return InDel.checkRight(this, 0, nextPos , sequence, errorMargin, rRNA)
                    else:
                        return InDel.checkLeft(this, index + 1, currentPos + 1, sequence, errorMargin, rRNA)
        else:
            nextPos = currentPos + (this.posList[-1]-this.posList[0]) + 1
            nextPos = nextPos if this.inserted != None else nextPos+1
            return this.checkRight(0, nextPos, sequence, errorMargin, rRNA)
        if (errorMargin < 3) and not(rRNA):
            if index == len(this.leftContext) - 1:
                nextPos = currentPos + (this.posList[-1]-this.posList[0]) + 1
                nextPos = nextPos if this.inserted != None else nextPos+1
                return InDel.checkRight(this, 0, nextPos, sequence, errorMargin + 1, rRNA)
            else:
                return InDel.checkLeft(this, index + 1, currentPos + 1, sequence, errorMargin + 1, rRNA)
        else:
            return False
    def checkRight(this, index, currentPos, sequence, errorMargin, rRNA):
        if (len(this.rightContext) != 0):
            for aa in this.rightContext[index]:
                if aa == sequence[currentPos]:
                    if index == len(this.rightContext) - 1:
                        return True
                    else:
                        return InDel.checkRight(this, index + 1, currentPos + 1, sequence, errorMargin, rRNA)
        else:
            return True
        if (errorMargin < 3) and not(rRNA):
            if index == len(this.rightContext) - 1:
                return True
            else:
                return InDel.checkRight(this, index + 1, currentPos + 1, sequence, errorMargin + 1, rRNA)
        else:
            return False
    def changeACT(this, i):
        lastPos = 0
        diff = 0
        for pos in this.posList:
            if lastPos > 0:
                diff += (pos-lastPos)
            lastPos = pos
            this.posACT.append(i+len(this.leftContext)+1+diff)
    def isValid(this):
        return len(this.posACT) != 0

class Insertion(InDel):
    def __init__(this, sequence, mtString, name):
        InDel.__init__(this, mtString, name, True)
        this.findACT(sequence)
    def findACT(this, sequence):
        this.posACT = []
        begin = this.posList[0] - 35
        end = this.posList[-1] + 25
        if begin < 0:
            begin = 0
            end = begin + 60
        elif end > (len(sequence) - (this.posList[-1]-this.posList[0]) - len(this.leftContext) - len(this.rightContext)):
            end = len(sequence) - (this.posList[-1]-this.posList[0]) - len(this.leftContext) - len(this.rightContext)
            begin = end - 60
        for i in range(begin, end+1):
            if this.checkLeft(0, i, sequence, 0):
                this.changeACT(i)
                break
    def condensedInfo(this):
        return (this.inserted, this.posACT, "+", "Ins:" + this.indel)
    def longIndel(this):
        return len(this.inserted) >= 4

class Deletion(InDel):
    def __init__(this, sequence, mtString, name, rRNA = False):
        InDel.__init__(this, mtString, name, False)
        this.findACT(sequence, rRNA)
    def changeACT(this, sequence, i):
        this.deletionACT = sequence[i+len(this.leftContext)]
        InDel.changeACT(this, i)
    def findACT(this, sequence, rRNA):
        this.posACT = []
        begin = this.posList[0] - 35
        end = this.posList[-1] + 26
        if begin < 0:
            begin = 0
            end = begin + 61
        elif end > (len(sequence) - 1 - (this.posList[-1]-this.posList[0]) - len(this.leftContext) - len(this.rightContext)):
            end = len(sequence) - 1 - (this.posList[-1]-this.posList[0]) - len(this.leftContext) - len(this.rightContext)
            begin = end - 61
        for i in range(begin, end+1):
            if this.checkLeft(0, i, sequence, 0):
                this.changeACT(sequence, i)
                break
    def isValid(this):
        return (len(this.posACT) != 0) and (this.deleted == this.deletionACT)
    def condensedInfo(this):
        return (this.deleted, this.posACT, "-", "Del:" + this.indel)