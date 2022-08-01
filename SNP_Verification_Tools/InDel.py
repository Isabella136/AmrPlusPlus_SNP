class InDel:
    def __init__(this, mtString, name, insertion):
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
        for char in contextList[0]:
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
    def checkLeft(this, index, currentPos, sequence, errorMargin):
        if (len(this.leftContext) != 0):
            for aa in this.leftContext[index]:
                if aa == sequence[currentPos]:
                    if index == len(this.leftContext) - 1:
                        return InDel.checkRight(this, 0, currentPos + len(this.posList), sequence, errorMargin)
                    else:
                        return InDel.checkLeft(this, index + 1, currentPos + 1, sequence, errorMargin)
        else:
            return this.checkRight(0, currentPos + len(this.posList), sequence, errorMargin)
        if errorMargin < 3:
            if index == len(this.leftContext) - 1:
                return InDel.checkRight(this, 0, currentPos + len(this.posList), sequence, errorMargin + 1)
            else:
                return InDel.checkLeft(this, index + 1, currentPos + 1, sequence, errorMargin + 1)
        else:
            return False
    def checkRight(this, index, currentPos, sequence, errorMargin):
        for aa in this.rightContext[index]:
            if aa == sequence[currentPos]:
                if index == len(this.rightContext) - 1:
                    return True
                else:
                    return InDel.checkRight(this, index + 1, currentPos + 1, sequence, errorMargin)
        if errorMargin < 3:
            if index == len(this.rightContext) - 1:
                return True
            else:
                return InDel.checkRight(this, index + 1, currentPos + 1, sequence, errorMargin + 1)
        else:
            return False
    def changeACT(this, i):
        lastPos = 0
        diff = 0
        for pos in this.posList:
            if lastPos > 0:
                diff += (pos-lastPos)
            lastPos = pos
            this.posACT.append(i+5+1+diff)
    def isValid(this):
        return this.posACT > -1

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
        elif end > (len(sequence) - len(this.posList) - len(this.leftContext) - len(this.rightContext)):
            end = len(sequence) - len(this.posList) - len(this.leftContext) - len(this.rightContext)
            begin = end - 60
        i = begin
        for x in sequence[begin:end]:
            if InDel.checkLeft(this, 0, i, sequence, 0):
                InDel.changeACT(this,i)
                break
            i+=1
    def condensedInfo(this):
        return (this.inserted, this.posACT, "+")

class Deletion(InDel):
    def __init__(this, sequence, mtString, name):
        InDel.__init__(this, mtString, name, False)
        this.findACT(sequence)
    def changeACT(this, sequence, i):
        this.deletionACT = sequence[i+5]
        InDel.changeACT(this, i)
    def findACT(this, sequence):
        this.posACT = []
        begin = this.posList[0] - 35
        end = this.posList[-1] + 26
        if begin < 0:
            begin = 0
            end = begin + 61
        elif end > (len(sequence) - 1 - len(this.posList) - len(this.leftContext) - len(this.rightContext)):
            end = len(sequence) - 1 - len(this.posList) - len(this.leftContext) - len(this.rightContext)
            begin = end - 61
        i = begin
        for x in sequence[begin:end]:
            if InDel.checkLeft(this, 0, i, sequence, 0):
                Deletion.changeACT(this, sequence, i)
                break
            i+=1
    def isValid(this):
        return (this.posACT > -1) and (this.deleted == this.deletionACT)
    def condensedInfo(this):
        return (this.deleted, this.posACT, "-")