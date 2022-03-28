from . import Gene
class SNP:
    def __init__(this, sequence, snpString):
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
        if sequence[this.posOG - 1] != this.wtOG:
            begin = 0
            end = len(sequence)
            if this.posOG - 1 > 30:
                begin = this.posOG - 31
            if this.posOG + 20 < len(sequence):
                end = this.posOG + 20
            i = begin
            temp = sequence[begin:end]
            for x in sequence[begin:end]:
                for laa0 in this.leftContext[0]:
                    if laa0 == x:
                        for laa1 in this.leftContext[1]:
                            if laa1 == sequence[i+1]:
                                for laa2 in this.leftContext[2]:
                                    if laa2 == sequence[i+2]:
                                        for laa3 in this.leftContext[3]:
                                            if laa3 == sequence[i+3]:
                                                for laa4 in this.leftContext[4]:
                                                    if laa4 == sequence[i+4]:
                                                        for raa0 in this.rightContext[0]:
                                                            if raa0 == sequence[i+6]:
                                                                for raa1 in this.rightContext[1]:
                                                                    if raa1 == sequence[i+7]:
                                                                        for raa2 in this.rightContext[2]:
                                                                            if raa2 == sequence[i+8]:
                                                                                for raa3 in this.rightContext[3]:
                                                                                    if raa3 == sequence[i+9]:
                                                                                        for raa4 in this.rightContext[4]:
                                                                                            if raa4 == sequence[i+10]:
                                                                                                this.wtACT = sequence[i+5]
                                                                                                this.posACT = i+5
                                                                                                break
                i += 1
            if this.posACT == -1:
                print("error")
        else:
            this.wtACT = this.wtOG
            this.posACT = this.posOG
