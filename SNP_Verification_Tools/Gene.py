import string
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
        this.geneTag = 'N'
        this.outputInfo = dict()
        this.additionalInfo = list()
        this.longIndel = 0
        this.currentReadNonstop = None

        infoList = infoString.split('|')
        for info in infoList:
            if "Nonstop" in info:
                if info[:4] == "Mult":
                    snpToAdd = SNP.SNP_Mis(this.translated, info[9:info.find(";")], this.name)
                    if (snpToAdd.isValid()):
                        this.isThereANonstop = (True, snpToAdd)
                else:
                    this.isThereANonstop = (True, None)
            elif info[:4] == "Mult":
                snpToAdd = SNP.SNP_Mult(this.translated, info[5:], this.name)
                if(snpToAdd.isValid()):
                    this.listOfMultSNPs.append(snpToAdd)
            elif info[:3] == "Mis":
                snpToAdd = SNP.SNP_Mis(this.translated, info[4:], this.name)
                if(snpToAdd.isValid()):
                    this.listOfMisSNPs.append(snpToAdd)
            elif info[:3] == "Del":
                indelToAdd = InDel.Deletion(this.translated, info[4:], this.name)
                if(indelToAdd.isValid()):
                    this.listOfDel.append(indelToAdd)
            elif info[:3] == "Ins":
                indelToAdd = InDel.Insertion(this.translated, info[4:], this.name)
                if(indelToAdd.isValid()):
                    this.listOfIns.append(indelToAdd)
            elif info[:5] == "Hyper":
                snpToAdd = SNP.SNP_Mult(this.translated, info[6:], this.name)
                if(snpToAdd.isValid()):
                    this.listOfHyperSNPs.append(snpToAdd)
                this.geneTag = 'H'
            elif info[:3] == "FS-":
                this.frameshiftInfo = info[3:]
                if (this.name != "MEG_6142") and ("miscellaneous" not in info):
                    this.geneTag = 'S' if this.name == "MEG_6094" else 'F'
            elif info[:7] == "NucMult":
                snpToAdd = SNP.SNP_Mult(this.sequence, info[8:], this.name)
                if(snpToAdd.isValid()):
                    this.listOfMultSNPs.append(snpToAdd)
            elif info[:6] == "NucDel":
                indelToAdd = InDel.Deletion(this.sequence, info[7:], this.name, True)
                if(indelToAdd.isValid()):
                    this.listOfDel.append(indelToAdd)
            elif info[:3] == "Nuc":
                snpToAdd = SNP.SNP_Mis(this.sequence, info[4:], this.name, True)
                if(snpToAdd.isValid()):
                    this.listOfMisSNPs.append(snpToAdd)
            elif info[:4] == "Must":
                this.listOfMusts = SNP.MustList(info[5:], this.name)
                this.geneTag = 'I'
            else: #temp[:3] == "Nonsense" 
                snpToAdd = SNP.SNP_Non(this.translated, info[9:], this.name)
                if(snpToAdd.isValid()):
                    this.listOfNonsenseSNPs.append(snpToAdd)
        if (this.geneTag == 'N') or (this.geneTag == 'F'):
            for i in range(0,12):
                this.outputInfo.update({i:0})
        elif (this.geneTag == 'I'):
            for i in range(0,10):
                this.outputInfo.update({i:0})
        elif (this.geneTag == 'S'):
            for i in range(0,14):
                this.outputInfo.update({i:0})
        else:
            for i in range(0,13):
                this.outputInfo.update({i:0})

    def currentReadHasNonstop(this):
        return this.currentReadNonstop
    def foundNonstop(this, foundNonstop):
        this.currentReadNonstop = foundNonstop
    def getOutputInfo(this):
        return this.outputInfo
    def clearOutputInfo(this):
        for key in this.outputInfo:
            this.outputInfo[key] = 0
        this.additionalInfo.clear()
        this.longIndel = 0
    def getGeneTag(this):
        return this.geneTag
    def addToOutputInfo(this, index):
        this.outputInfo[index] += 1
    def addDetails(this, read, info):
        if (type(info) == tuple):
            if 'Mult:' in info[-1]:
                deletionCount = 0
                for mt in info[0]:
                    if '-' == mt[2]:
                        deletionCount +=1
                if deletionCount >= 4:
                    this.longIndel += 1
            elif len(info[0]) >= 4:
                this.longIndel += 1
        if (len(this.additionalInfo) > 0) and (read is this.additionalInfo[-1][0]):
            temp = list(this.additionalInfo[-1])
            temp.append(info)
            this.additionalInfo[-1] = tuple(temp)
        else:
            this.additionalInfo.append((read, info))
    def getLastTupleInfo(this):
        return this.additionalInfo[-1]
    def redefineLastTupleInfo(this, read):
        if len(this.additionalInfo) == 0:
            this.additionalInfo.append((read.query_name, None))
        elif read is not this.additionalInfo[-1][0]:
            this.additionalInfo.append((read.query_name, None))
        else:
            toRedefine = this.additionalInfo.pop()
            toRedefine = list(toRedefine)
            toRedefine[0] = toRedefine[0].query_name
            for i in range(1,len(toRedefine)):
                if type(toRedefine[i]) == tuple:
                    toRedefine[i] = toRedefine[i][-1]
                elif toRedefine[i] == "hypersusceptible":
                    toRedefine[i] = "Hypersusceptible: " + toRedefine[i][-1][5:]
            this.additionalInfo.append(toRedefine)
    def mustSuppressFrameshift(this):
        if this.geneTag != 'S':
            return False
        elif len(this.additionalInfo) == 0:
            return False
        elif type(this.additionalInfo[-1][0]) == string:
            return False
        elif 'Suppressible C insert' not in this.additionalInfo[-1][1:]:
            return False
        return True
    def writeAdditionalInfo(this, file):
        header = {"Read": None}
        if not(this.rRna()):
            header.update({"FS till end": None})
            header.update({"Newly found nonsense": None})
            header.update({"12+bp indel": None})
            header.update({"12+bp frameshift": None})
        else:
            header.update({"12+bp indel": None})
        if this.geneTag == 'I':
            header.update({"All residues in query": None})
            header.update({"Some residues in query": None})
            header.update({"No residues in query": None})
            header.update({"Mutated residues in query": None})
        elif this.geneTag == 'S':
            header.update({"C insert + not SRTR": None})
            header.update({"C insert + not SRTRPR": None})
            header.update({"C insert followed by del/ins": None})
            header.update({"Suppressible C insert": None})
        condensedInfo = this.condensedHyperInfo()
        if len(condensedInfo) != 0:
            header.update({"Hypersusceptible:" + condensedInfo[-1][5:]: None})
        condensedInfo = this.condensedMisInDelInfo()
        if len(condensedInfo) != 0:
            for geneMtInfo in condensedInfo:
                header.update({geneMtInfo[-1]: None})
        condensedInfo = this.condensedMultInfo()
        if len(condensedInfo) != 0:
            for geneMtInfo in condensedInfo:
                header.update({geneMtInfo[-1]: None})
        condensedInfo = this.condensedNonInfo()
        if len(condensedInfo) != 0:
            for geneMtInfo in condensedInfo:
                header.update({geneMtInfo[-1]: None})
        if this.isThereANonstop[0] == True:
            snp = this.isThereANonstop[1]
            if snp == None:
                header.update({"Nonstop": None})
            else:
                header.update({"Nonstop + " + snp.condensedInfo()[-1]: None})
        
        def clearHeader():
            for key in header:
                header[key] = None

        comma = False
        for key in header:
            if comma:
                file.write(",")
            else:
                comma = True
            file.write(key)

        for read in this.additionalInfo:
            if read[1] == None: continue
            file.write("\n")
            header["Read"] = read[0]

            for info in read[1:]:
                if info == "All":
                    header["All residues in query"] = "T"
                elif info == "Some":
                    header["Some residues in query"] = "T"
                elif info == "NA":
                    header["No residues in query"] = "T"
                elif info == "Mutant":
                    header["Mutated residues in query"] = "T"
                elif info == "nonstop":
                    snp = this.isThereANonstop[1]
                    if snp == None:
                        header["Nonstop"] = "T"
                    else:
                        header["Nonstop + " + snp.condensedInfo()[-1]] = "T"
                elif "Newly found nonsense" in info:
                    header["Newly found nonsense"] = "Pos:" + info.split(":")[1]
                elif "12+bp indel" in info:
                    count = int(info.split(":")[1][1:]) - this.longIndel
                    if count > 0:
                        header["12+bp indel"] = "Count: " + str(count)
                elif "12+bp frameshift" in info:
                    header["12+bp frameshift"] = "Count:" + info.split(":")[1]
                elif "Hypersusceptible" in info:
                    header["Hypersusceptible:" + this.condensedHyperInfo()[-1][5:]] = "T"
                else:
                    header[info] = "T"
            comma = False
            for value in header.values():
                if comma:
                    file.write(",")
                else:
                    comma = True
                if value == None:
                    continue
                file.write(value)
            clearHeader()
            

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
        if len(this.listOfHyperSNPs) == 0:
            return []
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
    def getFrameshiftInfo(this):
        return this.frameshiftInfo
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