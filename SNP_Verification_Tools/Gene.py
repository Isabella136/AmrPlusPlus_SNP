from . import SNP
from . import InDel
from . import dnaTranslate


def addToList(info, sequence, name, mutation_type, mutable_list):
    rRNA = True if 'Nuclear' in mutation_type else False

    if 'Nonstop' in mutation_type:
        mutable_list[0] = True
        if 'Missense' in mutation_type:
            to_add = SNP.SNP_Mis(sequence, info, name)
            mutable_list[1] = to_add if to_add.isValid() else None

    elif 'Ntuple' in mutation_type:
        to_add = SNP.SNP_Mult(sequence, info, name)
        if to_add.isValid(): mutable_list.append(to_add)

    elif 'Missense' in mutation_type:
        to_add = SNP.SNP_Mis(sequence, info, name, rRNA)
        if to_add.isValid(): mutable_list.append(to_add)

    elif 'Deletion' in mutation_type:
        to_add = InDel.Deletion(sequence, info, name, rRNA)
        if to_add.isValid(): mutable_list.append(to_add)

    elif 'Insertion' in mutation_type:
        to_add = InDel.Insertion(sequence, info, name)
        if to_add.isValid(): mutable_list.append(to_add)

    elif 'Hypersusceptible' in mutation_type:
        to_add = SNP.SNP_Mult(sequence, info, name)
        if to_add.isValid(): mutable_list.append(to_add)

    elif 'Nonsense' in mutation_type:
        to_add = SNP.SNP_Non(sequence, info, name)
        if to_add.isValid(): mutable_list.append(to_add)

class Gene:
    def __init__(this, name, sequence):
        this.name = name.split('|')[0]
        this.full_name = name + "|RequiresSNPConfirmation"
        this.sequence = sequence.upper()

        #Create lists required by Gene class
        this.list_of_missenses = list()
        this.list_of_ntuples = list()     
        this.list_of_deletions = list()
        this.additional_info = list()
        this.output_info = list()
        this.tag = None
        
        this.long_indel = 0
        this.current_read_nonstop = None

    def currentReadHasNonstop(this):
        return this.current_read_nonstop
    def foundNonstop(this, found_nonstop):
        this.current_read_nonstop = found_nonstop
    def hasSpecialCase(this):
        pass
    def currentReadSpecial(this):
        return False
    def resetForNextRead(this):
        this.current_read_nonstop = None
    def getOutputInfo(this):
        return this.output_info
    def clearOutputInfo(this):
        for index in range(this.output_info):
            this.output_info[index] = 0
        this.additional_info.clear()
    def getGeneTag(this):
        return this.tag
    def addToOutputInfo(this, index):
        this.output_info[index] += 1
        if index == 0:
            this.long_indel = 0
    def addDetails(this, read, info):
        if (type(info) == tuple):
            if 'Mult:' in info[-1]:
                deletion_count = 0
                for mt in info[0]:
                    if '-' == mt[2]:
                        deletion_count +=1
                if deletion_count >= 4:
                    this.long_indel += 1
            elif len(info[0]) >= 4:
                this.long_indel += 1
        if (len(this.additional_info) > 0) and ((read is this.additional_info[-1][0]) or isinstance(read, str)):
            temp = list(this.additional_info[-1])
            temp.append(info)
            this.additional_info[-1] = tuple(temp)
        else:
            this.additional_info.append((read, info))
    def getLastTupleInfo(this):
        return this.additional_info[-1]
    def redefineLastTupleInfo(this, read):
        if len(this.additional_info) == 0:
            this.additional_info.append((read.query_name, None))
        elif read is not this.additional_info[-1][0]:
            this.additional_info.append((read.query_name, None))
        else:
            to_redefine = this.additional_info.pop()
            to_redefine = list(to_redefine)
            to_redefine[0] = to_redefine[0].query_name
            for i in range(1,len(to_redefine)):
                if type(to_redefine[i]) == tuple:
                    to_redefine[i] = to_redefine[i][-1]
                elif to_redefine[i] == "hypersusceptible":
                    to_redefine[i] = "Hypersusceptible: " + to_redefine[i][-1][5:]
            this.additional_info.append(to_redefine)
    def mustSuppressFrameshift(this):
        if this.geneTag != 'S':
            return False
        elif len(this.additional_info) == 0:
            return False
        elif type(this.additional_info[-1][0]) == str:
            return False
        elif 'Suppressible C insert' not in this.additional_info[-1][1:]:
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
        condensed_info = this.condensedHyperInfo()
        if len(condensed_info) != 0:
            header.update({"Hypersusceptible:" + condensed_info[-1][5:]: None})
        condensed_info = this.condensedMisInDelInfo()
        if len(condensed_info) != 0:
            for geneMtInfo in condensed_info:
                header.update({geneMtInfo[-1]: None})
        condensed_info = this.condensedMultInfo()
        if len(condensed_info) != 0:
            for geneMtInfo in condensed_info:
                header.update({geneMtInfo[-1]: None})
        condensed_info = this.condensedNonInfo()
        if len(condensed_info) != 0:
            for geneMtInfo in condensed_info:
                header.update({geneMtInfo[-1]: None})
        if this.nonstop_info[0] == True:
            snp = this.nonstop_info[1]
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

        for read in this.additional_info:
            if read[1] == None: continue
            file.write("\n")
            header["Read"] = read[0]

            for info in read[1:]:
                if info in ["res", "sus"]:
                    continue
                elif info == "All":
                    header["All residues in query"] = "T"
                elif info == "Some":
                    header["Some residues in query"] = "T"
                elif info == "NA":
                    header["No residues in query"] = "T"
                elif info == "Mutant":
                    header["Mutated residues in query"] = "T"
                elif info == "nonstop":
                    snp = this.nonstop_info[1]
                    if snp == None:
                        header["Nonstop"] = "T"
                    else:
                        header["Nonstop + " + snp.condensedInfo()[-1]] = "T"
                elif "Newly found nonsense" in info:
                    header["Newly found nonsense"] = "Pos:" + info.split(":")[1]
                elif "12+bp indel" in info:
                    if ((this.geneTag == 'F') and (read[-1] == "sus")) or ((this.geneTag != 'F') and (read[-1] == "res")):
                        count = int(info.split(":")[1][1:]) - this.long_indel
                        if count > 0:
                            header["12+bp indel"] = "Count: " + str(count)
                elif "12+bp frameshift" in info:
                    if ((this.geneTag == 'F') and (read[-1] == "sus")) or ((this.geneTag != 'F') and (read[-1] == "res")):
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
        condensed_info_list = []
        for snp in this.list_of_ntuples :
            condensed_info_list.append(snp.condensedInfo())
        return condensed_info_list
    def condensedHyperInfo(this):
        if len(this.listOfHyperSNPs) == 0:
            return []
        return this.listOfHyperSNPs[0].condensedInfo()
    def condensedMisInDelInfo(this):
        condensed_info_list = []
        for snp in this.list_of_missenses :
            condensed_info_list.append(snp.condensedInfo())
        for snp in this.list_of_insertions:
            condensed_info_list.append(snp.condensedInfo())
        for snp in this.list_of_deletions :
            condensed_info_list.append(snp.condensedInfo())
        return condensed_info_list
    def condensedNonInfo(this):
        condensed_info_list = []
        for snp in this.listOfNonsenseSNPs :
            condensed_info_list.append(snp.condensedInfo())
        return condensed_info_list
    def getNonstopInfo(this):
        return this.nonstop_info
    def getFrameshiftInfo(this):
        return this.frameshiftInfo
    def getFirstMustBetweenParams(this, begin, end):
        if this.listOfMusts == None:
            return None
        return this.listOfMusts.getFirstMustBetweenParams(begin, end)

    def getName(this):
        return this.name
    def getFullName(this):
        return this.full_name
    def ntSeqLength(this):
        return len(this.sequence)
    def aaOrNu(this):
        return this.listOfMusts.returnAaOrNu()
    def rRna(this):
        return False

class Normal(Gene):
    def __init__(this, name, sequence):
        Gene.__init__(this, name, sequence)
        this.output_info = [0]*12
        this.tag = 'N'

class rRNA(Normal):
    def __init__(this, name, sequence, info_string):
        Normal.__init__(this, name, sequence)
        infoList = info_string.split('|')
        for info in infoList:

            #Current gene is an rRNA and current variant has multiple mutations
            if info[:7] == "NucMult": addToList(info[8:], this.sequence, this.name, "Nuclear Ntuple", this.list_of_ntuples)

            #Current gene is an rRNA and current variant has a deletion
            elif info[:6] == "NucDel": addToList(info[7:], this.sequence, this.name, "Nuclear Deletion", this.list_of_deletions)

            #Current gene is an rRNA and current variant has a point mutation
            elif info[:3] == "Nuc": addToList(info[4:], this.sequence, this.name, "Nuclear Missense", this.list_of_missenses)

    def rRna(this):
        return True

class Protein(Normal):
    def __init__(this, name, sequence, info_string):
        Normal.__init__(this, name, sequence)
        this.translated = dnaTranslate(this.sequence, name)
        this.nonstop_info = [False, None]
        this.list_of_insertions = []
        this.list_of_nonsense = []

        #Will only be true for MEG_6142 in the case of a stop codon at pos 26
        this.special_FS_allowed = False
        this.read_has_special_FS = False

        infoList = info_string.split('|')
        for info in infoList:

            #If currently listed gene variant has a nonstop mutation
            if "Nonstop" in info:

                #If currently listed gene variant also has a missense mutation
                if info[:4] == "Mult":
                    snpToAdd = SNP.SNP_Mis(this.translated, info[9:info.find(";")], this.name)
                    if (snpToAdd.isValid()):
                        this.nonstop_info = (True, snpToAdd)

                #If currently listed gene variant doesn't need any other mutation
                else:
                    this.nonstop_info = (True, None)

            #If currently listed gene variant has multiple mutations
            elif info[:4] == "Mult":
                snpToAdd = SNP.SNP_Mult(this.translated, info[5:], this.name)
                if(snpToAdd.isValid()):
                    this.list_of_ntuples.append(snpToAdd)

            #If currently listed gene variant has a missense mutation
            elif info[:3] == "Mis":
                snpToAdd = SNP.SNP_Mis(this.translated, info[4:], this.name)
                if(snpToAdd.isValid()):
                    this.list_of_missenses.append(snpToAdd)

            #IF currently listed gene variant has a deletion
            elif info[:3] == "Del":
                indelToAdd = InDel.Deletion(this.translated, info[4:], this.name)
                if(indelToAdd.isValid()):
                    this.list_of_deletions.append(indelToAdd)

            #If currently listed gene variant has an insertion
            elif info[:3] == "Ins":
                indelToAdd = InDel.Insertion(this.translated, info[4:], this.name)
                if(indelToAdd.isValid()):
                    this.list_of_insertions.append(indelToAdd)

            #If currently listed gene variant has a frameshift
            elif info[:3] == "FS-":
                this.special_FS_allowed = True

            #If current gene variant has a nonsense mutation
            else:
                snpToAdd = SNP.SNP_Non(this.translated, info[9:], this.name)
                if(snpToAdd.isValid()):
                    this.listOfNonsenseSNPs.append(snpToAdd)

    def hasSpecialCase(this):
        this.read_has_special_FS = True
    def currentReadSpecial(this):
        return this.read_has_special_FS
    def resetForNextRead(this):
        this.current_read_nonstop = None
        this.read_has_special_FS = False

class Frameshift(Gene):
    def __init__(this, name, sequence, info_string):
        Gene.__init__(this, name, sequence)
        this.translated = dnaTranslate(this.sequence, name)
        this.frameshiftInfo = None
        this.list_of_insertions = []
        this.listOfNonsenseSNPs = []
        this.output_info = [0]*12
        this.tag = 'F'

        infoList = info_string.split('|')
        for info in infoList:

            #If currently listed gene variant has multiple mutations
            if info[:4] == "Mult":
                snpToAdd = SNP.SNP_Mult(this.translated, info[5:], this.name)
                if(snpToAdd.isValid()):
                    this.list_of_ntuples.append(snpToAdd)

            #If currently listed gene variant has a missense mutation
            elif info[:3] == "Mis":
                snpToAdd = SNP.SNP_Mis(this.translated, info[4:], this.name)
                if(snpToAdd.isValid()):
                    this.list_of_missenses.append(snpToAdd)

            #IF currently listed gene variant has a deletion
            elif info[:3] == "Del":
                indelToAdd = InDel.Deletion(this.translated, info[4:], this.name)
                if(indelToAdd.isValid()):
                    this.list_of_deletions.append(indelToAdd)

            #If currently listed gene variant has an insertion
            elif info[:3] == "Ins":
                indelToAdd = InDel.Insertion(this.translated, info[4:], this.name)
                if(indelToAdd.isValid()):
                    this.list_of_insertions.append(indelToAdd)

            #If currently listed gene variant has a frameshift
            elif info[:3] == "FS-":
                this.frameshiftInfo = info[3:]
                if (this.name != "MEG_6142"):
                    this.geneTag = 'S' if this.name == "MEG_6094" else 'F'

            #If current gene variant has a nonsense mutation
            else:
                snpToAdd = SNP.SNP_Non(this.translated, info[9:], this.name)
                if(snpToAdd.isValid()):
                    this.listOfNonsenseSNPs.append(snpToAdd)

class Suppressible(Gene):
    def __init__(this, name, sequence, info_string):
        Gene.__init__(this, name, sequence)
        this.translated = dnaTranslate(this.sequence, name)
        this.list_of_insertions = []
        this.listOfNonsenseSNPs = []
        this.output_info = [0]*14
        this.tag = 'S'

        infoList = info_string.split('|')
        for info in infoList[1:]:

            #If currently listed gene variant has multiple mutations
            if info[:4] == "Mult":
                snpToAdd = SNP.SNP_Mult(this.translated, info[5:], this.name)
                if(snpToAdd.isValid()):
                    this.list_of_ntuples.append(snpToAdd)

            #If currently listed gene variant has a missense mutation
            elif info[:3] == "Mis":
                snpToAdd = SNP.SNP_Mis(this.translated, info[4:], this.name)
                if(snpToAdd.isValid()):
                    this.list_of_missenses.append(snpToAdd)

            #IF currently listed gene variant has a deletion
            elif info[:3] == "Del":
                indelToAdd = InDel.Deletion(this.translated, info[4:], this.name)
                if(indelToAdd.isValid()):
                    this.list_of_deletions.append(indelToAdd)

            #If currently listed gene variant has an insertion
            elif info[:3] == "Ins":
                indelToAdd = InDel.Insertion(this.translated, info[4:], this.name)
                if(indelToAdd.isValid()):
                    this.list_of_insertions.append(indelToAdd)

            #If current gene variant has a nonsense mutation
            else:
                snpToAdd = SNP.SNP_Non(this.translated, info[9:], this.name)
                if(snpToAdd.isValid()):
                    this.listOfNonsenseSNPs.append(snpToAdd)

class Hypersusceptible(Gene):
    def __init__(this, name, sequence, info_string):
        Gene.__init__(this, name, sequence)
        this.translated = dnaTranslate(this.sequence, name)
        this.listOfHyperSNPs = []
        this.list_of_insertions = []
        this.listOfNonsenseSNPs = []
        this.output_info = [0]*13
        this.tag = 'H'

        infoList = info_string.split('|')
        for info in infoList:

            #If currently listed gene variant has multiple mutations
            if info[:4] == "Mult":
                snpToAdd = SNP.SNP_Mult(this.translated, info[5:], this.name)
                if(snpToAdd.isValid()):
                    this.list_of_ntuples.append(snpToAdd)

            #If currently listed gene variant has a missense mutation
            elif info[:3] == "Mis":
                snpToAdd = SNP.SNP_Mis(this.translated, info[4:], this.name)
                if(snpToAdd.isValid()):
                    this.list_of_missenses.append(snpToAdd)

            #IF currently listed gene variant has a deletion
            elif info[:3] == "Del":
                indelToAdd = InDel.Deletion(this.translated, info[4:], this.name)
                if(indelToAdd.isValid()):
                    this.list_of_deletions.append(indelToAdd)

            #If currently listed gene variant has an insertion
            elif info[:3] == "Ins":
                indelToAdd = InDel.Insertion(this.translated, info[4:], this.name)
                if(indelToAdd.isValid()):
                    this.list_of_insertions.append(indelToAdd)

            #If currently listed gene variant is hypersusceptible
            elif info[:5] == "Hyper":
                snpToAdd = SNP.SNP_Mult(this.translated, info[6:], this.name)
                if(snpToAdd.isValid()):
                    this.listOfHyperSNPs.append(snpToAdd)

            #If current gene variant has a nonsense mutation
            else:
                snpToAdd = SNP.SNP_Non(this.translated, info[9:], this.name)
                if(snpToAdd.isValid()):
                    this.listOfNonsenseSNPs.append(snpToAdd)

class Intrinsic(Gene):
    def __init__(this, name, sequence):
        Gene.__init__(this, name, sequence)
        this.listOfMusts = None
        this.output_info = [0]*10
        this.tag = 'I'


class rRNA(Intrinsic):
    def __init__(this, name, sequence, info_string):
        Gene.__init__(this, name, sequence)

        infoList = info_string.split('|')
        for info in infoList:

            #If current gene is an rRNA and current variant has multiple mutations
            if info[:7] == "NucMult":
                snpToAdd = SNP.SNP_Mult(this.sequence, info[8:], this.name)
                if(snpToAdd.isValid()):
                    this.list_of_ntuples.append(snpToAdd)

            #If current gene is an rRNA and current variant has a deletion
            elif info[:6] == "NucDel":
                indelToAdd = InDel.Deletion(this.sequence, info[7:], this.name, True)
                if(indelToAdd.isValid()):
                    this.list_of_deletions.append(indelToAdd)

            #If current gene is an rRNA and current variant has a point mutation
            elif info[:3] == "Nuc":
                snpToAdd = SNP.SNP_Mis(this.sequence, info[4:], this.name, True)
                if(snpToAdd.isValid()):
                    this.list_of_missenses.append(snpToAdd)

            #If current gene is intrisically resistant when it is current variant
            elif info[:4] == "Must":
                this.listOfMusts = SNP.MustList(info[5:], this.name)
                this.geneTag = 'I'
    def rRna(this):
        return True

class Protein(Intrinsic):
    def __init__(this, name, sequence, info_string):
        Gene.__init__(this, name, sequence)
        this.translated = dnaTranslate(this.sequence, name)
        this.list_of_insertions = []

        infoList = info_string.split('|')
        for info in infoList:

            #If currently listed gene variant has multiple mutations
            if info[:4] == "Mult":
                snpToAdd = SNP.SNP_Mult(this.translated, info[5:], this.name)
                if(snpToAdd.isValid()):
                    this.list_of_ntuples.append(snpToAdd)

            #If currently listed gene variant has a missense mutation
            elif info[:3] == "Mis":
                snpToAdd = SNP.SNP_Mis(this.translated, info[4:], this.name)
                if(snpToAdd.isValid()):
                    this.list_of_missenses.append(snpToAdd)

            #IF currently listed gene variant has a deletion
            elif info[:3] == "Del":
                indelToAdd = InDel.Deletion(this.translated, info[4:], this.name)
                if(indelToAdd.isValid()):
                    this.list_of_deletions.append(indelToAdd)

            #If currently listed gene variant has an insertion
            elif info[:3] == "Ins":
                indelToAdd = InDel.Insertion(this.translated, info[4:], this.name)
                if(indelToAdd.isValid()):
                    this.list_of_insertions.append(indelToAdd)

            #If current gene is intrisically resistant when it is current variant
            elif info[:4] == "Must":
                this.listOfMusts = SNP.MustList(info[5:], this.name)
                this.geneTag = 'I'
            
            #If current gene variant has a nonsense mutation
            else:
                snpToAdd = SNP.SNP_Non(this.translated, info[9:], this.name)
                if(snpToAdd.isValid()):
                    this.listOfNonsenseSNPs.append(snpToAdd)


