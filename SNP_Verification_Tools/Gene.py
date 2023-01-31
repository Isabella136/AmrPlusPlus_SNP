from . import SNP
from . import InDel
from . import dnaTranslate
import csv


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
    #Base class that contains same variables as an N-type rRNA
    def __init__(this, name, sequence):
        this.name = name.split('|')[0]
        this.full_name = name + "|RequiresSNPConfirmation"
        this.sequence = sequence.upper()

        #Create lists required by Gene class
        this.list_of_missenses = list()
        this.list_of_ntuples = list()     
        this.list_of_deletions = list()
        this.list_of_insertions = list()

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
        return False

    def createAdditionalInfoHeader(this, header):
        header.update({"12+bp indel": None})

    def updateAdditionalInfoHeader(this, header, info, read):
        if "12+bp indel" in info:
            if ((this.geneTag == 'F') and (read[-1] == "sus")) or ((this.geneTag != 'F') and (read[-1] == "res")):
                count = int(info.split(":")[1][1:]) - this.long_indel
                if count > 0:
                    header["12+bp indel"] = "Count: " + str(count)
        else:
            header[info] = "T"

    def writeAdditionalInfo(this, file):
        header = {"Read": None}
        this.createAdditionalInfoHeader(this, header)

        condensed_info = this.condensedMisInDelInfo()
        if len(condensed_info) != 0:
            for geneMtInfo in condensed_info:
                header.update({geneMtInfo[-1]: None})
        condensed_info = this.condensedMultInfo()
        if len(condensed_info) != 0:
            for geneMtInfo in condensed_info:
                header.update({geneMtInfo[-1]: None})
        
        def clearHeader():
            for key in header:
                header[key] = ''

        csvwriter = csv.writer(file, delimiter = ',')
        csvwriter.writerow(list(header.keys()))

        for read in this.additional_info:
            if read[1] == None: continue
            header["Read"] = read[0]
            for info in read[1:]:
                if info in ["res", "sus"]: continue
                else: this.updateAdditionalInfoHeader(header, info, read)
            csvwriter.writerow(list(header.values()))
            clearHeader()

    def aaSequence(this):
        return None
    def ntSequence(this):
        return this.sequence

    def condensedMultInfo(this):
        condensed_info_list = []
        for snp in this.list_of_ntuples :
            condensed_info_list.append(snp.condensedInfo())
        return condensed_info_list
    def condensedMisInDelInfo(this):
        condensed_info_list = []
        for snp in this.list_of_missenses :
            condensed_info_list.append(snp.condensedInfo())
        for snp in this.list_of_insertions:
            condensed_info_list.append(snp.condensedInfo())
        for snp in this.list_of_deletions :
            condensed_info_list.append(snp.condensedInfo())
        return condensed_info_list
    def getNonstopInfo(this):
        return None
    def getFrameshiftInfo(this):
        return None
    def getFirstMustBetweenParams(this, begin, end):
        return None

    def getName(this):
        return this.name
    def getFullName(this):
        return this.full_name
    def ntSeqLength(this):
        return len(this.sequence)
    def rRna(this):
        return True

class Protein(Gene):
    def __init__(this, name, sequence):
        this.list_of_nonsense = []
        this.translated = dnaTranslate(sequence, name)

    def condensedNonInfo(this):
        condensed_info_list = []
        for snp in this.list_of_nonsense :
            condensed_info_list.append(snp.condensedInfo())
        return condensed_info_list
    def aaSequence(this):
        return this.translated
    def rRna(this):
        return False
    def createAdditionalInfoHeader(this, header):
        header.update({"FS till end": None})
        header.update({"Newly found nonsense": None})
        header.update({"12+bp indel": None})
        header.update({"12+bp frameshift": None})
        condensed_info = this.condensedNonInfo()
        if len(condensed_info) != 0:
            for geneMtInfo in condensed_info:
                header.update({geneMtInfo[-1]: None})
    def updateAdditionalInfoHeader(this, header, info, read):
        if "Newly found nonsense" in info:
            header["Newly found nonsense"] = "Pos:" + info.split(":")[1]
        elif "12+bp frameshift" in info:
            if ((this.geneTag == 'F') and (read[-1] == "sus")) or ((this.geneTag != 'F') and (read[-1] == "res")):
                header["12+bp frameshift"] = "Count:" + info.split(":")[1]
        else: Gene.updateAdditionalInfoHeader(this, header, info, read)

class Normal(Gene):
    def __init__(this, name, sequence):
        Gene.__init__(this, name, sequence)
        this.output_info = [0]*12
        this.tag = 'N'

class NormalrRNA(Normal):
    def __init__(this, name, sequence, info_string):
        Normal.__init__(this, name, sequence)
        infoList = info_string.split('|')
        for info in infoList:
            #Current gene is an rRNA and current variant has multiple mutations
            if info[:7] == "NucMult":   addToList(info[8:], this.sequence, this.name, "Nuclear Ntuple", this.list_of_ntuples)
            #Current gene is an rRNA and current variant has a deletion
            elif info[:6] == "NucDel":  addToList(info[7:], this.sequence, this.name, "Nuclear Deletion", this.list_of_deletions)
            #Current gene is an rRNA and current variant has a point mutation
            elif info[:3] == "Nuc":     addToList(info[4:], this.sequence, this.name, "Nuclear Missense", this.list_of_missenses)

            else:   raise NotImplementedError("{} identified as N-type rRNA but contains unrecognized variant".format(this.name))

    def createAdditionalInfoHeader(this, header):
        Gene.createAdditionalInfoHeader(this, header)

class NormalProtein(Protein, Normal):
    def __init__(this, name, sequence, info_string):
        Normal.__init__(this, name, sequence)
        this.nonstop_info = [False, None]

        #Will only be true for MEG_6142 in the case of a stop codon at pos 26
        this.special_FS_allowed = False
        this.read_has_special_FS = False

        infoList = info_string.split('|')
        for info in infoList:

            #If currently listed gene variant has a nonstop mutation
            if "Nonstop" in info:
                #If currently listed gene variant also has a missense mutation
                newInfo = info[9:info.find(";")] if info[:4] == "Mult" else None
                mutation_type = "Nonstop Missense" if info[:4] == "Mult" else "Nonstop"
                addToList(newInfo, this.translated, this.name, mutation_type, this.nonstop_info)

            #If currently listed gene variant has multiple mutations
            elif info[:4] == "Mult":    addToList(info[5:], this.translated, this.name, "Ntuple", this.list_of_ntuples)
            #If currently listed gene variant has a missense mutation
            elif info[:3] == "Mis":     addToList(info[4:], this.translated, this.name, "Missense", this.list_of_missenses)
            #IF currently listed gene variant has a deletion
            elif info[:3] == "Del":     addToList(info[4:], this.translated, this.name, "Deletion", this.list_of_deletions)
            #If currently listed gene variant has an insertion
            elif info[:3] == "Ins":     addToList(info[4:], this.translated, this.name, "Insertion", this.list_of_insertions)
            #If current gene variant has a nonsense mutation
            elif info[:8] == "Nonsense":addToList(info[9:], this.translated, this.name, "Nonsense", this.list_of_nonsense)
            #If currently listed gene variant has a frameshift
            elif info[:3] == "FS-":     this.special_FS_allowed = True

            else:   raise NotImplementedError("{} identified as N-type protein but contains unrecognized variant".format(this.name))

    def updateAdditionalInfoHeader(this, header, info, read):
        if info == "nonstop":
            snp = this.nonstop_info[1]
            if snp == None:
                header["Nonstop"] = "T"
            else:
                header["Nonstop + " + snp.condensedInfo()[-1]] = "T"
        else: Protein.updateAdditionalInfoHeader(this, header, info, read)
    def hasSpecialCase(this):
        this.read_has_special_FS = True
    def currentReadSpecial(this):
        return this.read_has_special_FS
    def resetForNextRead(this):
        this.current_read_nonstop = None
        this.read_has_special_FS = False
    def getNonstopInfo(this):
        return this.nonstop_info
    def createAdditionalInfoHeader(this, header):
        Protein.createAdditionalInfoHeader(this, header)
        if this.nonstop_info[0] == True:
            snp = this.nonstop_info[1]
            if snp == None:
                header.update({"Nonstop": None})
            else:
                header.update({"Nonstop + " + snp.condensedInfo()[-1]: None})


class Frameshift(Protein, Gene):
    def __init__(this, name, sequence, info_string):
        Gene.__init__(this, name, sequence)
        this.frameshift_info = None
        
        this.output_info = [0]*12
        this.tag = 'F'

        infoList = info_string.split('|')
        for info in infoList:
            #If currently listed gene variant has multiple mutations
            if info[:4] == "Mult":      addToList(info[5:], this.translated, this.name, "Ntuple", this.list_of_ntuples)
            #If currently listed gene variant has a missense mutation
            elif info[:3] == "Mis":     addToList(info[4:], this.translated, this.name, "Missense", this.list_of_missenses)
            #IF currently listed gene variant has a deletion
            elif info[:3] == "Del":     addToList(info[4:], this.translated, this.name, "Deletion", this.list_of_deletions)
            #If currently listed gene variant has an insertion
            elif info[:3] == "Ins":     addToList(info[4:], this.translated, this.name, "Insertion", this.list_of_insertions)
            #If current gene variant has a nonsense mutation
            elif info[:8] == "Nonsense":addToList(info[9:], this.translated, this.name, "Nonsense", this.list_of_nonsense)
            #If currently listed gene variant has a frameshift
            elif info[:3] == "FS-":     this.frameshift_info = info[3:]

            else:   raise NotImplementedError("{} identified as F-type protein but contains unrecognized variant".format(this.name))

    def getFrameshiftInfo(this):
        return this.frameshift_info

    def createAdditionalInfoHeader(this, header):
        Protein.createAdditionalInfoHeader(this, header)

class Suppressible(Protein, Gene):
    def __init__(this, name, sequence, info_string):
        Gene.__init__(this, name, sequence)
        this.output_info = [0]*14
        this.tag = 'S'

        infoList = info_string.split('|')
        for info in infoList[1:]:
            #If currently listed gene variant has multiple mutations
            if info[:4] == "Mult":      addToList(info[5:], this.translated, this.name, "Ntuple", this.list_of_ntuples)
            #If currently listed gene variant has a missense mutation
            elif info[:3] == "Mis":     addToList(info[4:], this.translated, this.name, "Missense", this.list_of_missenses)
            #IF currently listed gene variant has a deletion
            elif info[:3] == "Del":     addToList(info[4:], this.translated, this.name, "Deletion", this.list_of_deletions)
            #If currently listed gene variant has an insertion
            elif info[:3] == "Ins":     addToList(info[4:], this.translated, this.name, "Insertion", this.list_of_insertions)
            #If current gene variant has a nonsense mutation
            elif info[:8] == "Nonsense":addToList(info[9:], this.translated, this.name, "Nonsense", this.list_of_nonsense)

            else:   raise NotImplementedError("{} identified as S-type protein but contains unrecognized variant".format(this.name))
    
    def mustSuppressFrameshift(this):
        if len(this.additional_info) == 0:
            return False
        elif type(this.additional_info[-1][0]) == str:
            return False
        elif 'Suppressible C insert' not in this.additional_info[-1][1:]:
            return False
        return True

    def createAdditionalInfoHeader(this, header):
        Protein.createAdditionalInfoHeader(this, header)
        header.update({"C insert + not SRTR": None})
        header.update({"C insert + not SRTRPR": None})
        header.update({"C insert followed by del/ins": None})
        header.update({"Suppressible C insert": None})

class Hypersusceptible(Protein, Gene):
    def __init__(this, name, sequence, info_string):
        Gene.__init__(this, name, sequence)
        this.output_info = [0]*13
        this.tag = 'H'

        infoList = info_string.split('|')
        for info in infoList:
            #If currently listed gene variant has multiple mutations
            if info[:4] == "Mult":      addToList(info[5:], this.translated, this.name, "Ntuple", this.list_of_ntuples)
            #If currently listed gene variant has a missense mutation
            elif info[:3] == "Mis":     addToList(info[4:], this.translated, this.name, "Missense", this.list_of_missenses)
            #IF currently listed gene variant has a deletion
            elif info[:3] == "Del":     addToList(info[4:], this.translated, this.name, "Deletion", this.list_of_deletions)
            #If currently listed gene variant has an insertion
            elif info[:3] == "Ins":     addToList(info[4:], this.translated, this.name, "Insertion", this.list_of_insertions)
            #If current gene variant has a nonsense mutation
            elif info[:8] == "Nonsense":addToList(info[9:], this.translated, this.name, "Nonsense", this.list_of_nonsense)
            #If currently listed gene variant is hypersusceptible
            elif info[:5] == "Hyper":   addToList(info[6:], this.translated, this.name, "Hypersusceptible", this.list_of_hypersusceptible_snps)

            else:   raise NotImplementedError("{} identified as H-type protein but contains unrecognized variant".format(this.name))

    def condensedHyperInfo(this):
        if len(this.list_of_hypersusceptible_snps) == 0:
            return list()
        return this.list_of_hypersusceptible_snps[0].condensedInfo()

    def createAdditionalInfoHeader(this, header):
        Protein.createAdditionalInfoHeader(this, header)
        condensed_info = this.condensedHyperInfo()
        if len(condensed_info) != 0:
            header.update({"Hypersusceptible:" + condensed_info[-1][5:]: None})
            this.condensed_hyper_info = condensed_info[-1][5:]

    def updateAdditionalInfoHeader(this, header, info, read):
        if "Hypersusceptible" in info:
            header["Hypersusceptible:" + this.condensedHyperInfo()[-1][5:]] = "T"
        else:
            Protein.updateAdditionalInfoHeader(this, header, info, read)

class Intrinsic(Gene):
    def __init__(this, name, sequence):
        Gene.__init__(this, name, sequence)
        this.list_of_musts = None
        this.output_info = [0]*10
        this.tag = 'I'

    def getFirstMustBetweenParams(this, begin, end):
        return this.list_of_musts.getFirstMustBetweenParams(begin, end)

    def createAdditionalInfoHeader(this, header):
        header.update({"All residues in query": None})
        header.update({"Some residues in query": None})
        header.update({"No residues in query": None})
        header.update({"Mutated residues in query": None})

    def updateAdditionalInfoHeader(this, header, info, read):
        if info == "All":       header["All residues in query"] = "T"
        elif info == "Some":    header["Some residues in query"] = "T"
        elif info == "NA":      header["No residues in query"] = "T"
        elif info == "Mutant":  header["Mutated residues in query"] = "T"


class IntrinsicrRNA(Intrinsic):
    def __init__(this, name, sequence, info_string):
        Gene.__init__(this, name, sequence)

        infoList = info_string.split('|')
        for info in infoList:
            #Current gene is an rRNA and current variant has multiple mutations
            if info[:7] == "NucMult":   addToList(info[8:], this.sequence, this.name, "Nuclear Ntuple", this.list_of_ntuples)
            #Current gene is an rRNA and current variant has a deletion
            elif info[:6] == "NucDel":  addToList(info[7:], this.sequence, this.name, "Nuclear Deletion", this.list_of_deletions)
            #Current gene is an rRNA and current variant has a point mutation
            elif info[:3] == "Nuc":     addToList(info[4:], this.sequence, this.name, "Nuclear Missense", this.list_of_missenses)
            #If current gene is intrisically resistant when it is current variant
            elif info[:4] == "Must":    this.list_of_musts = SNP.MustList(info[5:], this.name)

            else:   raise NotImplementedError("{} identified as I-type rRNA but contains unrecognized variant".format(this.name))

    def createAdditionalInfoHeader(this, header):
        Intrinsic.createAdditionalInfoHeader(this, header)
        Gene.createAdditionalInfoHeader(this, header)


class IntrinsicProtein(Protein, Intrinsic):
    def __init__(this, name, sequence, info_string):
        Gene.__init__(this, name, sequence)

        infoList = info_string.split('|')
        for info in infoList:
            #If currently listed gene variant has multiple mutations
            if info[:4] == "Mult":      addToList(info[5:], this.translated, this.name, "Ntuple", this.list_of_ntuples)
            #If currently listed gene variant has a missense mutation
            elif info[:3] == "Mis":     addToList(info[4:], this.translated, this.name, "Missense", this.list_of_missenses)
            #IF currently listed gene variant has a deletion
            elif info[:3] == "Del":     addToList(info[4:], this.translated, this.name, "Deletion", this.list_of_deletions)
            #If currently listed gene variant has an insertion
            elif info[:3] == "Ins":     addToList(info[4:], this.translated, this.name, "Insertion", this.list_of_insertions)
            #If current gene variant has a nonsense mutation
            elif info[:8] == "Nonsense":addToList(info[9:], this.translated, this.name, "Nonsense", this.list_of_nonsense)
            #If current gene is intrisically resistant when it is current variant
            elif info[:4] == "Must":    this.list_of_musts = SNP.MustList(info[5:], this.name)

            else:   raise NotImplementedError("{} identified as I-type protein but contains unrecognized variant".format(this.name))

    def createAdditionalInfoHeader(this, header):
        Intrinsic.createAdditionalInfoHeader(this, header)
        Protein.createAdditionalInfoHeader(this, header)


