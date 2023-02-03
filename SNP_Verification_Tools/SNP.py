from . import establishContext

class SNP:
    def __init__(this, mutable_mutant_string, name):
        this.name = name
        mutant_string = mutable_mutant_string[0]
        this.snp = mutant_string
        this.wildtype_base = mutant_string[:1]
        mutant_string = mutant_string[1:]
        this.position = 0
        while (mutant_string[0].isdigit()):
            this.position = this.position*10 + int(mutant_string[0])
            mutant_string = mutant_string[1:]
        mutable_mutant_string[0] = mutant_string
        this.left_context = list()
        this.right_context = list()

    # Because of discrepancies between megares sequence and source sequence, it is necessary to check
    # the actual wild-type residue/base and the actual positon using the context
    def findACT(this, sequence, rRNA):
        this.wildtype_base_ACT = ""
        this.position_ACT = -1
        # A 4-residue deletion in megares sequence requires hard-coding this scenario
        if (this.name == "MEG_4057") and (this.position == 347):
            this.wildtype_base_ACT = this.wildtype_base
            this.position_ACT = this.position
        else:
            # Actual mutation needs to be thirty base pairs/residues 
            # away from recorded mutation position
            begin = this.position - 30 - len(this.left_context)
            end = this.position + 31  - len(this.left_context)
            if end > (len(sequence) - 1 - len(this.right_context) - len(this.left_context)):
                end = len(sequence) - 1 - len(this.right_context) - len(this.left_context)
                begin = end - 61
            elif begin < 0:
                begin = 0
                end = 61
            for i in range(begin, end+1):
                if this.checkLeft(0, i, sequence, 0, rRNA):
                    this.changeACT(sequence, i)
                    break
    def checkLeft(this, index, currentPos, sequence, errorMargin, rRNA):
        if (len(this.left_context) != 0):
            for aa in this.left_context[index]:
                if aa == sequence[currentPos]:
                    if index == len(this.left_context) - 1:
                        return this.checkRight(0, currentPos + 2, sequence, errorMargin, rRNA)
                    else:
                        return this.checkLeft(index + 1, currentPos + 1, sequence, errorMargin, rRNA)
        else:
            return this.checkRight(0, currentPos + 2, sequence, errorMargin, rRNA)
        if (errorMargin < 3) and not(rRNA):
            if index == len(this.left_context) - 1:
                return this.checkRight(0, currentPos + 2, sequence, errorMargin + 1, rRNA)
            else:
                return this.checkLeft(index + 1, currentPos + 1, sequence, errorMargin + 1, rRNA)
        else:
            return False
    def checkRight(this, index, currentPos, sequence, errorMargin, rRNA):
        if (len(this.right_context) != 0):
            for aa in this.right_context[index]:
                if aa == sequence[currentPos]:
                    if index == len(this.right_context) - 1:
                        return True
                    else:
                        return this.checkRight(index + 1, currentPos + 1, sequence, errorMargin, rRNA)
        else:
            return True
        if (errorMargin < 3) and not(rRNA):
            if index == len(this.right_context) - 1:
                return True
            else:
                return this.checkRight(index + 1, currentPos + 1, sequence, errorMargin + 1, rRNA)
        else:
            return False
    def isValid(this):
        return this.position_ACT > -1

class Missense(SNP):
    def __init__(this, sequence, mutant_string, name, rRNA = False):        # What mutant_string looks like:
        mutable_mutant_string = [mutant_string]                                 # W123MUTANT_ABCDE_FGHIJ
        SNP.__init__(this, mutable_mutant_string, name)
        mutant_string = mutable_mutant_string[0]                                # MUTANT_ABCDE_FGHIJ
        this.mutant_list = list(mutant_string[:mutant_string.find('_')])
        mutant_string = mutant_string[mutant_string.find('_')+1:]               # ABCDE_FGHIJ
        establishContext(this, mutant_string, this.left_context, this.right_context)
        SNP.findACT(this, sequence, rRNA)
    def condensedInfo(this):
        wild_type = this.wildtype_base
        if (this.wildtype_base != this.wildtype_base_ACT):
            wild_type += this.wildtype_base_ACT
        return (wild_type, this.position_ACT, this.mutant_list, "Mis:" + this.snp)
    def changeACT(this, sequence, i):
        try: #if sequence contains mt
            this.mutant_list.index(sequence[i+len(this.leftleft_contextContext)])
            this.wildtype_base_ACT = this.wildtype_base
            this.position_ACT = i+len(this.left_context)+1
        except ValueError: 
            this.wildtype_base_ACT = sequence[i+len(this.left_context)]
            this.position_ACT = i+len(this.left_context)+1

class Nonsense(SNP):
    def __init__(this, sequence, mutant_string, name):                      # What mutant_string looks like:
        mutable_mutant_string = [mutant_string]                                 # W123*_ABCDE_FGHIJ
        SNP.__init__(this, mutable_mutant_string, name)
        mutant_string = mutable_mutant_string[0]                                # *_ABCDE_FGHIJ
        mutant_string = mutant_string[2:]                                       # ABCDE_FGHIJ
        establishContext(this, mutant_string, this.left_context, this.right_context)
        SNP.findACT(this, sequence, False)
    def condensedInfo(this):
        wild_type = this.wildtype_base
        if (this.wildtype_base != this.wildtype_base_ACT):
            wild_type += this.wildtype_base_ACT
        return (wild_type, this.position_ACT, "*", "Nonsense:" + this.snp)
    def changeACT(this, sequence, i):
        this.wildtype_base_ACT = sequence[i+len(this.left_context)]
        this.position_ACT = i+len(this.left_context)+1