from . import establishContext
from . import lookThroughContext

# Superclass for Missense and Nonsense
# 
# Contains the following local variable:
#   name:                           Has the form MEG_XXXX
#   snp:                            String rendition of the information gathered from SNPinfo.fasta
#   wildtype_base:                  Original base (or residue) as gathered from either the literature or other databases
#   position:                       Position of the mutation
#   left_context:                   List of the 5 bases or residues that are located to the left of the mutation
#   right_context:                  List of the 5 bases or residues that are located to the right of the mutation
#   wildtype_base_ACT:              The original base (or residue) according to the MEGARes sequence
#   position_ACT:                   The position of the mutation according to the MEGARes sequence

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

    # Finds the original base/residue and the position of the mutation based on the MEGARes sequence
    def findACT(this, sequence, rRNA):
        this.wildtype_base_ACT = ""
        this.position_ACT = -1

        # A 4-residue deletion in megares sequence requires hard-coding this scenario
        if (this.name == "MEG_4057") and (this.position == 347):
            this.wildtype_base_ACT = this.wildtype_base
            this.position_ACT = this.position
        
        # Position of the mutation in MEGARes needs to be thirty base pairs/residues away from recorded mutation position
        else:
            begin = this.position - 30 - len(this.left_context)
            end = this.position + 31  - len(this.left_context)
            if end > (len(sequence) - 2 - len(this.right_context) - len(this.left_context)):
                end = len(sequence) - 2 - len(this.right_context) - len(this.left_context)
                begin = end - 61
            elif begin < 0:
                begin = 0
                end = 61
            for sequence_index in range(begin, end+1):
                if lookThroughContext(sequence_index, sequence, this.left_context, this.right_context, rRNA=rRNA):
                    this.changeACT(sequence, sequence_index)
                    break
    def isValid(this):
        return this.position_ACT > -1



# Subclass that is reserved for missense mutations
#
# Other local variables:
#   mutant_list:                    List of all residues or bases that can replace the original residue/base and confer resistance.

class Missense(SNP):
    def __init__(this, sequence, mutant_string, name, rRNA = False):        # What mutant_string looks like:
        mutable_mutant_string = [mutant_string]                                 # W123MUTANT_ABCDE_FGHIJ
        SNP.__init__(this, mutable_mutant_string, name)
        mutant_string = mutable_mutant_string[0]                                # MUTANT_ABCDE_FGHIJ
        this.mutant_list = list(mutant_string[:mutant_string.find('_')])
        mutant_string = mutant_string[mutant_string.find('_')+1:]               # ABCDE_FGHIJ
        establishContext(mutant_string, this.left_context, this.right_context)
        SNP.findACT(this, sequence, rRNA)
    def condensedInfo(this):
        wild_type = this.wildtype_base
        if (this.wildtype_base != this.wildtype_base_ACT):
            wild_type += this.wildtype_base_ACT
        return (wild_type, this.position_ACT, this.mutant_list, "Mis:" + this.snp)
    def changeACT(this, sequence, i):
        try: #if sequence contains mt
            this.mutant_list.index(sequence[i+len(this.left_context)])
            this.wildtype_base_ACT = this.wildtype_base
            this.position_ACT = i+len(this.left_context)+1
        except ValueError: 
            this.wildtype_base_ACT = sequence[i+len(this.left_context)]
            this.position_ACT = i+len(this.left_context)+1



# Subclass that is reserved for nonsense mutations

class Nonsense(SNP):
    def __init__(this, sequence, mutant_string, name):                      # What mutant_string looks like:
        mutable_mutant_string = [mutant_string]                                 # W123*_ABCDE_FGHIJ
        SNP.__init__(this, mutable_mutant_string, name)
        mutant_string = mutable_mutant_string[0]                                # *_ABCDE_FGHIJ
        mutant_string = mutant_string[2:]                                       # ABCDE_FGHIJ
        establishContext(mutant_string, this.left_context, this.right_context)
        SNP.findACT(this, sequence, False)
    def condensedInfo(this):
        wild_type = this.wildtype_base
        if (this.wildtype_base != this.wildtype_base_ACT):
            wild_type += this.wildtype_base_ACT
        return (wild_type, this.position_ACT, "*", "Nonsense:" + this.snp)
    def changeACT(this, sequence, i):
        this.wildtype_base_ACT = sequence[i+len(this.left_context)]
        this.position_ACT = i+len(this.left_context)+1