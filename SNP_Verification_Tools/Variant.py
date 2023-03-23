from . import SNP, InDel
from . import establishContext


class MutatedVariant:

    # Description of input:
    # 
    #   sequence:                       Gene sequence; will be used for context finding
    #   variant_string:                 Variant information as found in SNPInfo.fasta
    #   variant_type:                   Information on type of mutation(s) present in the variant; is equal to 
    #                                   either 'Ntuple', 'Missense', 'Deletion', 'Insertion', 'Nonsense'
    #   name:                           Has the form MEG_XXXX
    #   rRNA:                           Boolean that is True if and only if the gene is an rRNA.

    def __init__(this, sequence, variant_string, variant_type, name, rRNA = False):
        this.variant = variant_string
        this.list_of_mutations = list()                     # Used if this is an N-tuple mutation
        this.single_mutation = None                         # Used if this is a single mutation

        this.deletion_count = 0                             # Makes sure that known variant with long deletion or insertion
        this.insertion_count = 0                            # won't be marked as special in detailed_output

        if 'Ntuple' in variant_type:
            prelim_list_of_mutations = variant_string.split(";")
            for info in prelim_list_of_mutations:
                if (info[:3] == "Mis") or (info[:3] == "Nuc"):
                    this.list_of_mutations.append(SNP.Missense(sequence, info[4:], name, rRNA))
                elif info[:3] == "Ins":
                    this.list_of_mutations.append(InDel.Insertion(sequence, info[4:], name, rRNA))
                elif info[:3] == "Del":
                    this.list_of_mutations.append(InDel.Deletion(sequence, info[4:], name, rRNA))
                    this.deletion_count += 1
                else:
                    raise NotImplementedError("{} contains N-tuple with unrecognized mutation type".format(name))
        elif 'Missense' in variant_type:
            this.single_mutation = SNP.Missense(sequence, this.variant, name, rRNA)
        elif 'Deletion' in variant_type:
            this.single_mutation = InDel.Deletion(sequence, this.variant, name, rRNA)
        elif 'Insertion' in variant_type:
            this.single_mutation = InDel.Insertion(sequence, this.variant, name, rRNA)
            this.insertion_count = this.single_mutation.getInsertionCount()
        elif 'Nonsense' in variant_type:
            this.single_mutation = SNP.Nonsense(sequence, this.variant, name)    
            
    def condensedInfo(this):
        if len(this.list_of_mutations) == 0:
            return this.single_mutation.condensedInfo()
        toReturn = []
        for snp in this.list_of_mutations:
            toReturn.append(snp.condensedInfo())
        return (toReturn, "Mult:" + this.variant)
    def isValid(this):
        if len(this.list_of_mutations) == 0:
            return this.single_mutation.isValid()
        for snp in this.list_of_mutations:
            if not(snp.isValid()):
                return False
        return True
    def longIndel(this):
        return this.deletion_count >= 4 or this.insertion_count >= 4

# Wild-type base or residue that is required for intrinsic resistance
class Must:
    def __init__(this, wildtype_string, name):
        this.name = name
        this.wildtype_string = wildtype_string
        this.wildtype_base = wildtype_string[:1]
        this.wildtype_pos = int(wildtype_string[1:wildtype_string.find('_')])
        this.left_context = list()
        this.right_context = list()
        establishContext(wildtype_string, this.left_context, this.right_context)
        this.next = None
    def condensedInfo(this):
        return (this.wildtype_base, this.wildtype_pos)
    def getPos(this):
        return this.wildtype_pos
    def defineNext(this,next_must):
        this.next = next_must
    def getNext(this):
        return this.next
    def getWt(this):
        return this.wildtype_base

# Original gene allele that already provides AMR; i.e. gene is intrinsically resistant.
# Considered as a SNPConfirmation gene because this allele contains key bases or residues that must be present for AMR.
# This means that a SNP or an InDel at any of those bases/residues can potentially lead to a loss in resistance.

class IntrinsicVariant:
    def __init__(this, variant_string, name, rRNA = False):
        this.variant = variant_string
        this.list_of_musts = dict()
        variant_info_list = this.variant.split(';')
        for i in range(len(variant_info_list)): 
            variant_info_list[i] = variant_info_list[i][4 if rRNA else 6:]
            to_add =  Must(variant_info_list[i], name)
            if len(this.list_of_musts) > 0:
                list(this.list_of_musts.values())[-1].defineNext(to_add)
            else:
                this.first_pos = to_add.getPos()
            this.list_of_musts[to_add.getPos()] = to_add
        this.last_pos = list(this.list_of_musts.keys())[-1]

    # Description of input:
    #   begin:                      Position of first base/residue of read aligned to megares
    #   end:                        Position of last base/residue of read aligned to megares
    # 
    # Description of potential output:
    #   If no 'Must" was found:     None
    #   If 'Must' was found:        Tuple of that contained the following two elements:
    #                                   1.  First 'Must' object located between begin and end
    #                                   2.  Boolean indicating whether it is the first 'Must' of the variant

    def getFirstMustBetweenParams(this, begin, end):
        if (end < this.first_pos) or (begin > this.last_pos):
            return None
        elif (end >= this.first_pos) and (begin <= this.first_pos): 
            return (this.list_of_musts[this.first_pos], True)
        for pos in this.list_of_musts.keys():
            if (begin <= pos) and (end >= pos):
                return (this.list_of_musts[pos], False)
        return None

