from . import SNP, InDel

class MutatedVariant:
    def __init__(this, sequence, variant_string, variant_type, name, rRNA):
        this.variant = variant_string
        this.list_of_mutations = list()         #If this is an N-tuple mutation
        this.single_mutation = None             #If this is a single mutation
        this.deletionCount = 0                  #Makes sure that known variant with long deletion 
                                                #won't be marked as special in detailed_output
        if 'Ntuple' in variant_type:
            prelim_list_of_mutations = variant_string.split(";")
            for info in prelim_list_of_mutations:
                if (info[:3] == "Mis") or (info[:3] == "Nuc"):
                    this.list_of_mutations.append(SNP.SNP_Mis(sequence, info[4:], name, rRNA))
                elif info[:3] == "Ins":
                    this.list_of_mutations.append(InDel.Insertion(sequence, info[4:], name, rRNA))
                elif info[:3] == "Del":
                    this.list_of_mutations.append(InDel.Deletion(sequence, info[4:], name, rRNA))
                    this.deletionCount += 1
                else:
                    raise NotImplementedError("{} contains N-tuple with unrecognized mutation type".format(name))
        elif 'Missense' in variant_type:
            this.single_mutation = SNP.SNP_Mis(sequence, this.variant, name, rRNA)
        elif 'Deletion' in variant_type:
            this.single_mutation = InDel.Deletion(sequence, this.variant, name, rRNA)
        elif 'Insertion' in variant_type:
            this.single_mutation = InDel.Insertion(sequence, this.variant, name, rRNA)
        elif 'Nonsense' in variant_type:
            this.single_mutation = SNP.SNP_Non(sequence, this.variant, name)    
            
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
        return this.deletionCount >= 4

class IntrinsicVariant:
    def __init__(this, sequence, variant_string, variant_type, name, rRNA):