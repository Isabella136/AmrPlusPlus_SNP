from . import establishContext
from . import lookThroughContext

# Superclass for Insertiona and Deletion
# 
# Contains the following local variable:
#   name:                           Has the form MEG_XXXX
#   position_list:                  Positions of all deletions or insertions
#   left_context:                   List of the 5 bases or residues that are located to the left of the mutation
#   right_context:                  List of the 5 bases or residues that are located to the right of the mutation
#   wildtype_base_ACT:              The original base (or residue) according to the MEGARes sequence
#   position_ACT:                   The position of the mutation according to the MEGARes sequence

class InDel:
    def __init__(this, mutant_string, name):
        this.name = name
        position_string = mutant_string[:mutant_string.find("_")]
        mutant_string = mutant_string[mutant_string.find("_")+1:]
        this.position_list = position_string.split("/")
        this.position_list = [int(i) for i in this.position_list]
        this.left_context = list()
        this.right_context = list()
        establishContext(mutant_string, this.left_context, this.right_context)

    def changeACT(this, sequence_index):
        previous = this.position_list[0]                            # Previous indel postion
        difference_from_first = 0                                   # Difference from first indel position
        for current in this.position_list:                          # Current indel position
            difference_from_first += (current - previous)
            previous = current
            this.position_list_ACT.append(1 +                       # Left-most position changed from 0 to 1
                sequence_index +                                    # Index of first base/residue in left_context
                len(this.left_context) +                            # Length of left_context
                difference_from_first)                              # Difference from first indel position
    def isValid(this):
        return len(this.position_list_ACT) != 0

class Insertion(InDel):
    def __init__(this, sequence, mutant_string, name, rRNA = False):
        this.indel = mutant_string
        this.inserted = list()
        while not(mutant_string[0].isdigit()):
            this.inserted.append(mutant_string[0])
            mutant_string = mutant_string[1:]
        InDel.__init__(this, mutant_string, name)
        this.findACT(sequence, rRNA)
    def findACT(this, sequence, rRNA):
        this.position_list_ACT = []
        # Actual insertion needs to be thirty base pairs/residues 
        # away from recorded insertion position
        begin = this.position_list[0] - 30 - len(this.left_context)
        end = this.position_list[-1] + 30 - len(this.left_context)
        if begin < 0:
            begin = 0
            end = begin + 60
        elif end > (len(sequence) - 1 - (this.position_list[-1]-this.position_list[0]) - len(this.left_context) - len(this.right_context)):
            end = len(sequence) - 1 - (this.position_list[-1]-this.position_list[0]) - len(this.left_context) - len(this.right_context)
            begin = end - 60
        for sequence_index in range(begin, end+1):
            if lookThroughContext(sequence_index, sequence, this.left_context, this.right_context, rRNA=rRNA, inserted = True, position_list=this.position_list):
                this.changeACT(sequence_index)
                break
    def condensedInfo(this):
        return (this.inserted, this.position_list_ACT, "+", "Ins:" + this.indel)
    def getInsertionCount(this):
        return len(this.inserted)

class Deletion(InDel):
    def __init__(this, sequence, mutant_string, name, rRNA = False):
        this.indel = mutant_string
        this.deleted = mutant_string[0]
        InDel.__init__(this, mutant_string[1:], name)
        this.findACT(sequence, rRNA)
    def changeACT(this, sequence, sequence_index):
        this.deletionACT = sequence[sequence_index+len(this.left_context)]
        InDel.changeACT(this, sequence_index)
    def findACT(this, sequence, rRNA):
        this.position_list_ACT = []
        # Actual deletion needs to be thirty base pairs/residues 
        # away from recorded deletion position
        begin = this.position_list[0] - 30 - len(this.left_context)
        end = this.position_list[-1] + 31  - len(this.left_context)
        if begin < 0:
            begin = 0
            end = begin + 61
        elif end > (len(sequence) - 2 - (this.position_list[-1]-this.position_list[0]) - len(this.left_context) - len(this.right_context)):
            end = len(sequence) - 2 - (this.position_list[-1]-this.position_list[0]) - len(this.left_context) - len(this.right_context)
            begin = end - 61
        for sequence_index in range(begin, end+1):
            if lookThroughContext(sequence_index, sequence, this.left_context, this.right_context, rRNA=rRNA, position_list=this.position_list):
                this.changeACT(sequence, sequence_index)
                break
    def isValid(this):
        return (len(this.position_list_ACT) != 0) and (this.deleted == this.deletionACT)
    def condensedInfo(this):
        return (this.deleted, this.position_list_ACT, "-", "Del:" + this.indel)