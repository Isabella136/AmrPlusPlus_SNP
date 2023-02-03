from . import establishContext

class InDel:
    def __init__(this, mutant_string, name, insertion):
        this.indel = mutant_string
        this.name = name
        this.inserted = list()
        this.deleted = None
        if insertion:
            while not(mutant_string[0].isdigit()):
                this.inserted.append(mutant_string[0])
                mutant_string = mutant_string[1:]
        else:
            this.deleted = mutant_string[0]
        mutant_string = mutant_string[1:]
        position_string = mutant_string[:mutant_string.find("_")]
        mutant_string = mutant_string[mutant_string.find("_")+1:]
        this.position_list = position_string.split("/")
        this.position_list = [int(i) for i in this.position_list]
        this.left_context = list()
        this.right_context = list()
        establishContext(this, mutant_string, this.left_context, this.right_context)
    def checkLeft(this, current_list_index, current_sequence_index, sequence, error_margin, rRNA = False):
        if (len(this.left_context) != 0):
            for aa in this.left_context[current_list_index]:
                if aa == sequence[current_sequence_index]:
                    if current_list_index == len(this.left_context) - 1:
                        nextPos = current_sequence_index + (this.position_list[-1]-this.position_list[0]) + 1
                        nextPos = nextPos if this.inserted != None else nextPos+1
                        return InDel.checkRight(this, 0, nextPos , sequence, error_margin, rRNA)
                    else:
                        return InDel.checkLeft(this, current_list_index + 1, current_sequence_index + 1, sequence, error_margin, rRNA)
        else:
            nextPos = current_sequence_index + (this.position_list[-1]-this.position_list[0]) + 1
            nextPos = nextPos if this.inserted != None else nextPos+1
            return this.checkRight(0, nextPos, sequence, error_margin, rRNA)
        if (error_margin < 3) and not(rRNA):
            if current_list_index == len(this.left_context) - 1:
                nextPos = current_sequence_index + (this.position_list[-1]-this.position_list[0]) + 1
                nextPos = nextPos if this.inserted != None else nextPos+1
                return InDel.checkRight(this, 0, nextPos, sequence, error_margin + 1, rRNA)
            else:
                return InDel.checkLeft(this, current_list_index + 1, current_sequence_index + 1, sequence, error_margin + 1, rRNA)
        else:
            return False
    def checkRight(this, current_list_index, current_sequence_index, sequence, error_margin, rRNA):
        if (len(this.right_context) != 0):
            for aa in this.right_context[current_list_index]:
                if aa == sequence[current_sequence_index]:
                    if current_list_index == len(this.right_context) - 1:
                        return True
                    else:
                        return InDel.checkRight(this, current_list_index + 1, current_sequence_index + 1, sequence, error_margin, rRNA)
        else:
            return True
        if (error_margin < 3) and not(rRNA):
            if current_list_index == len(this.right_context) - 1:
                return True
            else:
                return InDel.checkRight(this, current_list_index + 1, current_sequence_index + 1, sequence, error_margin + 1, rRNA)
        else:
            return False
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
    def __init__(this, sequence, mtString, name, rRNA = False):
        InDel.__init__(this, mtString, name, True)
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
        elif end > (len(sequence) - (this.position_list[-1]-this.position_list[0]) - len(this.left_context) - len(this.right_context)):
            end = len(sequence) - (this.position_list[-1]-this.position_list[0]) - len(this.left_context) - len(this.right_context)
            begin = end - 60
        for sequence_index in range(begin, end+1):
            if this.checkLeft(0, sequence_index, sequence, 0, rRNA):
                this.changeACT(sequence_index)
                break
    def condensedInfo(this):
        return (this.inserted, this.position_list_ACT, "+", "Ins:" + this.indel)
    def getInsertionCount(this):
        return len(this.inserted)

class Deletion(InDel):
    def __init__(this, sequence, mtString, name, rRNA = False):
        InDel.__init__(this, mtString, name, False)
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
        elif end > (len(sequence) - 1 - (this.position_list[-1]-this.position_list[0]) - len(this.left_context) - len(this.right_context)):
            end = len(sequence) - 1 - (this.position_list[-1]-this.position_list[0]) - len(this.left_context) - len(this.right_context)
            begin = end - 61
        for sequence_index in range(begin, end+1):
            if this.checkLeft(0, sequence_index, sequence, 0, rRNA):
                this.changeACT(sequence, sequence_index)
                break
    def isValid(this):
        return (len(this.position_list_ACT) != 0) and (this.deleted == this.deletionACT)
    def condensedInfo(this):
        return (this.deleted, this.position_list_ACT, "-", "Del:" + this.indel)