mycoplasmataceaeList = ["MEG_3185", "MEG_3240", "MEG_5328", "MEG_5330", "MEG_5331"]

aa = {
    "TTT" : 'F',    "TCT" : 'S',    "TAT" : 'Y',    "TGT" : 'C',
    "TTC" : 'F',    "TCC" : 'S',    "TAC" : 'Y',    "TGC" : 'C',
    "TTA" : 'L',    "TCA" : 'S',    "TAA" : '*',    "TGA" : '*',
    "TTG" : 'L',    "TCG" : 'S',    "TAG" : '*',    "TGG" : 'W',
    "CTT" : 'L',    "CCT" : 'P',    "CAT" : 'H',    "CGT" : 'R',
    "CTC" : 'L',    "CCC" : 'P',    "CAC" : 'H',    "CGC" : 'R',
    "CTA" : 'L',    "CCA" : 'P',    "CAA" : 'Q',    "CGA" : 'R',
    "CTG" : 'L',    "CCG" : 'P',    "CAG" : 'Q',    "CGG" : 'R',
    "ATT" : 'I',    "ACT" : 'T',    "AAT" : 'N',    "AGT" : 'S',
    "ATC" : 'I',    "ACC" : 'T',    "AAC" : 'N',    "AGC" : 'S',
    "ATA" : 'I',    "ACA" : 'T',    "AAA" : 'K',    "AGA" : 'R',
    "ATG" : 'M',    "ACG" : 'T',    "AAG" : 'K',    "AGG" : 'R',
    "GTT" : 'V',    "GCT" : 'A',    "GAT" : 'D',    "GGT" : 'G',
    "GTC" : 'V',    "GCC" : 'A',    "GAC" : 'D',    "GGC" : 'G',
    "GTA" : 'V',    "GCA" : 'A',    "GAA" : 'E',    "GGA" : 'G',
    "GTG" : 'V',    "GCG" : 'A',    "GAG" : 'E',    "GGG" : 'G'
}

def resistant(name, increment, argInfoDict):
    argInfo = argInfoDict.get(name, False)
    if (argInfo == False):
        argInfoDict.update({name:(increment, 1)})
    else:
        temp = list(argInfo)
        temp[0] += increment
        temp[1] += 1
        argInfo = tuple(temp)
        argInfoDict.update({name:argInfo})

def disregard(name, argInfoDict):
    resistant(name, 0, argInfoDict)

def dnaTranslate(dna, name):
    toReturn = ""
    for i in range(0, len(dna), 3):
        if dna[i:i+3].find('N') == -1:
            if (name in mycoplasmataceaeList) and (dna[i:i+3] == "TGA"):
                toReturn += 'W'
            else:
                toReturn += aa[dna[i:i+3]]
        else:
            toReturn += 'N'
    return toReturn

def reverseTranslation(amino_acid):
    if amino_acid == '*':
        return "TAA"
    for codon, aAcid in aa.items():
        if amino_acid == aAcid:
            return codon
    return ""

# Description of input:
#   context_list:                   comes from the 5 bases/residues either left or right of the mutation 
#   current_list_index:             index of current base/residue in context_list (starts at 0 when called from lookThroughContext())
#   current_sequence_index:         index of corresponding base/residue in sequence
#   sequence:                       sequence of gene whose variant is being checked
#   error_margin:                   number of mismatch found between sequence and context 
#                                   (starts at 0 when called from lookThroughContext() the first time)
#
# Returns tuple of last sequence index looked at by method and final tally of error_margin for left context

def checkContext(context_list, current_list_index, current_sequence_index, sequence, error_margin):
    if (len(context_list) == current_list_index):
        return (current_sequence_index-1, error_margin)
    for aa in context_list[current_list_index]:
        if aa == sequence[current_sequence_index]:
            return checkContext(context_list, current_list_index + 1, current_sequence_index + 1, sequence, error_margin)
    return checkContext(context_list, current_list_index + 1, current_sequence_index + 1, sequence, error_margin + 1)


# Description of input (in addition to what was described for checkContext):
#   deleted:                        True if current mutation is a deletion; otherwise False
#   position_list:                  If current mutation is an insertion or deletion of a repeated residue, 
#                                   then multiple positions can be called as having the inserted/deleted residue by SAMtools
#   rRNA:                           True if gene is rRNA; otherwise False
#
#
# Returns True if 1. it is an rRNA and context has perfect match, or 2. it is a protein, and it has at least a 70% match

def lookThroughContext(current_sequence_index, sequence, left_context, right_context, inserted = False, position_list = None, rRNA = False):
    current_sequence_index, error_margin = checkContext(left_context, 0, current_sequence_index, sequence, 0)
    if (error_margin/(len(left_context)+len(right_context)) > (3/10)) or (rRNA and error_margin > 0):
        return False
    current_sequence_index += 1 if position_list == None else (position_list[-1] - position_list[0] + 1)
    current_sequence_index += 0 if inserted else 1
    current_sequence_index, error_margin = checkContext(right_context, 0, current_sequence_index, sequence, error_margin)
    if (error_margin/(len(left_context)+len(right_context)) > (3/10)) or (rRNA and error_margin > 0):
        return False
    return True

# Description of input:
#   context_string:                 Retrieved from SNPinfo.fasta, provides up to 5 bases/residues to the left and to the right of the mutation
#   left_context:                   Mutable list that belongs to Must, SNP, or InDel
#   right_context:                  Mutable list that belongs to Must, SNP, or InDel

def establishContext(context_string, left_context, right_context):
    context_list = context_string.split("_")
    go_to_next = True
    for char in context_list[0]:
        if char == '[':
            go_to_next = False
            left_context.append([])
        elif char == ']':
            go_to_next = True
        else:
            if go_to_next:
                left_context.append([char])
            else:
                left_context[-1].append(char)
    for char in context_list[1]:
        if char == '[':
            go_to_next = False
            right_context.append([])
        elif char == ']':
            go_to_next = True
        else:
            if go_to_next:
                right_context.append([char])
            else:
                right_context[-1].append(char)