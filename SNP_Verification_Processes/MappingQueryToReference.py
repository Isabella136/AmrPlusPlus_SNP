from SNP_Verification_Tools import dnaTranslate

def extendCigar(cigar):
    toReturn = ""
    count = 0
    for c in cigar:
        if c.isdigit():
            count = count * 10 + int(c)
        else:
            while count > 0:
                toReturn = toReturn + c
                count -= 1
    return toReturn

def MapQueryToReference(rRna, read, name):
    cigar = extendCigar(read.cigarstring)
    aligned_pair = read.get_aligned_pairs()
    nt_alignment_map = mapCigarToAlignment(cigar, aligned_pair, rRna)
    querySeq = read.query_sequence
    mapOfInterest = transformNtAlignmentMap(nt_alignment_map)
    start = nt_alignment_map[0][1]
    i = 1
    while start == None:
        start = nt_alignment_map[i][1]
        i += 1
    end = nt_alignment_map[len(nt_alignment_map)-1][1]
    i = 2
    while end == None:
        end = nt_alignment_map[len(nt_alignment_map)-i][1]
        i += 1
    trimmedQuerySequence = querySeq[start:end+1]
    seqOfInterest = trimmedQuerySequence
    if not(rRna):
        aaQuerySequence = dnaTranslate(trimmedQuerySequence, name)
        aa_alignment_map = aaAlignment(nt_alignment_map)
        seqOfInterest = aaQuerySequence
        mapOfInterest = aa_alignment_map
    return (seqOfInterest, mapOfInterest)

def mapCigarToAlignment(cigar, aligned_pair, rRna):
    alignment_map = []          # list of tuples that consist of opcode, query pos, ref pos, ref nuc shift, and pre ref nuc shift
    queryLength = 0
    refLength = 0
    index = 0                   # query seq index
    shift = 0                   # number that consist of shift between ref and query, determined by equation: previous total insertions - previous total deletions 
    prevShift = 0
    splitIndex = -1             # trimmed query seq index; starts counting from index in codon
    translatable = False
    for op in cigar:
        if op == "S":
            index += 1
        elif (op == "M") | (op == "X") | (op == "="):
            prevShift = shift
            if not(translatable) and not(rRna):
                if splitIndex < 0:
                    splitIndex = aligned_pair[index][1] % 3
                else:
                    splitIndex += 1
                if (splitIndex % 3) != 0:
                    index += 1
                else:
                    translatable = True
                    queryLength += 1
                    refLength += 1
                    alignment_map.append(("M", aligned_pair[index][0], aligned_pair[index][1], shift, prevShift))
                    index += 1
            else:
                if not(translatable):
                    translatable = True
                queryLength += 1
                refLength += 1
                alignment_map.append(("M", aligned_pair[index][0], aligned_pair[index][1], shift, prevShift))
                index += 1
        elif op == "I":
            prevShift = shift
            shift += 1
            if not(translatable) and not(rRna):
                if splitIndex < 0:
                    temp = aligned_pair[index][1]
                    i = 1
                    while temp == None:
                        temp = aligned_pair[index+i][1]
                        i += 1
                    splitIndex = temp % 3
                else:
                    splitIndex += 1
                if (splitIndex % 3) != 0:
                    index += 1
                else:
                    translatable = True
                    queryLength += 1
                    alignment_map.append((op, aligned_pair[index][0], aligned_pair[index][1], shift, prevShift))
                    index += 1
            else:
                if not(translatable):
                    translatable = True
                queryLength += 1
                alignment_map.append((op, aligned_pair[index][0], aligned_pair[index][1], shift, prevShift))
                index += 1
        elif op == "D":
            prevShift = shift
            shift -= 1
            if not(translatable) and not(rRna):
                if ((splitIndex % 3) != 2) and (splitIndex != -1):
                    index += 1
                else:
                    translatable = True
                    refLength += 1
                    alignment_map.append((op, aligned_pair[index][0], aligned_pair[index][1], shift, prevShift))
                    index += 1
            else:
                if not(translatable):
                    translatable = True
                refLength += 1
                alignment_map.append((op, aligned_pair[index][0], aligned_pair[index][1], shift, prevShift))
                index += 1
        else:
            continue
    if not(rRna):
        while (queryLength % 3) != 0:
            poped = alignment_map.pop()
            if poped[0] != "D" : queryLength -= 1
    return alignment_map

def transformNtAlignmentMap(nt_alignment_map):
    new_nt_alignment_map = {}
    ntQueryIndex = 0
    for nt in nt_alignment_map:
        if nt[0] == "I":
            ntQueryIndex += 1 
        elif nt[0] == "D":
            new_nt_alignment_map.update({nt[2]:(None,)})   
        else:
            new_nt_alignment_map.update({nt[2]:(ntQueryIndex,)})
            ntQueryIndex += 1    
    return new_nt_alignment_map

def aaAlignment(nt_alignment_map):
    aa_alignment_map = {

    }
    ntRefIndex = nt_alignment_map[0][2]
    i = 1
    while ntRefIndex == None:
        ntRefIndex = nt_alignment_map[i][2]
        i+=1
    ntQueryIndex = 0
    aaQueryIndex = 0
    insertCount = 0
    prevAaShift = None    
    firstAlignment = None
    thirdAlignment = 0
    mapIndex = 0
    inbetween = None
    deleteCount = 0
    hasDeletion = False
    def addToMapDeletion(index, toAdd):
        aa = aa_alignment_map.get(index-1, False)
        if aa == False:
            aa_alignment_map.update({index-1:(inbetween, )})
        else:
            temp = list(aa)
            temp.append(inbetween)
            aa = tuple(temp)
            aa_alignment_map.update({index-1:aa})
        addToMap(index, toAdd)
    def addToMapInbetween(index, toAdd):
        aa = aa_alignment_map.get(index, False)
        if aa == False:
            if (inbetween != None):
                aa_alignment_map.update({index:(inbetween, toAdd, )})
            else:
                aa_alignment_map.update({index:(toAdd,)})
        else:
            temp = list(aa)
            if (inbetween != None):
                temp.append(inbetween)
            temp.append(toAdd)
            aa = tuple(temp)
            aa_alignment_map.update({index:aa})
    def addToMap(index, toAdd):
        aa = aa_alignment_map.get(index, False)
        if aa == False:
            aa_alignment_map.update({index:(toAdd,)})
        else:
            temp = list(aa)
            temp.append(toAdd)
            aa = tuple(temp)
            aa_alignment_map.update({index:aa})
    for nt in nt_alignment_map:
        if (ntQueryIndex % 3) == 0:
            if deleteCount == 0: #1/2D3M
                hasDeletion = False
        if (ntRefIndex % 3) == 0:
            if nt[0] == "I":
                insertCount += 1
                deleteCount = 0
                if (nt[3] % 3) == 0:  
                    if insertCount == 3: #3I
                        if inbetween == None:
                            addToMap(int(ntRefIndex/3)-1, thirdAlignment)
                        inbetween = thirdAlignment
                    else: #1I...1M2I or 2I...2M1I
                        addToMapInbetween(int(ntRefIndex/3)-1, thirdAlignment)
                        inbetween = None
                    insertCount = 0
                if prevAaShift == None:
                    prevAaShift = nt[4]
                ntQueryIndex += 1
            elif nt[0] == "D":
                ntRefIndex += 1
                if prevAaShift == None:
                    prevAaShift = nt[4]
                insertCount = 0
                deleteCount += 1
                if deleteCount == 3:
                    addToMap(int(ntRefIndex/3)-1, '-')
                    inbetween = '-'
                hasDeletion = True
            else: #nt[0] == "M"
                deleteCount = 0
                ntQueryIndex += 1
                ntRefIndex += 1
                insertCount = 0
                if prevAaShift == None:
                    prevAaShift = nt[4]
                if (ntQueryIndex % 3) == 0: #2I1M3M
                    if (prevAaShift == 0):
                        inbetween = None
                    elif inbetween == None:
                        if not(hasDeletion):
                            addToMap(int(ntRefIndex/3)-1, thirdAlignment)
                        if nt[3] < 0:
                            addToMap(int(ntRefIndex/3), thirdAlignment)
                        inbetween = thirdAlignment
                    else:
                        inbetween = thirdAlignment
                firstAlignment = aaQueryIndex
        elif (ntRefIndex % 3) == 1:
            if nt[0] == "I":
                deleteCount = 0
                ntQueryIndex += 1
                if (ntQueryIndex % 3) == 0: #1I2M2M1I
                    if prevAaShift == 0:
                        inbetween = None
                    elif inbetween == None:
                        addToMap(int(ntRefIndex/3)-1, thirdAlignment)
                        inbetween = thirdAlignment
                    else:
                        inbetween = thirdAlignment
            elif nt[0] == "D":
                ntRefIndex += 1
                deleteCount += 1
                if deleteCount == 3:
                    addToMap(int(ntRefIndex/3)-1, '-')
                    inbetween = '-'
                hasDeletion = True
            else: #nt[0] == "M"
                deleteCount = 0
                ntQueryIndex += 1
                ntRefIndex += 1
                if firstAlignment == None:
                    firstAlignment = aaQueryIndex
                if prevAaShift == None:
                    prevAaShift = nt[4]
                if (ntQueryIndex % 3) == 0: #1I2M3M
                    if (prevAaShift == 0):
                        inbetween = None
                    elif inbetween == None:
                        if not(hasDeletion):
                            addToMap(int(ntRefIndex/3)-1, thirdAlignment)
                        if nt[3] < 0:
                            addToMap(int(ntRefIndex/3), thirdAlignment)
                        inbetween = thirdAlignment
                    else:
                        inbetween = thirdAlignment
        else: #(ntRefIndex % 3) == 2
            if nt[0] == "I":
                deleteCount = 0
                ntQueryIndex += 1
            elif nt[0] == "D":
                hasDeletion = True
                deleteCount += 1
                ntRefIndex += 1
                if firstAlignment != None:
                    if (prevAaShift % 3) == 0: #1M2D, 2M1D
                        addToMap(int(nt[2]/3), thirdAlignment)
                else: #3D
                    addToMap(int(nt[2]/3), '-')
                    if (prevAaShift % 3) != 0:
                        inbetween = '-'
                prevAaShift = None
                firstAlignment = None
            else: #nt[0] == "M"
                deleteCount = 0
                ntQueryIndex += 1
                ntRefIndex += 1
                if firstAlignment == None:
                    firstAlignment = aaQueryIndex
                if prevAaShift == None:
                    prevAaShift = nt[4]
                if firstAlignment == thirdAlignment:
                    if ((prevAaShift % 3) == 0) | ((nt[4] % 3) == 0): #3M, 1D1I2M, 1D/2I...2D1M, 2D/1I...1D2M, 2D3M, 1D3M
                        if inbetween != None:
                            if hasDeletion:
                                if inbetween != '-':
                                    inbetween = thirdAlignment
                                    addToMapDeletion(int(nt[2]/3), thirdAlignment)
                                else:
                                    addToMapInbetween(int(nt[2]/3), thirdAlignment)
                            else:
                                addToMapInbetween(int(nt[2]/3), thirdAlignment)
                            inbetween = None
                        else:
                            addToMap(int(nt[2]/3), thirdAlignment)
                elif (prevAaShift % 3) == 0:
                    if (nt[4] % 3) == 0: #3Mand3xI
                        if inbetween != None: #inbetween == "-"
                            aa_alignment_map.update({int(nt[2]/3):(inbetween, firstAlignment, thirdAlignment)})
                            inbetween = None
                        else:
                            aa_alignment_map.update({int(nt[2]/3):(firstAlignment, thirdAlignment)})
                    else: #1/2Iand3M
                        addToMapInbetween(int(nt[2]/3), firstAlignment)
                        inbetween = None
                elif (nt[4] % 3) == 0: #2D/1I...1M1D1M
                    if inbetween != None:
                        if hasDeletion:
                            if inbetween != '-':
                                inbetween = thirdAlignment
                                addToMapDeletion(int(nt[2]/3), thirdAlignment)
                            else:
                                addToMapInbetween(int(nt[2]/3), thirdAlignment)
                        else:
                            addToMapInbetween(int(nt[2]/3), thirdAlignment)
                        inbetween = None
                    else:
                        addToMap(int(nt[2]/3), thirdAlignment)
                prevAaShift = None
                firstAlignment = None

        aaQueryIndex = int(ntQueryIndex / 3)
        thirdAlignment = aaQueryIndex
        mapIndex += 1
    return aa_alignment_map    

