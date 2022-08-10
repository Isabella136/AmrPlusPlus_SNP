from SNP_Verification_Tools import dnaTranslate

def extendCigar(cigar):                                         #example: transforms 3M2I3M to MMMIIMMM
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
    cigar = extendCigar(read.cigarstring)                       #example: transforms 3M2I3M to MMMIIMMM
    aligned_pair = read.get_aligned_pairs()
    nt_alignment_map = mapCigarToAlignment(cigar, aligned_pair, rRna)   #map has been trimmed to not have broken codons if not rRNA
    querySeq = read.query_sequence
    mapOfInterest = transformNtAlignmentMap(nt_alignment_map)   #makes ntAlignmentMap have same format as aa_alignment_map (useful for rRNA)
    start = nt_alignment_map[0][1]
    i = 1
    while start == None:                                        #trimms query sequence in order to not have broken codons during translation
        start = nt_alignment_map[i][1]                          #for rRNA, the mapCigarToAlignment returns full map, so sequence won't be trimmed
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
        aa_alignment_map = aaAlignment(nt_alignment_map)        #for 'F' type of genes, seqOfInterest may extend past mapOfInterest
        seqOfInterest = aaQuerySequence                         #due to frameshift that extend past the end of query sequence
        mapOfInterest = aa_alignment_map
    return (seqOfInterest, mapOfInterest)

def mapCigarToAlignment(cigar, aligned_pair, rRna):             #map has been trimmed to not have broken codons if not rRNA
    alignment_map = []                                          #list of tuples: (opcode, query pos, ref pos, ref nuc shift, and pre ref nuc shift)
    queryLength = 0                                             #length of query sequence that will be kept
    refLength = 0                                               #length of reference sequence
    index = 0                                                   #query seq index
    shift = 0                                                   #previous total insertions - previous total deletions 
    prevShift = 0
    splitIndex = -1                                             #trimmed query seq index; starts counting from index in codon
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

def transformNtAlignmentMap(nt_alignment_map):                  #makes ntAlignmentMap have same format as aa_alignment_map (useful for rRNA)
    new_nt_alignment_map = {}                                   #format of map is {refNt, (correspondingQueryNt1, correspondingQueryNt2)}
    ntQueryIndex = 0
    ntRefIndex = nt_alignment_map[0][2]
    for nt in nt_alignment_map:
        if nt[0] == "I":
            if ntRefIndex not in new_nt_alignment_map:
                new_nt_alignment_map.update({ntRefIndex:tuple()})
            temp = list(new_nt_alignment_map[ntRefIndex])
            temp.append(ntQueryIndex)
            new_nt_alignment_map[ntRefIndex] = tuple(temp)
            temp = list(new_nt_alignment_map[ntRefIndex-1])
            temp.append(ntQueryIndex-1)
            new_nt_alignment_map[ntRefIndex-1] = tuple(temp)
            ntQueryIndex += 1 
        elif nt[0] == "D":
            if nt[2] not in new_nt_alignment_map:
                new_nt_alignment_map.update({nt[2]:tuple()})
            temp = list(new_nt_alignment_map[nt[2]])
            temp.append(None)
            new_nt_alignment_map[nt[2]] = tuple(temp)   
            ntRefIndex += 1
        else:
            if nt[2] not in new_nt_alignment_map:
                new_nt_alignment_map.update({nt[2]:tuple()})
            temp = list(new_nt_alignment_map[nt[2]])
            temp.append(ntQueryIndex)
            new_nt_alignment_map[nt[2]] = tuple(temp)  
            ntRefIndex += 1
            ntQueryIndex += 1    
    return new_nt_alignment_map

def aaAlignment(nt_alignment_map):                              #returns translated version of nt_alignment_map
    aa_alignment_map = {}                                       #format of map is {refCodon, (correspondingQueryCodon1, correspondingQueryCodon2)}
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
    inbetween = None                                            #if in frameshift, inbetween should either be '-' or previous query codon
    deleteCount = 0
    hasDeletion = False

    def addToMapDeletion(index, toAdd):                         #maps queryCodon in inbetween to previous referenceCodon
        aa = aa_alignment_map.get(index-1, False)               #currently, toAdd and inbetween are the same
        if aa == False:
            aa_alignment_map.update({index-1:(inbetween, )})
        else:
            temp = list(aa)
            temp.append(inbetween)
            aa = tuple(temp)
            aa_alignment_map.update({index-1:aa})
        addToMap(index, toAdd)

    def addToMapInbetween(index, toAdd):                        #adds both queryCodon and inbetweenQueryCodon to same referenceCodon
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

    def addToMap(index, toAdd):                                 #maps queryCodon to referenceCodon
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
                        if inbetween == None:                                       #second codon in MMM|III|...|MMM
                            addToMap(int(ntRefIndex/3)-1, thirdAlignment)           #maps to first in ABC|ABC
                        inbetween = thirdAlignment
                    else: #1I...1M2I or 2I...2M1I                                   #both codon in IMM|MII
                        addToMapInbetween(int(ntRefIndex/3)-1, thirdAlignment)      #maps to ABC
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
                    addToMap(int(ntRefIndex/3)-1, '-')                              #for codon MDDDMM, '-' mapped to first in ABC|ABC
                    inbetween = '-'                                                 #codon was previously mapped in else/elif D
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
                        inbetween = None                                            #codon IIM
                    elif inbetween == None:
                        if not(hasDeletion):                                        #second codon in IIM|MMM
                            addToMap(int(ntRefIndex/3)-1, thirdAlignment)           #maps to first in ABC|A...
                        if nt[3] < 0:
                            addToMap(int(ntRefIndex/3), thirdAlignment)             #first codon in MDMM|MM...
                        inbetween = thirdAlignment                                  #maps to last in ABC|ABC
                    else:
                        inbetween = thirdAlignment                                  #last codon in IIM|MMM|MMM
                firstAlignment = aaQueryIndex
        elif (ntRefIndex % 3) == 1:
            if nt[0] == "I":
                deleteCount = 0
                ntQueryIndex += 1
                if (ntQueryIndex % 3) == 0: #1I2M2M1I
                    if prevAaShift == 0:                                            #codon MII
                        inbetween = None
                    elif inbetween == None:                                         #second codon in IMM|MMI
                        addToMap(int(ntRefIndex/3)-1, thirdAlignment)               #maps to first codon in ABC|A...
                        inbetween = thirdAlignment
                    else:                                                           #last codon in IMM|MMM|MMI
                        inbetween = thirdAlignment
            elif nt[0] == "D":
                ntRefIndex += 1
                deleteCount += 1
                if deleteCount == 3:                                                #for codon MMDDDM, '-' mapped to first in ABC|ABC
                    addToMap(int(ntRefIndex/3)-1, '-')                              #codon was previously mapped in else/elif D
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
                        if not(hasDeletion):                                        #second codon in IMM|MMM
                            addToMap(int(ntRefIndex/3)-1, thirdAlignment)           #maps to first in ABC|AB...
                        if nt[3] < 0:
                            addToMap(int(ntRefIndex/3), thirdAlignment)             #first codon in MDDMM|M...
                        inbetween = thirdAlignment                                  #maps to last in ABC|ABC
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
                    if (prevAaShift % 3) == 0: #1M2D, 2M1D                          #first codon in MDDMM|M... or just MDDDMM
                        addToMap(int(nt[2]/3), thirdAlignment)                      #maps to first in ABC|ABC
                else: #3D
                    addToMap(int(nt[2]/3), '-')                                     #'-' maps to ABC
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
                                if inbetween != '-':                                #last codon in MIM|MMM|MDMM or in MDMM|DDMMM
                                    inbetween = thirdAlignment                      #maps to last two in ABC|ABC|ABC
                                    addToMapDeletion(int(nt[2]/3), thirdAlignment)
                                else:                                               #codon such as MDDDMM and '-'
                                    addToMapInbetween(int(nt[2]/3), thirdAlignment) #maps to second in ABC|ABC

                            else:                                                   #???????????????
                                addToMapInbetween(int(nt[2]/3), thirdAlignment)     #???????????????
                            inbetween = None
                        else:                                                       #simple codon MMM or DMIM or MDMM
                            addToMap(int(nt[2]/3), thirdAlignment)                  #maps to ABC
                elif (prevAaShift % 3) == 0:
                    if (nt[4] % 3) == 0: #3Mand3xI
                        if inbetween != None: #inbetween == "-"                                                     #both codon in MMI|DDDIIM and '-'
                            aa_alignment_map.update({int(nt[2]/3):(inbetween, firstAlignment, thirdAlignment)})     #maps to last codon in ABC|ABC 
                            inbetween = None
                        else:                                                                                       #codon such as IMM|IIM
                            aa_alignment_map.update({int(nt[2]/3):(firstAlignment, thirdAlignment)})                #maps to ABC
                    else: #1/2Iand3M
                        addToMapInbetween(int(nt[2]/3), firstAlignment)                                             #codon such as MII|MM...
                        inbetween = None                                                                            #maps to ABC
                elif (nt[4] % 3) == 0: #2D/1I...1M1D1M
                    if inbetween != None:
                        if hasDeletion:
                            if inbetween != '-':                                    #last codon in MIM|MMM|MMDM
                                inbetween = thirdAlignment                          #maps to last two in ABC|ABC|ABC
                                addToMapDeletion(int(nt[2]/3), thirdAlignment)
                            else:                                                   #last codon in MIM|MMDDDDM and '-'
                                addToMapInbetween(int(nt[2]/3), thirdAlignment)     #maps to last codon in ABC|ABC
                        else:
                            addToMapInbetween(int(nt[2]/3), thirdAlignment)         #codon MDDDMM and '-'
                        inbetween = None                                            #maps to last codon in ABC|ABC

                    else:                                                           #???????????????
                        addToMap(int(nt[2]/3), thirdAlignment)                      #???????????????
                prevAaShift = None
                firstAlignment = None

        aaQueryIndex = int(ntQueryIndex / 3)
        thirdAlignment = aaQueryIndex
        mapIndex += 1
    return aa_alignment_map    

