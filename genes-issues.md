# Disclaimer
This file is a work in progress

# Genes Issues - Kargva/CARD
## Problems with reference sequence
### MEG_3594
Reference sequence in CARD is from a different species than the reference sequence in article listed
- CARD sequence is S. lincolnensis, which produced lincomycin; its lmrA is intrinsically resistant to lincomycin
- Article sequence is B. subtilis which can acquire resistance in its lmrA gene through SNP 
- SNP information for MEG_3594 is still included in SNPInfo database but is never used by the SNP_Verification program
  - Information is based on article sequence
### MEG_3994
Despite being an rRNA sequence, the reference nucleotide sequence retrieved from CARD contains amino acid residues
### MEG_5779
Reference sequence in CARD is different from sequence in MEGARes
- CARD sequence is a phosphatidyltransferase called pgsA
- MEGARes sequence is a capsular polyglutamate synthetase called capA which is also known as pgsA in NCBI
- SNP information for MEG_5779 is still included in SNPInfo database but is never used by the SNP_Verification program
  - Info is based on CARD sequence

## Problems with mutations
### MEG_413
Mutations listed in CARD that are not present in articles referenced: 
- R45C

### MEG_414 [^1]
Single mutations listed in CARD that are actually part of an *n*-tuple mutation:
- Y114F - part of the double mutation Y114F, V165I
- V165I - part of the double mutation Y114F, V165I

Mutations missing from CARD:
- Transposase insertion

### MEG_1187
Mutations missing from CARD:
- R105S

### MEG_1731
Single mutations listed in CARD that are actually part of an *n*-tuple mutation:
- +AII14-16 - part of double mutation +AII14-16, -NFQ74-76

### MEG_1732
Mutations listed wrong in CARD:
- K59T - is actually a K59N mutation

Mutations missing from CARD:
- N13S

### MEG_2710
Mutations that are susceptible or neutral:
- L413P
- E504Q
- D1024N

Single mutations listed in CARD that are actually part of an *n*-tuple mutation:
- R507G - part of double mutation M306I, R507G
- R471P - part of double mutation D299E, R471P
- R469P - part of double mutation D299E, R469P
- I465D - part of double mutation D299E, I465D
- P446H - part of double mutation D299E, P446H
- P397Q - part of quadruple mutation M306V, E368A, S380R, P397Q
- S380R - part of quadruple mutation M306V, E368A, S380R, P397Q
- E378A - part of double mutation D299E, E378A
- E368A - part of quadruple mutation M306V, E368A, S380R, P397Q

Mutations listed wrong in CARD:
- A314G,Y322C - is actually a A313G,Y319C double mutation

Mutations missing from CARD:
- A439T
- H1002R
- V282G
- F285L
- D328H

*N*-tuple mutations that actually include at least one single mutation that can confer resistance by itself
- A313G,Y319C - Y319C is already listed as a resistance-confering single mutation

### MEG_2711
Mutations that are suscepible or neutral: [^2]
- Q998R
- T610K
- F1012S

### MEG_2712
Mutations listed in CARD that are not present in articles referenced: 
- V287F
- H285Y
- A247P

Mutations that are susceptible or neutral:
- T270I

Mutations missing from CARD:
- E305D 
- I406V

*N*-tuple mutations that actually include at least one single mutation that can confer resistance by itself
- A244T, G288W, V303G - not a single mutation, but G288W, V303G is already listd as a resistance-conferring double mutation
- A247P, T270I, I297T - I297T is already listed as a resistance-conferring single mutation
- A247P, I297L, W326R - I297L is already listed as a resistance-conferring single mutation
- T270I, I297T - I297T is already listed as a resistance-conferring single mutation
- T270I, I297L, W326R - I297L is already listed as a resistance-conferring single mutation
- I297L, W326R - I297L is already listed as a resistance-conferring single mutation

[^1]: Although the PMID is in progress under this particular entry, the original article was found in the MEG_6045 entry
[^2]: All mutations in CARD are not resistance-conferring; therefore the gene is not present in the SNPInfo database
