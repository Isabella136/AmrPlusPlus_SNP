# Disclaimer
This file is a work in progress

# Genes Issues - Kargva/CARD
## Problems with reference sequence
<details>
  <summary><h3>MEG_3594</h3></summary>
  Reference sequence in CARD is from a different species than the reference sequence in article listed
  <ul>
    <li>CARD sequence is S. lincolnensis, which produced lincomycin; its lmrA is intrinsically resistant to lincomycin</li>
    <li>Article sequence is B. subtilis which can acquire resistance in its lmrA gene through SNP</li>
    <li>SNP information for MEG_3594 is still included in SNPInfo database but is never used by the SNP_Verification program</li>
    <ul>
      <li>Information is based on article sequence</li>
    </ul>
  </ul>
</details>
<details>
  <summary><h3>MEG_3994</h3></summary>
  Despite being an rRNA sequence, the reference nucleotide sequence retrieved from CARD contains amino acid residues
</details>
<details>
  <summary><h3>MEG_5779</h3></summary>
  Reference sequence in CARD is different from sequence in MEGARes
  <ul>
    <li>CARD sequence is a phosphatidyltransferase called pgsA</li>
    <li>MEGARes sequence is a capsular polyglutamate synthetase called capA which is also known as pgsA in NCBI</li>
    <li>SNP information for MEG_5779 is still included in SNPInfo database but is never used by the SNP_Verification program</li>
    <ul>
      <li>Info is based on CARD sequence</li>
    </ul>
  </ul>
</details>

## Problems with mutations
<details>
  <summary><h3>MEG_413</h3></summary>
  Mutations listed in CARD that are not present in articles referenced: 
  <ul>
    <li>R45C</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_414<sup><a href="#footnote1" id="reference1">1</a></sup></h3></summary> 
  Single mutations listed in CARD that are actually part of an *n*-tuple mutation:
  <ul>
    <li>Y114F - part of the double mutation Y114F, V165I</li>
    <li>V165I - part of the double mutation Y114F, V165I</li>
  </ul>
  Mutations missing from CARD:
  <ul>
    <li>Transposase insertion</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_1187</h3></summary>
  Mutations missing from CARD:
  <ul>
    <li>R105S</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_1731</h3></summary>
  Single mutations listed in CARD that are actually part of an *n*-tuple mutation:
  <ul>
    <li>+AII14-16 - part of double mutation +AII14-16, -NFQ74-76</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_1732</h3></summary>
  Mutations listed wrong in CARD:
  <ul>
    <li>K59T - is actually a K59N mutation</li>
  </ul>
  Mutations missing from CARD:
  <ul>
    <li>N13S</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_2710</h3></summary>
  Mutations that are susceptible or neutral:
  <ul>
    <li>L413P</li>
    <li>E504Q</li>
    <li>D1024N</li>
  </ul>
  Single mutations listed in CARD that are actually part of an *n*-tuple mutation:
  <ul>
    <li>R507G - part of double mutation M306I, R507G</li>
    <li>R471P - part of double mutation D299E, R471P</li>
    <li>R469P - part of double mutation D299E, R469P</li>
    <li>I465D - part of double mutation D299E, I465D</li>
    <li>P446H - part of double mutation D299E, P446H</li>
    <li>P397Q - part of quadruple mutation M306V, E368A, S380R, P397Q</li>
    <li>S380R - part of quadruple mutation M306V, E368A, S380R, P397Q</li>
    <li>E378A - part of double mutation D299E, E378A</li>
    <li>E368A - part of quadruple mutation M306V, E368A, S380R, P397Q</li>
  </ul>
  Mutations listed wrong in CARD:
  <ul>
    <li>A314G,Y322C - is actually a A313G,Y319C double mutation</li>
  </ul>
  Mutations missing from CARD:
  <ul>
    <li>N13S</li>
  </ul>
  Mutations listed wrong in CARD:
  <ul>
    <li>K59T - is actually a K59N mutation</li>
  </ul>
  Mutations missing from CARD:
  <ul>
    <li>N13S</li>
  </ul>
</details>


Mutations missing from CARD:
- A439T
- H1002R
- V282G
- F285L
- D328H

*N*-tuple mutations that actually include at least one single mutation that can confer resistance by itself:
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

Single mutations listed in CARD that are actually part of an *n*-tuple mutation:
- G288V - part of triple mutation G288V, M310K, Y327N 
- G288W - part of double mutation G288W, V303G 
- Y296H - part of quadruple mutation T270I, Y296H, G308D, G325S 
- Y296S - part of double mutation Y296S, R302G
- M300R - part of quadruple mutation G272S, H285Y, M300R, A307T
- R302G - part of double mutation Y296S, R302G
- V303G - part of double mutation G288W, V303G 
- A307T - part of quadruple mutation G272S, H285Y, M300R, A307T
- G308D - part of quadruple mutation T270I, Y296H, G308D, G325S 
- Y309N - part of double mutation V287F, Y309N 
- M310K - part of triple mutation G288V, M310K, Y327N 
- G325S - part of quadruple mutation T270I, Y296H, G308D, G325S 
- W326R - part of double mutation I297L, W326R
- Y327N - part of triple mutation G288V, M310K, Y327N 

Mutations missing from CARD:
- E305D 
- I406V

*N*-tuple mutations that actually include at least one single mutation that can confer resistance by itself:
- A244T, G288W, V303G - not a single mutation, but G288W, V303G is already listd as a resistance-conferring double mutation
- A247P, T270I, I297T - I297T is already listed as a resistance-conferring single mutation
- A247P, I297L, W326R - I297L is already listed as a resistance-conferring single mutation
- T270I, I297T - I297T is already listed as a resistance-conferring single mutation
- I297L, W326R - I297L is already listed as a resistance-conferring single mutation

### MEG_2866

---

<sup id="footnote1">
  1.Although the PMID is in progress under this particular entry, the original article was found in the MEG_6045 entry <a href="#reference1">↩</a>
</sup>

[^2]: All mutations in CARD are not resistance-conferring; therefore the gene is not present in the SNPInfo database
