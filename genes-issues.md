# Disclaimer
This file is a work in progress

# Genes Issues - Kargva/CARD
## Problems with reference sequence
<details>
  <summary><h3>MEG_3586</h3></summary>
  Reference sequence in CARD is a DNA-binding response regulator in LuxR family, not liaR
</details>
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
  Single mutations listed in CARD that are actually part of an <em>n</em>-tuple mutation:
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
  Single mutations listed in CARD that are actually part of an <em>n</em>-tuple mutation:
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
  Single mutations listed in CARD that are actually part of an <em>n</em>-tuple mutation:
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
    <li>S244T - actually comes from a different gene</li>
  </ul>
  Mutations missing from CARD:
  <ul>
    <li>A439T</li>
    <li>H1002R</li>
    <li>V282G</li>
    <li>F285L</li>
    <li>D328H</li>
  </ul>
  <em>N</em>-tuple mutations that actually include at least one single mutation that can confer resistance by itself:
  <ul>
    <li>A313G,Y319C - Y319C is already listed as a resistance-confering single mutation</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_2711<sup><a href="#footnote2" id="reference2">2</a></sup></h3></summary> 
  Mutations that are suscepible or neutral:
  <ul>
    <li>Q998R</li>
    <li>T610K</li>
    <li>F1012S</li>
  </ul>
</details>


<details>
  <summary><h3>MEG_2712</h3></summary>
  Mutations listed in CARD that are not present in articles referenced: 
    <ul>
    <li>V287F</li>
    <li>H285Y</li>
    <li>A247P</li>
  </ul>
  Mutations that are susceptible or neutral:
  <ul>
    <li>T270I</li>
  </ul>
  Single mutations listed in CARD that are actually part of an <em>n</em>-tuple mutation:
  <ul>
    <li>G288V - part of triple mutation G288V, M310K, Y327N </li>
    <li>G288W - part of double mutation G288W, V303G </li>
    <li>Y296H - part of quadruple mutation T270I, Y296H, G308D, G325S </li>
    <li>Y296S - part of double mutation Y296S, R302G</li>
    <li>M300R - part of quadruple mutation G272S, H285Y, M300R, A307T</li>
    <li>R302G - part of double mutation Y296S, R302G</li>
    <li>V303G - part of double mutation G288W, V303G </li>
    <li>A307T - part of quadruple mutation G272S, H285Y, M300R, A307T</li>
    <li>G308D - part of quadruple mutation T270I, Y296H, G308D, G325S </li>
    <li>Y309N - part of double mutation V287F, Y309N </li>
    <li>M310K - part of triple mutation G288V, M310K, Y327N </li>
    <li>G325S - part of quadruple mutation T270I, Y296H, G308D, G325S </li>
    <li>W326R - part of double mutation I297L, W326R</li>
    <li>Y327N - part of triple mutation G288V, M310K, Y327N </li>
  </ul>
  Mutations missing from CARD:
  <ul>
    <li>E305D</li>
    <li>I406V</li>
  </ul>
  <em>N</em>-tuple mutations that actually include at least one single mutation that can confer resistance by itself:
  <ul>
    <li>A244T, G288W, V303G - not a single mutation, but G288W, V303G is already listd as a resistance-conferring double mutation</li>
    <li>A247P, T270I, I297T - I297T is already listed as a resistance-conferring single mutation</li>
    <li>A247P, I297L, W326R - I297L is already listed as a resistance-conferring single mutation</li>
    <li>L251R, A254G, T270I  - both L251R and A254G are individually listed as resistance-conferring single mutations</li>
    <li>T270I, G288W, V303G - not a single mutation, but G288W, V303G is already listd as a resistance-conferring double mutation</li>
    <li>T270I, I297T - I297T is already listed as a resistance-conferring single mutation</li>
    <li>I297L, W326R - I297L is already listed as a resistance-conferring single mutation</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_2866</h3></summary>
  Mutations listed in CARD that are not present in articles referenced: 
    <ul>
    <li>Deletion of nucleotide 65</li>
    <li>Insertion at nucleoide 811</li>
  </ul>
  Mutations missing from CARD:
  <ul>
    <li>Q246STOP</li>
    <li>GC deleted at nucleotides 1322-1323</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_2934<sup><a href="#footnote3" id="reference3">3</a></sup></h3></summary>
  Mutations listed in CARD that are not present in articles referenced: 
    <ul>
    <li>F25I</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3065</h3></summary>
  Mutations removed due to conflicting evidence:
  <ul>
    <li>F441Y - Position 441 had neither an F or a Y</li>
    <li>T387I, E449K- Neither position had their respsective wild type or mutant</li>
  </ul>
  <em>N</em>-tuple mutations that actually include at least one single mutation that can confer resistance by itself:
  <ul>
    <li>A67T, P406L - P406L is already listed as a resistance-conferring single mutation</li>
    <li>A70V, A160V, H457Y - H457Y is already listed as a resistance-conferring single mutation</li>
    <li>V90I, H457Q, L461K, A655V - V90I, H457Q and L461K are all already listed as individual resistance-conferring single mutations</li>
    <li>G452C, R659L - G452C is already listed as a resistance-conferring single mutation</li>
    <li>H457Y, S416F - H457Y is already listed as a resistance-conferring single mutation</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3135</h3></summary>
  Mutations that are susceptible or neutral:
  <ul>
    <li>P84L</li>
    <li>G37V</li>
    <li>E92A</li>
  </ul>
  Single mutations listed in CARD that are actually part of an <em>n</em>-tuple mutation:
  <ul>
    <li>W45C - part of triple mutation W45C, W148R</li>
    <li>W148R - part of double mutation W45C, W148R</li>
  </ul>
  Mutations listed wrong in CARD:
  <ul>
    <li>N52T - is actually a N51T mutation</li>
    <li>L49F - is actually a L49P mutation</li>
  </ul>
  Mutations missing from CARD:
  <ul>
    <li>V65A</li>
    <li>E92D</li>
    <li>L128S</li>
    <li>D132V</li>
    <li>G insertion at nt 352</li>
    <li>GCGCCGAGGAG deletion from nt 353</li>
    <li>G deletion at nt 102</li>
    <li>G deletion at nt 351</li>
    <li>C deletion at nt 115</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3143</h3></summary>
  Mutations that are susceptible or neutral:
  <ul>
    <li>F3I</li>
    <li>L27F</li>
    <li>A100V</li>
    <li>V213I</li>
    <li>G352D</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3144</h3></summary>
  Mutations missing from CARD:
  <ul>
    <li>insertion of ISEcp1 (1661 bp) between 766T and 767G</li>
    <li>insertion of ISEcp1 (1880 bp) between 552T and 553G</li>
    <li>20 bp deletion from nucleotide 602 to 621</li>
  </ul>
  <em>N</em>-tuple mutations that actually include at least one single mutation that can confer resistance by itself:
  <ul>
    <li>E448K, Q444E, E443Q, L297F - E448K is already listed as a resistance-conferring single mutation</li>
    <li>E448K, G302D - E448K is already listed as a resistance-conferring single mutation</li>
    <li>E448K, G33R - E448K is already listed as a resistance-conferring single mutation</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3177</h3></summary>
  Single mutations listed in CARD that are actually part of an <em>n</em>-tuple mutation:
  <ul>
    <li>D95N - part of triple mutation S91F, D95N</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3178<sup><a href="#footnote4" id="reference4">4</a></sup></h3></summary>
  Mutations listed in CARD that are not present in articles referenced:
  <ul>
    <li>S83L</li>
  </ul>
  Single mutations listed in CARD that are actually part of an <em>n</em>-tuple mutation:
  <ul>
    <li>N57K - part of double mutation N57K, S83L</li>
    <li>H80P - part of double mutation H80P, S83L</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3180<sup><a href="#footnote5" id="reference5">5</a></sup></h3></summary>
  Mutations listed in CARD that are not present in articles referenced:
  <ul>
    <li>G88A</li>
    <li>D89V</li>
    <li>D94V</li>
  </ul>
  Mutations that are susceptible or neutral:
  <ul>
    <li>T80A</li>
    <li>A90V - by itself</li>
    <li>A90G</li>
    <li>G247S - by itself</li>
  </ul>
  Single mutations listed in CARD that are actually part of an <em>n</em>-tuple mutation:
  <ul>
    <li>P102H - part of double mutation A90V, P102H</li>
  </ul>
  Mutations missing from CARD:
  <ul>
    <li>T80A, A90G - is hypersusceptible<sup><a href="#footnote6" id="reference6">6</a></sup></li>
  </ul>
  <em>N</em>-tuple mutations that actually include at least one single mutation that can confer resistance by itself:
  <ul>
    <li>A90V, D94G - D94G is already listed as a resistance-conferring single mutation</li>
    <li>MEG_3180 A90V, MEG_3243 D472H- MEG_3243 D472H is already listed as a resistance-conferring single mutation</li>
    <li>MEG_3180 G247S, MEG_3243 D500N - MEG_3243 D500N is already listed as a resistance-conferring single mutation</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3181</h3></summary>
  Mutations listed wrong in CARD:
  <ul>
    <li>A91T - is actually a A91V mutation</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3182</h3></summary>
  Mutations missing from CARD:
  <ul>
    <li>S101W</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3185</h3></summary>
  Mutations missing from CARD:
  <ul>
    <li>D99N</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3187</h3></summary>
  Mutations listed in CARD that are not present in articles referenced:
  <ul>
    <li>D87Y</li>
    <li>D87H</li>
  </ul>
  Single mutations listed in CARD that are actually part of an <em>n</em>-tuple mutation:
  <ul>
    <li>S83I - part of quintuple mutation S83I, I112V, L127M, A128S, K154R</li>
  </ul>
  <em>N</em>-tuple mutations that actually include at least one single mutation that can confer resistance by itself:
  <ul>
    <li>S83L, D87Y - S83L is already listed as a resistance-conferring single mutation</li>
    <li>S83L, D87N - S83L and D87N are already listed as individual resistance-conferring single mutations</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3240</h3></summary>
  Mutations missing from CARD:
  <ul>
    <li>R523Q</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3241<sup><a href="#footnote7" id="reference7">7</a></sup></h3></summary>
  Mutations that are susceptible or neutral:
  <ul>
    <li>S416A</li>
  </ul>
  Mutations listed wrong in CARD:
  <ul>
    <li>K444F - is actually a L444F mutation</li>
  </ul>
  Mutations listed wrong in KARGVA:
  <ul>
    <li>E467V - is actually a E466V mutation</li>
    <li>E467K - is actually a E466K mutation</li>
    <li>R388P - is actually a R389P mutation</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3242</h3></summary>
  Mutations listed wrong in KARGVA:
  <ul>
    <li>S464A - is actually a S463A mutation</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3243<sup><a href="#footnote8" id="reference8">8</a></sup></h3></summary>
  Mutations that are susceptible or neutral:
  <ul>
    <li>V340L</li>
  </ul>
  Mutations listed wrong in CARD:
  <ul>
    <li>N510D - is actually a N538D mutation</li>
    <li>D472H  - is actually a D500H mutation</li>
  </ul>
  <em>N</em>-tuple mutations that actually include at least one single mutation that can confer resistance by itself:
  <ul>
    <li>MEG_3180 G247S, MEG_3243 D500N - MEG_3243 D500N is already listed as a resistance-conferring single mutation</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3244</h3></summary>
  Mutations listed in CARD that are not present in articles referenced:
  <ul>
    <li>R136I</li>
    <li>R136E</li>
    <li>R136G</li>
    <li>R136L</li>
  </ul>
  Mutations missing from CARD:
  <ul>
    <li>G164V</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3245<sup><a href="#footnote9" id="reference9">9</a></sup></h3></summary>
  Mutations that are susceptible or neutral:
  <ul>
    <li>S128L</li>
  </ul>
  Single mutations listed in CARD that are actually part of an <em>n</em>-tuple mutation:
  <ul>
    <li>I56S - part of double mutation I56S, R144S</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3246<sup><a href="#footnote10" id="reference10">10</a></sup></h3></summary>
  Mutations listed wrong in KARGVA:
  <ul>
    <li>G125S - is actually a G124S mutation</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3429</h3></summary>
  Mutations listed in CARD that are not present in articles referenced:
  <ul>
    <li>I21T</li>
  </ul>
  Mutations that are susceptible or neutral:
  <ul>
    <li>G3G</li>
  </ul>
  Mutations listed wrong in CARD:
  <ul>
    <li>V78A - is actually a S94A mutation</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3430</h3></summary>
  Mutations missing from CARD:
  <ul>
    <li>5 bp deletion from nt 282 to 286</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3445</h3></summary>
  Mutations that are susceptible or neutral:
  <ul>
    <li>G312S</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3446</h3></summary>
  Mutations removed due to conflicting evidence:
  <ul>
    <li>A234G - G is the wild-type for this organism
      <ul>
        <li>When looking at the catalase/peroxidase domain in other organisms, G is the wild-type as well</li>
        <li>Additionally, all instances of A234G mutations were accompanied with nonsense mutations</li>
      </ul>
    </li>
    <li>A431V - V is the wild-type for this organism
      <ul>
        <li>CARD also has V431A as a SNP</li>
        <li>When looking at the catalase/peroxidase domain in other organisms; V is the wild-type as well</li>
        <li>However, A431V is a single mutation, while V431A is in a double mutation with G490S
          <ul>
            <li>G490S is also a single mutation</li>
          </ul>
        </li>
      </ul>
    </li>
  </ul>
  Mutations listed in CARD that are not present in articles referenced:
  <ul>
    <li>G593D</li>
  </ul>
  Mutations that are susceptible or neutral:
  <ul>
    <li>R463L</li>
  </ul>
  Single mutations listed in CARD that are actually part of an <em>n</em>-tuple mutation:
  <ul>
    <li>S17T - part of triple mutation S17T, T insertion at nt 277<sup><a href="#footnote11" id="reference11">11</a></sup>, and S315T</li>
    <li>V68G - part of triple mutation V68G, A234G, E454*<sup><a href="#footnote11" id="reference11">11</a></sup></li>
    <li>Y98C - part of triple mutation Y98C, A349T, R463L</li>
    <li>G123E - part of double mutation G123E, G299S</li>
    <li>P131R - part of double mutation P131R, D194Y </li>
    <li>D194Y - part of double mutation P131R, D194Y </li>
    <li>D194G - part of triple mutation D194G, S315T, M624V</li>
    <li>G299S - part of double mutation G123E, G299S</li>
    <li>A349T - part of triple mutation Y98C, A349T, R463L</li>
    <li>D387H - part of double mutation S315T, D387H</li>
    <li>M624V - part of triple mutation D194G, S315T, M624V</li>
  </ul>
  Mutations listed wrong in CARD:
  <ul>
    <li>D36E - is actually a D63E mutation</li>
    <li>G300W - is actually a W300G mutation</li>
    <li>A717P - is actually a A716P mutation</li>
  </ul>
  Mutations missing from CARD:
  <ul>
    <li>K46*<sup><a href="#footnote11" id="reference11">11</a></sup>, A234G</li>
    <li>R128P</li>
    <li>N138H</li>
    <li>N138T</li>
    <li>L148R</li>
    <li>A172T</li>
    <li>E195K</li>
    <li>W198*</li>
    <li>W204R</li>
    <li>E217G</li>
    <li>A264V</li>
    <li>H270Q</li>
    <li>G285C</li>
    <li>Y304S</li>
    <li>T308P</li>
    <li>W321G</li>
    <li>L336P</li>
    <li>W341S</li>
    <li>T344P</li>
    <li>D381G</li>
    <li>K414N</li>
    <li>K433Q, Q434A, T435D</li>
    <li>V473*</li>
    <li>G490S - by itself</li>
    <li>G491D</li>
    <li>G491V</li>
    <li>W505S</li>
    <li>R515C</li>
    <li>D735A</li>
    <li>Deletion of amino acids W191 and E192 </li>
    <li>Deletion of nucleotides 478 and 479</li>
    <li>Deletion of nucleotide 371</li>
    <li>Single bp insertion in codon 304</li>
    <li>Single bp insertion in codon 521</li>
    <li>8bp deletion from codon 217 to 219</li>
    <li>22bp deletion from codon 289 to codon 296</li>
    <li>GGTC insertion in codon 518</li>
    <li>C deletion at nt 6</li>
    <li>A deletion at nt 12</li>
    <li>C deletion at nt 33</li>
    <li>G deletion at nt 78</li>
    <li>G deletion at nt 193</li>
    <li>G deletion at nt 201</li>
    <li>G deletion at nt 320</li>
    <li>T insertion at nt 344</li>
    <li>G deletion at nt 1245</li>
    <li>G insertion at nt 1317</li>
    <li>A deletion at nt 1384</li>
    <li>G deletion at nt 1453</li>
    <li>G deletion at nt 1472</li>
    <li>A deletion at nt 1477</li>
  </ul>
  <em>N</em>-tuple mutations that actually include at least one single mutation that can confer resistance by itself:
  <ul>
    <li>Y98C, A349T, R463L - Not a single mutation, but because R463L is a neutral mutation, only Y98C and A349T confer resistance here
    <li>D194G, S315T, M624V - S315T is already listed as a resistance-conferring single mutation</li>
    <li>S315T, D387H - S315T is already listed as a resistance-conferring single mutation</li>
    <li>V431A, G490S - G490S is already listed as a resistance-conferring single mutation</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3588</h3></summary>
  Mutations listed wrong in CARD:
  <ul>
	  <li>I177 insertion and deletion - actually located between positions 168 and 171
      <ul>
        <li>Information is based on the fact that Ile repeat is located at those positions</li>
      </ul>
    </li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3589</h3></summary>
  Mutations listed in CARD that are not present in articles referenced:
  <ul>
	  <li>A180T</li>
  </ul>
  Single mutations listed in CARD that are actually part of an <em>n</em>-tuple mutation:
  <ul>
	  <li>H264Q - part of double mutation T120A, H264Q</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3591</h3></summary>
  Mutations missing from CARD:
  <ul>
	  <li>L27I</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3594</h3></summary>
  Single mutations listed in CARD that are actually part of an <em>n</em>-tuple mutation:
  <ul>
	  <li>Q52P - part of double mutation Q52P, Ter189S</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3626</h3></summary>
  Mutations missing from CARD:
  <ul>
    <li>Q234*</li>
    <li>2-bp deletion at nt 76</li>
    <li>445-bp deletion at nt 364</li>
    <li>30-bp deletion at nt 391</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3627</h3></summary>
  Mutations missing from CARD:
  <ul>
    <li>Single base deletion at nt 135</li>
    <li>84-bp deletion at nt 858</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_3931</h3></summary>
  Mutations that are susceptible or neutral:
  <ul>
    <li>V104A</li>
  </ul>
  Single mutations listed in CARD that are actually part of an <em>n</em>-tuple mutation:
  <ul>
    <li>V73A - part of double mutation V73A, L270Q</li>
    <li>L270Q - part of double mutation V73A, L270Q</li>
  </ul>
  Mutations missing from CARD:
  <ul>
    <li>First 809bp deletion</li>
    <li>C deletion at nt 293</li>
    <li>8-bp deletion at nt 710</li>
    <li>30-bp deletion at nt 927</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_4057<sup><a href="#footnote12" id="reference12">12</a></sup></h3></summary>
  Mutations removed due to wild-type not being in sequence:
  <ul>
    <li>V351E</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_4092</h3></summary>
  Mutations listed wrong in CARD:
  <ul>
    <li>D116C - is actually intrinsically resistant at D116</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_4094</h3></summary>
  Mutations that are susceptible or neutral:
  <ul>
    <li>V65L</li>
    <li>G257D</li>
    <li>D278E</li>
    <li>E291D</li>
    <li>Q362R</li>
    <li>T396N</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_4095</h3></summary>
  Mutations listed in CARD that are not present in articles referenced:
  <ul>
    <li>C115S</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_4109</h3></summary>
  Single mutations listed in CARD that are actually part of an <em>n</em>-tuple mutation:
  <ul>
    <li>A186T - part of double mutation G71E, A186T</li>
    <li>S209R - part of double mutation G71E, S209R</li>
  </ul>
  Mutations missing from CARD:
  <ul>
    <li>E72K, D75H, A85D</li>
    <li>T100K, L110P</li>
    <li>P105A, P122L, S141N</li>
    <li>E204G, H205N</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_4110</h3></summary>
  Mutations missing from CARD:
  <ul>
    <li>T158I</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_4130</h3></summary>
  Mutations that are susceptible or neutral:
  <ul>
    <li>V18A</li>
  </ul>
  Mutations missing from CARD:
  <ul>
    <li>A226E</li>
    <li>N316K</li>
    <li>G339A</li>
    <li>G insertion at nt 969</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_4131</h3></summary>
  Mutations missing from CARD:
  <ul>
    <li>A insertion at nt 272</li>
    <li>T insertion at nt 439</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_4289</h3></summary>
  Mutations that are susceptible or neutral:
  <ul>
    <li>R154A</li>
    <li>R154D</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_4296</h3></summary>
  Mutations listed wrong in CARD:
  <ul>
    <li>Q142X - is actually a Q142* mutation</li>
  </ul>
  Mutations missing from CARD:
  <ul>
    <li>W339*</li>
    <li>14-bp deletion at nt 433</li>
    <li>A insertion at nt 457</li>
    <li>17-bp deletion at nt 482</li>
    <li>T insertion at nt 673</li>
    <li>2-bp AT insertion at nt 1113</li>
    <li>2-bp deletion at nt 1163</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_5325</h3></summary>
  Mutations listed wrong in CARD:
  <ul>
    <li>S80L - is actually a S87L mutation</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_5328</h3></summary>
  Mutations missing from CARD:
  <ul>
    <li>D82N</li>
    <li>E87K</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_5331</h3></summary>
  Mutations that are susceptible or neutral:
  <ul>
    <li>P62S</li>
    <li>A67S</li>
    <li>A69T</li>
    <li>S83N</li>
    <li>A119V</li>
  </ul>
  Mutations listed wrong in CARD:
  <ul>
    <li>A66T - is actually a A69T mutation, which is susceptible</li>
  </ul>
  Mutations missing from CARD:
  <ul>
    <li>D87V</li>
    <li>K97R</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_5332</h3></summary>
  Single mutations listed in CARD that are actually part of an <em>n</em>-tuple mutation:
  <ul>
    <li>A85T - part of double mutation with an S80 mutation</li>
    <li>D111H - part of double mutation with an S80 mutation</li>
    <li>S129P - part of double mutation with an S80 mutation</li>
  </ul>
  Mutations missing from CARD:
  <ul>
    <li>S80I</li>
    <li>S80R</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_5333</h3></summary>
  Single mutations listed in CARD that are actually part of an <em>n</em>-tuple mutation:
  <ul>
    <li>E84G - part of double mutation S80I, E84G</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_5394</h3></summary>
  Mutations missing from CARD:
  <ul>
    <li>L416F</li>
    <li>I444F</li>
    <li>L445H</li>
    <li>S458T</li>
    <li>E460D</li>
    <li>I464F</li>
    <li>I529L</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_5401</h3></summary>
  Mutations missing from CARD:
  <ul>
    <li>V78I</li>
    <li>H168Y</li>
    <li>D346 insertion</li>
    <li>A501V</li>
    <li>InsertHere</li>
    <li>InsertHere</li>
  </ul>
</details>

---

<sup id="footnote1">
  1. Although the PMID is in progress under this particular entry, the original article was found in the MEG_6045 entry <a href="#reference1">↩</a>
</sup><br>
<sup id="footnote2">
  2. All mutations in CARD are not resistance-conferring; therefore the gene is not present in the SNPInfo database <a href="#reference2">↩</a>
</sup><br>
<sup id="footnote3">
  3. No mutant at 25 for Streptococcus pyogenes (this mutation exists in E.coli). Many mutations are listed in figure 4 of article referenced in CARD (PMID: 9593127), but it is unclear which ones actually confer resistance. Therefore, the gene is not present in the SNPInfo database <a href="#reference3">↩</a>
</sup><br>
<sup id="footnote4">
  4. Some of the mutant strains are susceptible to norfloxacin but resistant to nalidixic acid <a href="#reference4">↩</a>
</sup><br>
<sup id="footnote5">
  5. The CARD entry features a double loci mutation with MEG_3243; some of the information listed under MEG_3180 came from one of MEG_3243's referenced article: PMID: 22761889 <a href="#reference5">↩</a>
</sup><br>
<sup id="footnote6">
  6. More explanations on hypersusceptible mutations <a href="https://github.com/Isabella136/AmrPlusPlus_SNP/blob/main/data/SNPInfoGuide.md#hyper-why-nonresistant-snps-are-included">here </a><a href="#reference6">↩</a>
</sup><br>
<sup id="footnote7">
  7. The organism used to determine the mutations was not the same as the one in CARD/MEGARes, so a pairwise alignment between the two sequences was made to figure out the SNPs' positions based on the CARD/MEGARes sequence. However, the original positions are used in this file when referring to mutations <a href="#reference7">↩</a>
</sup><br>
<sup id="footnote8">
  8. The sequence in CARD didn't match the sequence in MEGARes, so a pairwise alignment between the two sequences was made to figure out the SNPs' positions based on the CARD/MEGARes sequence. However, the original positions are used in this file when referring to mutations <a href="#reference8">↩</a>
</sup><br>
<sup id="footnote9">
  9. The SNPs included doesn’t exactly match the sequence in CARD; you need to subtract 1 from the SNPs’ position to find the corresponding positions for the sequence in CARD <a href="#reference9">↩</a>
</sup><br>
<sup id="footnote10">
  10. The sequence used to find the SNPs is an older version of the one in CARD/MEGARes, so a pairwise alignment between the two sequences was made to figure out the SNPs' positions based on the CARD/MEGARes sequence. However, the original positions are used in this file when referring to mutations <a href="#reference10">↩</a>
</sup><br>
<sup id="footnote11">
  11. In an <em>n</em>-tuple mutation, nonsense and frameshift mutations take precedence and are therefore considered single mutations <a href="#reference11">↩</a>
</sup><br>
<sup id="footnote12">
  12. The CARD sequence contains a 4-aa deletion from positions 349 to 352. All SNPs listed in this entry (minus F57S and I418N) are based on the articles’ sequence which doesn’t contain the 4-aa deletion <a href="#reference12">↩</a>
</sup><br>
