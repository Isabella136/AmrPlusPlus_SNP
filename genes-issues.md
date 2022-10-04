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
  Mutations listed wrong in CARD:
  <ul>
    <li>A314G,Y322C - is actually a A313G,Y319C double mutation</li>
    <li>S244T - actually comes from a different gene</li>
  </ul>
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
  <em>N</em>-tuple mutations that actually include at least one single mutation that can confer resistance by itself:
  <ul>
    <li>D299E, E378A - D299E is already listed as a resistance-confering single mutation</li>
    <li>D299E, P446H - D299E is already listed as a resistance-confering single mutation</li>
    <li>D299E, I465D - D299E is already listed as a resistance-confering single mutation</li>
    <li>D299E, I465D - D299E is already listed as a resistance-confering single mutation</li>
    <li>D299E, R469P - D299E is already listed as a resistance-confering single mutation</li>
    <li>D299E, R471P - D299E is already listed as a resistance-confering single mutation</li>
    <li>M306I, R507G - M306I is already listed as a resistance-confering single mutation</li>
    <li>M306V, E368A, S380R, P397Q - M306V is already listed as a resistance-confering single mutation</li>
    <li>A313G, Y319C - Y319C is already listed as a resistance-confering single mutation</li>
  </ul>
  Mutations missing from CARD:
  <ul>
    <li>A439T</li>
    <li>H1002R</li>
    <li>V282G</li>
    <li>F285L</li>
    <li>D328H</li>
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
  Mutations missing from CARD:
  <ul>
    <li>E305D</li>
    <li>I406V</li>
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
  Mutations listed wrong in CARD:
  <ul>
    <li>N52T - is actually a N51T mutation</li>
    <li>L49F - is actually a L49P mutation</li>
  </ul>
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
  <em>N</em>-tuple mutations that actually include at least one single mutation that can confer resistance by itself:
  <ul>
    <li>S91F, D95N - S91F is already listed as a resistance-conferring single mutation</li>
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
  <em>N</em>-tuple mutations that actually include at least one single mutation that can confer resistance by itself:
  <ul>
    <li>N57K, S83L - S83L is already listed as a resistance-conferring single mutation</li>
    <li>H80P, S83L - S83L is already listed as a resistance-conferring single mutation</li>
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
  <em>N</em>-tuple mutations that actually include at least one single mutation that can confer resistance by itself:
  <ul>
    <li>A90V, D94G - D94G is already listed as a resistance-conferring single mutation</li>
    <li>MEG_3180 A90V, MEG_3243 D472H- MEG_3243 D472H is already listed as a resistance-conferring single mutation</li>
    <li>MEG_3180 G247S, MEG_3243 D500N - MEG_3243 D500N is already listed as a resistance-conferring single mutation</li>
  </ul>
  Mutations missing from CARD:
  <ul>
    <li>T80A, A90G - is hypersusceptible<sup><a href="#footnote6" id="reference6">6</a></sup></li>
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
  Mutations that are susceptible or neutral:
  <ul>
    <li>S416A</li>
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
  Mutations listed wrong in CARD:
  <ul>
    <li>N510D - is actually a N538D mutation</li>
    <li>D472H  - is actually a D500H mutation</li>
  </ul>
  Mutations that are susceptible or neutral:
  <ul>
    <li>V340L</li>
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
  <em>N</em>-tuple mutations that actually include at least one single mutation that can confer resistance by itself:
  <ul>
    <li>I56S, R144S - R144S is already listed as a resistance-conferring single mutation</li>
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
  Mutations listed wrong in CARD:
  <ul>
    <li>V78A - is actually a S94A mutation</li>
  </ul>
  Mutations listed in CARD that are not present in articles referenced:
  <ul>
    <li>I21T</li>
  </ul>
  Mutations that are susceptible or neutral:
  <ul>
    <li>G3G</li>
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
  Mutations listed wrong in CARD:
  <ul>
    <li>D36E - is actually a D63E mutation</li>
    <li>G300W - is actually a W300G mutation</li>
    <li>A717P - is actually a A716P mutation</li>
  </ul>
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
  <em>N</em>-tuple mutations that actually include at least one single mutation that can confer resistance by itself:
  <ul>
    <li>Y98C, A349T, R463L - Not a single mutation, but because R463L is a neutral mutation, only Y98C and A349T confer resistance here
    <li>D194G, S315T, M624V - S315T is already listed as a resistance-conferring single mutation</li>
    <li>S315T, D387H - S315T is already listed as a resistance-conferring single mutation</li>
    <li>V431A, G490S - G490S is already listed as a resistance-conferring single mutation</li>
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
  <em>N</em>-tuple mutations that actually include at least one single mutation that can confer resistance by itself:
  <ul>
    <li>T120A, H264Q - T120A is already listed as a resistance-conferring single mutation</li>
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
  <em>N</em>-tuple mutations that actually include at least one single mutation that can confer resistance by itself:
  <ul>
    <li>G71E, A186T - G71E is already listed as a resistance-conferring single mutation</li>
    <li>G71E, S209R - G71E is already listed as a resistance-conferring single mutation</li>
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
  Mutations listed wrong in CARD:
  <ul>
    <li>A66T - is actually a A69T mutation, which is susceptible</li>
  </ul>
  Mutations that are susceptible or neutral:
  <ul>
    <li>P62S</li>
    <li>A67S</li>
    <li>A69T</li>
    <li>S83N</li>
    <li>A119V</li>
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
  <em>N</em>-tuple mutations that actually include at least one single mutation that can confer resistance by itself:
  <ul>
    <li>S80I - S80I, E84G is already listed as a resistance-conferring single mutation</li>
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
  Mutations listed in CARD that are not present in articles referenced:
  <ul>
    <li>N512Y</li>
  </ul>
	Single mutations listed in CARD that are actually part of an <em>n</em>-tuple mutation:
  <ul>
    <li>I312M  - part of double mutation I312M, G545S</li>
    <li>V316T - part of double mutation V316T, G545S</li>
  </ul>
  <em>N</em>-tuple mutations that actually include at least one single mutation that can confer resistance by itself:
  <ul>
    <li>G545S - I312M, G545S is already listed as a resistance-conferring single mutation</li>
    <li>G545S - V316T, G545S is already listed as a resistance-conferring single mutation</li>
  </ul>
  Mutations missing from CARD:
  <ul>
    <li>D346 insertion</li>
    <li>A501V</li>
    <li>A516G</li>
    <li>E538G</li>
    <li>P551S</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_5406</h3></summary>
  Mutations listed wrong in CARD:
  <ul>
    <li>T445A - is actually a T451A mutation</li>
    <li>E475G - is actually a E481G mutation</li>
    <li>T488A - is actually a T494A mutation</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_5407</h3></summary>
  Mutations that are susceptible or neutral:
  <ul>
    <li>M339F</li>
    <li>M400T</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_5780</h3></summary>
  Mutations listed wrong in CARD:
  <ul>
    <li>V59N - is actually a V59T mutation</li>
    <li>K66R - is actually a K75N mutation</li>
  </ul>
  Single mutations listed in CARD that are actually part of an <em>n</em>-tuple mutation:
  <ul>
    <li>K75N - part of double mutation K75N, GE76 insertion</li>
    <li>G76 insertion - part of double mutation K75N, GE76 insertion</li>
    <li>E77 insertion - part of double mutation K75N, GE76 insertion</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_5784<sup><a href="#footnote13" id="reference13">13</a></sup></h3></summary>
  Mutations that are susceptible or neutral:
  <ul>
    <li>A110V</li>
    <li>R117L</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_5785</h3></summary>
  Mutations that are susceptible or neutral:
  <ul>
    <li>V260G - by itself</li>
  </ul>
  Single mutations listed in CARD that are actually part of an <em>n</em>-tuple mutation:
  <ul>
    <li>K123Q - part of double mutation K123Q, V260G</li>
  </ul>
  Mutations missing from CARD:
  <ul>
    <li>A143V</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_5803</h3></summary>
  Mutations listed wrong in CARD:
  <ul>
    <li>M1S - is actually a M1I mutation</li>
    <li>A8G - is actually a D8G mutation</li>
    <li>Y68D - is actually a W68D mutation</li>
    <li>L72P - is actually a L172P mutation</li>
    <li>A140S - is actually a R140S mutation</li>
    <li>R157W - is actually a V157W mutation</li>
  </ul>
  Mutations listed in CARD that are not present in articles referenced:
  <ul>
    <li>H59D</li>
    <li>W68D</li>
    <li>V157W</li>
  </ul>
  Mutations that are susceptible or neutral:
  <ul>
    <li>Q10*</li>
    <li>D12A</li>
    <li>D12N</li>
    <li>T47A</li>
    <li>P54L</li>
    <li>H57D</li>
    <li>Y103*</li>
    <li>Y103S</li>
    <li>H137R</li>
    <li>T142M - by itself</li>
    <li>A171V</li>
  </ul>
  Single mutations listed in CARD that are actually part of an <em>n</em>-tuple mutation:
  <ul>
    <li>R140S - part of double mutation with an 8 bp deletion at nt 446<sup><a href="#footnote14" id="reference14">14</a></sup></li>
  </ul>
  <em>N</em>-tuple mutations that actually include at least one single mutation that can confer resistance by itself:
  <ul>
    <li>V130G, GG insertion at nt 420 - V130G is already listed as a resistance-conferring single mutation</li>
  </ul>
  Mutations missing from CARD:
  <ul>
    <li>D8A</li>
    <li>F13SL</li>
    <li>C14*</li>
    <li>L27R</li>
    <li>H43P</li>
    <li>V44G</li>
    <li>K48E</li>
    <li>H51P</li>
    <li>D53N, CACTG insertion at nt 349<sup><a href="#footnote14" id="reference14">14</a></sup></li>
    <li>P54Q</li>
    <li>P54R</li>
    <li>P54S</li>
    <li>H57P</li>
    <li>H57Y</li>
    <li>P62L</li>
    <li>P62Q</li>
    <li>P62T</li>
    <li>Y64D</li>
    <li>Y64*</li>
    <li>H71P</li>
    <li>T72P</li>
    <li>E91*</li>
    <li>F94C</li>
    <li>F94L</li>
    <li>F94S</li>
    <li>Y95*</li>
    <li>K96R</li>
    <li>K96T</li>
    <li>G97A</li>
    <li>G97D</li>
    <li>T100A</li>
    <li>A102P</li>
    <li>S104R</li>
    <li>T114P</li>
    <li>L120R</li>
    <li>V125G</li>
    <li>V130G</li>
    <li>G132A</li>
    <li>G132C</li>
    <li>G132D</li>
    <li>A134G</li>
    <li>T135A, T142M</li>
    <li>T135N</li>
    <li>D136G</li>
    <li>H137D</li>
    <li>R140P</li>
    <li>T142A</li>
    <li>A143G</li>
    <li>R148S</li>
    <li>N149V</li>
    <li>T153N</li>
    <li>R154W</li>
    <li>L156P</li>
    <li>G162V</li>
    <li>M175I</li>
    <li>M175T</li>  
    <li>RSM120-122 insertion</li>
    <li>RSM130-132 insertion</li>
    <li>RSM131-133 insertion</li>
    <li>MWS132-134 insertion</li>
    <li>8 bp deletion at start</li>
    <li>11-bp deletion at start</li>
    <li>TCATCG deletion at nt 14</li>
    <li>C deletion at nt 28</li>
    <li>GACT insertion at nt 37</li>
    <li>C insertion at nt 44</li>
    <li>G insertion at nt 52</li>
    <li>C deletion at nt 59</li>
    <li>G deletion at nt 70</li>
    <li>G deletion at nt 71</li>
    <li>G deletion at nt 77</li>
    <li>C deletion at nt 84</li>
    <li>C deletion at nt 104</li>
    <li>G insertion at nt 136</li>
    <li>T deletion at nt 150</li>
    <li>80 bp deletion at nt 151</li>
    <li>A deletion at nt 158</li>
    <li>C deletion at nt 161</li>
    <li>C deletion at nt 186</li>
    <li>4 bp insertion at nt 185</li>
    <li>A insertion at nt 186</li>
    <li>A insertion at nt 192</li>
    <li>A insertion at nt 193</li>
    <li>TATCAGG insertion at nt 193</li>
    <li>68 bp deletion at nt 195</li>
    <li>CGCATTGCCG insertion at nt 218</li>
    <li>G insertion at nt 222</li>
    <li>T deletion at nt 224</li>
    <li>24 bp deletion at nt 265</li>
    <li>T insertion at nt 287</li>
    <li>T insertion at nt 288</li>
    <li>33 bp insertion at nt 288</li>
    <li>G deletion at nt 290</li>
    <li>T deletion at nt 291</li>
    <li>A insertion at nt 301</li>
    <li>G deletion at nt 301</li>
    <li>TACAG deletion at nt 307</li>
    <li>C insertion at nt 338</li>
    <li>C deletion at 341</li>
    <li>AG insertion at nt 368</li>
    <li>18 bp insertion at nt 368</li>
    <li>11 bp deletion at nt 379</li>
    <li>GG deletion at nt 381</li>
    <li>AG insertion at nt 382</li>
    <li>ATG deletion at nt 386 - D129 deletion</li>
    <li>ATGT deletion at nt 386</li>
    <li>G insertion at nt 391</li>
    <li>GG insertion at nt 391</li>
    <li>G insertion at nt 392</li>
    <li>GG insertion at nt 392</li>
    <li>G insertion at nt 393</li>
    <li>GG insertion at nt 393</li>
    <li>17-bp deletion at nt 395</li>
    <li>T insertion at nt 397</li>
    <li>T deletion at nt 398</li>
    <li>C insertion at nt 403</li>
    <li>CC insertion at nt 403</li>
    <li>G deletion at nt 406</li>
    <li>C insertion at nt 408</li>
    <li>A insertion at nt 408</li>
    <li>104-bp deletion at nt 409</li>
    <li>G insertion at nt 414</li>
    <li>TG deletion at nt 416</li>
    <li>G insertion at nt 417</li>
    <li>G insertion at nt 420</li>
    <li>CAGACGGCGCCAG insertion at nt 423</li>
    <li>GG insertion at nt 428</li>
    <li>G deletion at nt 443</li>
    <li>8 pb deletion at nt 446</li>
    <li>T deletion at nt 452</li>
    <li>T insertion at nt 465</li>
    <li>C insertion at nt 475</li>
    <li>TGAC insertion at nt 480</li>
    <li>GT deletion at nt 487</li>
    <li>C insertion at nt 493</li>
    <li>CG insertion at nt 501</li>
    <li>C deletion at nt 512</li>
    <li>C deletion at nt 514</li>
    <li>CG insertion at nt 516</li>
    <li>5 bp insertion at nt 518</li>
    <li>C insertion at nt 532</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_6045<sup><a href="#footnote15" id="reference15">15</a></sup></h3></summary>
  Mutations missing from CARD:
  <ul>
    <li>E53*</li>
    <li>W89*</li>
    <li>R108*</li>
    <li>transposase insertion</li>
    <li>integrase insertion</li>
    <li>ACAAAGCGAT deletion</li>
    <li>CTCGACGTCGGCCAT deletion</li>
    <li>CACAAAGCGAT deletion</li>
    <li>GC insertion</li>
    <li>G deletion</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_6046<sup><a href="#footnote16" id="reference16">16</a></sup></h3></summary>
  Single mutations listed in CARD that are actually part of an <em>n</em>-tuple mutation:
  <ul>
    <li>E160D - part of double mutation with a 15-bp deletion<sup><a href="#footnote17" id="reference17">17</a></sup></li>
    <li>M84I - part of double mutation T18P, M84I</li>
  </ul>
  <em>N</em>-tuple mutations that actually include at least one single mutation that can confer resistance by itself:
  <ul>
    <li>T18P, M84I - T18P is already listed as a resistance-conferring single mutation</li>
  </ul>
  Mutations missing from CARD:
  <ul>
    <li>Q19*</li>
    <li>E160*</li>
    <li>10 bp deletion</li>
    <li>C insertion at nt 461</li>
  </ul>
</details>

<details>
  <summary><h3>MEG_6090</h3></summary>
  Mutations listed wrong in CARD<sup><a href="#footnote18" id="reference18">18</a></sup>:
  <ul>
    <li>V146F - is actually a V170F mutation</li>
    <li>V176F (from fifth article) - is actually a V170F mutation</li>
    <li>V176F (from 24th article) - is actually a V198F mutation</li>
    <li>A381V - is actually a A286V mutation</li>
    <li>S450L - is actually a S369L mutation</li>
    <li>Q490H - is actually a Q409H mutation</li>
    <li>Q438K - is actually a Q432K mutation</li>
    <li>D441Y - is actually a D435Y mutation</li>
    <li>D441V - is actually a D435V mutation</li>
    <li>S447Q - is actually a S441Q mutation</li>
    <li>H451D - is actually a H445D mutation</li>
    <li>H451R - is actually a H445R mutation</li>
    <li>H451C - is actually a H445C mutation</li>
    <li>H451Y - is actually a H445Y mutation</li>
    <li>S456W - is actually a S450W mutation</li>
    <li>S456L - is actually a S450L mutation</li>
    <li>E504A - is actually a E423A mutation</li>
    <li>F505L - is actually a F424L mutation</li>
    <li>G507S - is actually a G426S mutation</li>
    <li>G507D - is actually a G426D mutation</li>
    <li>T508H - is actually a T427H mutation</li>
    <li>T508P - is actually a T427P mutation</li>
    <li>T508A - is actually a T427A mutation</li>
    <li>T508S - is actually a T427S mutation</li>
    <li>T508N - is actually a T427N mutation</li>
    <li>S509R - is actually a S428R mutation</li>
    <li>S509Q - is actually a S428Q mutation</li>
    <li>L511M - is actually a L430M mutation</li>
    <li>L511P - is actually a L430P mutation</li>
    <li>L511R - is actually a L430R mutation</li>
    <li>L511V - is actually a L430V mutation</li>
    <li>L511P,M515I - is actually a L430P,M434I mutation</li>
    <li>L511R,D516V - is actually a L430R,D435V mutation</li>
    <li>L511P,S512T,D516V - is actually a L430P,S431T,D435V mutation</li>   
    <li>S512I - is actually a S431I mutation</li>
    <li>S512T - is actually a S431T mutation</li>
    <li>S512R - is actually a S431R mutation</li>
    <li>S512N - is actually a S431N mutation</li>
    <li>-S512,-Q513,-F514,H526Q - is actually a -S431,-Q432-F433,H445Q mutation</li>  
    <li>Q513P - is actually a Q432P mutation</li>
    <li>Q513K - is actually a Q432K mutation</li>
    <li>Q513L - is actually a Q432L mutation</li>
    <li>Q513E - is actually a Q432E mutation</li>
    <li>Q513D - is actually a Q432D mutation</li>
    <li>Q513R - is actually a Q432R mutation</li>
    <li>-Q513 - is actually a -Q432 mutation</li>
    <li>Q513STOP - is actually a Q432STOP mutation</li>
    <li>Q513H,-514F,-515M,-516D - is actually a Q432H,-433F,-434M,-435D mutation</li>
    <li>-F514 - is actually a -F433 mutation</li>
    <li>M515I - is actually a M434I mutation</li>
    <li>M515V - is actually a M434V mutation</li>
    <li>-M515 - is actually a -M434 mutation</li>
    <li>-515M,-516D,-517Q,-518N - is actually a -434M,-435D,-436Q,-437N mutation</li>
    <li>D516V - is actually a D435V mutation</li>
    <li>D516T - is actually a D435T mutation</li>
    <li>D516Y - is actually a D435Y mutation</li>
    <li>D516G - is actually a D435G mutation</li>
    <li>D516N - is actually a D435N mutation</li>
    <li>D516H - is actually a D435H mutation</li>
    <li>D516K - is actually a D435K mutation</li>
    <li>-D516 - is actually a -D435 mutation</li>
    <li>D516E,S522L - is actually a D435E,S441L mutation</li>
    <li>D516Y,L511R - is actually a D435Y,L430R mutation</li>
    <li>T516I,G523W,D525Y - is actually a T444I,G442W,D435Y mutation</li>
    <li>Q517L - is actually a Q436L mutation</li>
    <li>Q517H - is actually a Q436H mutation</li>  
    <li>N518T - is actually a N437T mutation</li>
    <li>D518H - is actually a N437H mutation</li>
    <li>N518H - is actually a N437H mutation</li>
    <li>N518I - is actually a N437I mutation</li>
    <li>-N518 - is actually a -N437 mutation</li>
    <li>N519K - is actually a N438K mutation</li>
    <li>-519N,S522L,S531L - is actually a 438N,S441L,S450L mutation</li>
    <li>P520T - is actually a P439T mutation</li>
    <li>L521P - is actually a L452P mutation</li>
    <li>L521M - is actually a L440M mutation</li>
    <li>S522L - is actually a S441L mutation</li>
    <li>S522W - is actually a S441W mutation</li>
    <li>S522Q - is actually a S441Q mutation</li>
    <li>S522STOP - is actually a S441STOP mutation</li>
    <li>G523A - is actually a G442A mutation</li>
    <li>L524S - is actually a L443S mutation</li>
    <li>L524W,T525P,H526Q,-527K - is actually a L443W,T444P,H445Q,-446K mutation</li>
    <li>H526C - is actually a H445C mutation</li>
    <li>H526Y - is actually a H445Y mutation</li>
    <li>H526D - is actually a H445D mutation</li>
    <li>H526S - is actually a H445S mutation</li>
    <li>H526R - is actually a H445R mutation</li>
    <li>H526L - is actually a H445L mutation</li>
    <li>H526N - is actually a H445N mutation</li>
    <li>H526T - is actually a H445T mutation</li>
    <li>H526E - is actually a H445E mutation</li>
    <li>H526G - is actually a H445G mutation</li>
    <li>H526F - is actually a H445F mutation</li>
    <li>H526Q - is actually a H445Q mutation</li>
    <li>H526P - is actually a H445P mutation</li>
    <li>-H526 - is actually a -H445 mutation</li>
    <li>H526S,M515V - is actually a H445S,M434V mutation</li>
    <li>H526S,P535H - is actually a H445S,P454H mutation</li>
    <li>H526Y,E541G - is actually a H445Y,E460G mutation</li> 
    <li>H526D,E541G,S553A - is actually a H445D,E460G,S472A mutation</li>
    <li>H526D,S531 - is actually a H445D,S450 mutation</li>
    <li>H526P,K527Q - is actually a H445P,K446Q mutation</li>
    <li>K527Q - is actually a K446Q mutation</li>
    <li>K527N - is actually a K446N mutation</li>
    <li>R528P - is actually a R447P mutation</li>
    <li>R528H - is actually a R447H mutation</li>
    <li>R529Q - is actually a R448Q mutation</li>
    <li>S531L - is actually a S450L mutation</li>
    <li>S531W - is actually a S450W mutation</li>
    <li>S531G - is actually a S450G mutation</li>
    <li>S531C - is actually a S450C mutation</li>
    <li>S531F - is actually a S450F mutation</li>
    <li>S531L,S622A - is actually a S450L,S541A mutation</li>
    <li>S531L,F514V - is actually a S450L,F433V mutation</li>
    <li>S531L,H526C - is actually a S450L,H445C mutation</li>
    <li>L533R - is actually a L452R mutation</li>
    <li>L533P - is actually a L452P mutation</li>
    <li>L538R - is actually a L457R mutation</li>
    <li>L538P - is actually a L457P mutation</li>
    <li>L545M - is actually a L464M mutation</li>
    <li>E562G,P564L - is actually a E481G,P483L mutation</li>    
    <li>L571V - is actually a L490V mutation</li>
    <li>I572F - is actually a I491F mutation</li>
    <li>S574L - is actually a S493L mutation</li>
    <li>R633C - is actually a R552C mutation</li>
    <li>E672D - is actually a E592D mutation</li>
  </ul>
  Mutations removed due to conflicting evidence:
  <ul>
    <li>S450L,S541A - position 541 doesn't have a serine as its wild-type, and S450L is already a resistance-conferring single mutation</li>
  </ul>
  Mutations listed in CARD that are not present in articles referenced:
  <ul>
    <li>Q432D</li>
    <li>H445E</li>
  </ul>
  Mutations that are susceptible or neutral:
  <ul>
    <li>S369L</li>
    <li>L440P</li>
    <li>S493L</li>
  </ul>
  Single mutations listed in CARD that are actually part of an <em>n</em>-tuple mutation:
  <ul>
    <li>E423A - part of an unspecified double mutation</li>
    <li>F424L - part of double mutation F424L, D435Y</li>
    <li>G426D - part of double mutation G426D, S450L; triple mutation G426D, N437Y, L457P; and quadruple mutation G426D, Q432*<sup><a href="#footnote19" id="reference19">19</a></sup>, H445Y, S450F</li>
    <li>G426S - part of double mutation G426S, T427P </li>
    <li>T427H - part of triple mutations T427H, H445Q, S450F and T427H, G442A, H445R</li>
    <li>T427P - part of double mutation G426S, T427P and triple mutation T427P, D435H, H445R</li>
    <li>T427A - part of double mutation T427A, G442A,</li>
    <li>T427S - part of triple mutation T427S, S431R, H445D</li>
    <li>S428Q - part of double mutation S428Q, S450L</li>
    <li>L430M - part of triple mutation L430M, S431R, H445D</li>
    <li>S431I - part of double mutation S431I, D435G </li>
    <li>S431T - part of double mutation L430P, S431T; triple mutation L430P, S431T, D435V; and quadruple mutation S431T, Q432N, H445F, S450L</li>
    <li>S431R - part of double mutation S431R, D435V and triple mutations L430M, S431R, H445D and T427S, S431R, H445D</li>
    <li>M434I - part of double mutations M434I, D435Y; M434I, L452P; and L430P, M434I</li>
    <li>D435H - part of triple mutation T427P, D435H, H445R</li>
    <li>D435K - part of triple mutation D435K, S450L, P454H</li>
    <li>D435N - part of an unspecified double mutation</li>
    <li>N437H - part of double mutation N437H, S450L</li>
    <li>N438K - part of quadruple mutation L430V, Q432E, N438K, H445R</li>
    <li>H445S - part of double mutations H445S, K446Q; M434V, H445S; and H445S, P454H</li>
    <li>H445F - part of quadruple mutation S431T, Q432N, H445F, S450L</li>
    <li>H445Q - part of double mutations L430P, H445Q and H445Q, S450C; triple mutation T427H, H445Q, S450F; and quadruple mutations L443W, T444P, H445Q, K446- and H445Q, S431-, Q432-, F433-</li>
    <li>K446Q - part of double mutations H445S, K446Q and H445P, K446Q</li>
    <li>R447H - part of double mutation R447H, S450W</li>
    <li>R448Q - part of double mutation D435G, R448Q</li>
    <li>S450C - part of double mutation H445Q, S450C</li>
    <li>L452R - part of double mutation H445C, L452R </li>
    <li>L457P - part of triple mutation G426D, N437Y, L457P</li>
  </ul>
  <em>N</em>-tuple mutations that actually include at least one single mutation that can confer resistance by itself:
  <ul>
    <li>F424L, D435Y - D435Y is already listed as a resistance-conferring single mutation</li>
    <li>G426D, S450L - S450L is already listed as a resistance-conferring single mutation</li>
    <li>T427H, H445Q, S450F - S450F is already listed as a resistance-conferring single mutation</li>
    <li>T427H, G442A, H445R - G442A and H445R are already listed as individual resistance-conferring single mutations</li>
    <li>T427P, D435H, H445R - H445R is already listed as a resistance-conferring single mutation</li>
    <li>T427A, G442A - G442A is already listed as a resistance-conferring single mutation</li>
    <li>T427S, S431R, H445D - H445D is already listed as a resistance-conferring single mutation</li>
    <li>S428Q, S450L - S450L is already listed as a resistance-conferring single mutation</li>
    <li>L430P, S431T - L430P is already listed as a resistance-conferring single mutation</li>
    <li>L430P, S431T, D435V - L430P and D435V are already listed as individual resistance-conferring single mutations</li>
    <li>L430R, D435V - L430R and D435V are already listed as individual resistance-conferring single mutations</li>
    <li>L430P, M434I - L430P is already listed as a resistance-conferring single mutation</li>
    <li>L430M, S431R, H445D - H445D is already listed as a resistance-conferring single mutation</li>
    <li>S431I, D435G - D435G is already listed as a resistance-conferring single mutation</li>
    <li>Q432H, -433F, -434M, -435D - not a single mutation, but -433F,-434M,-435D is already listed as a resistance-conferring triple mutation</li>
    <li>-434M, -435D, -436Q, -437N - -436Q and -437N are already listed as individual resistance-conferring single mutations</li>
    <li>D435Y, L430R - L430R and D435Y are already listed as individual resistance-conferring single mutations</li>
    <li>D435E, S441L - S441L is already listed as a resistance-conferring single mutation</li>
    <li>D435Y, G442W, T444I - D435Y is already listed as a resistance-conferring single mutation</li>
    <li>-438N, S441L, S450L - -438N, S441L and S450L are already listed as individual resistance-conferring single mutations</li>
    <li>H445D, E460G, S472A - H445D is already listed as a resistance-conferring single mutation</li>
    <li>H445D, S450 - H445D is already listed as a resistance-conferring single mutation</li>
    <li>H445P, K446Q - H445P and K446Q are already listed as individual resistance-conferring single mutations</li>
    <li>H445Y, E460G - H445Y is already listed as a resistance-conferring single mutation</li>
    <li>H445S, M434V - M434V is already listed as a resistance-conferring single mutation</li>
    <li>S450L, H445C - H445C and S450L are already listed as individual resistance-conferring single mutations</li>
    <li>S450L, F433V - S450L is already listed as a resistance-conferring single mutation</li>
    
    
    
    
    <li>S431T, Q432N, H445F, S450L - S450L is already listed as a resistance-conferring single mutation</li>
    <li>InsertHere - InsertHere is already listed as a resistance-conferring single mutation</li>
    <li>InsertHere - InsertHere is already listed as a resistance-conferring single mutation</li>
    <li>InsertHere - InsertHere is already listed as a resistance-conferring single mutation</li>
    <li>InsertHere - InsertHere is already listed as a resistance-conferring single mutation</li>
    <li>InsertHere - InsertHere is already listed as a resistance-conferring single mutation</li>
    <li>InsertHere - InsertHere is already listed as a resistance-conferring single mutation</li>
  </ul>
  Mutations missing from CARD:
  <ul>
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
<sup id="footnote13">
  13. All mutations in CARD are not resistance-conferring; therefore the gene is not present in the SNPInfo database <a href="#reference13">↩</a>
</sup><br>
<sup id="footnote14">
  14. In an <em>n</em>-tuple mutation, nonsense and frameshift mutations take precedence and are therefore considered single mutations <a href="#reference14">↩</a>
</sup><br>
<sup id="footnote15">
  15. Positions of insertions and deletions not provided <a href="#reference15">↩</a>
</sup><br>
<sup id="footnote16">
  16. For some insertions and deletions, positions are not provided <a href="#reference16">↩</a>
</sup><br>
<sup id="footnote17">
  17. In an <em>n</em>-tuple mutation, nonsense and frameshift mutations take precedence and are therefore considered single mutations <a href="#reference17">↩</a>
</sup><br>
<sup id="footnote18">
  18.Almost all mutations recovered from the literature are based on E. coli positions <a href="#reference18">↩</a>
</sup><br>
<ul>
  <sup><li>The exception includes all but the last three SNPs recovered from the fifth article</li></sup>
</ul>
<ul>
  <sup><li>However, the M. tuberculosis reference sequence used in that article has 6 extra amino acids at the beginning of the sequence compare to CARD's reference sequence</li></sup>
</ul>
<ul>
  <sup><li>The positions that are used in this file is based on the M.tubercolosis CARD reference sequence</li></sup>
</ul><br>
<sup id="footnote19">
  19. In an <em>n</em>-tuple mutation, nonsense and frameshift mutations take precedence and are therefore considered single mutations <a href="#reference19">↩</a>
</sup><br>
