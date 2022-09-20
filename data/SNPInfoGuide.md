# SNP Info Guide - Understanding the SNP Info annotation scheme
## Table of Contents
- [What does *Nuc* mean?](https://github.com/Isabella136/AmrPlusPlus_SNP/blob/main/data/SNPInfoGuide.md#what-does-nuc-mean)
- [The *Mis* notation, or the original SNP](https://github.com/Isabella136/AmrPlusPlus_SNP/blob/main/data/SNPInfoGuide.md#the-mis-notation-or-the-original-snp)
- [*Nonsense*, an unexpected stop codon](https://github.com/Isabella136/AmrPlusPlus_SNP/blob/main/data/SNPInfoGuide.md#nonsense-an-unexpected-stop-codon)
- [*Nonstop* and and its different reasons for existing](https://github.com/Isabella136/AmrPlusPlus_SNP/blob/main/data/SNPInfoGuide.md#nonstop-and-its-different-reasons-for-existing)
- [*Hyper*: why nonresistant SNPs are included](https://github.com/Isabella136/AmrPlusPlus_SNP/blob/main/data/SNPInfoGuide.md#hyper-why-nonresistant-snps-are-included)
- [Straying away from SNPs with *Ins* and *Del*](https://github.com/Isabella136/AmrPlusPlus_SNP/blob/main/data/SNPInfoGuide.md#straying-away-from-snps-with-ins-and-del)
- [*Mult* because not all mutations come in ones](https://github.com/Isabella136/AmrPlusPlus_SNP/blob/main/data/SNPInfoGuide.md#mult-because-not-all-mutations-come-in-ones)
- [So what about *Must*?](https://github.com/Isabella136/AmrPlusPlus_SNP/blob/main/data/SNPInfoGuide.md#so-what-about-must)
- [*FS*, the new label that requires confirmation](https://github.com/Isabella136/AmrPlusPlus_SNP/blob/main/data/SNPInfoGuide.md#fs-the-new-label-that-requires-confirmation)

## What does *Nuc* mean?
| Accession Number | Type  | Class                 | Mechanism                                              | Group  | Mutations                |                         |                        |                      |
|------------------|-------|-----------------------|--------------------------------------------------------|--------|--------------------------|-------------------------|------------------------|----------------------|
| MEG_1            | Drugs | Aminoglycosides       | Aminoglycoside-resistant_16S_ribosomal_subunit_protein | A16S   | Nuc:A1192G_GGATG_CGTCA   | Nuc:C1193GT_GATGA_GTCAA | Nuc:G1194C_ATGAC_TCAAG |                      |
| MEG_1730         | Drugs | Lipopeptides          | Daptomycin-resistant_mutant                            | CLS    | Mis:L52F_AWLLV_VFLPL     | Mis:T33N_LLFAF_IIFME    | Mis:F60S_LPLFG_ILYLL   | Mis:A23V_SIFIG_FILNL |
| MEG_3740         | Drugs | Multi-drug_resistance | MDR_23S_rRNA_mutation                                  | MDR23S | Nuc:A2071CGT_GACGG_AAGAC | Nuc:A2072GC_ACGGA_AGACC |                        |                      |

*Nuc* stands for nucleic:
- Although the .sam input contains DNA reads, most mutations in protein-coding genes such as MEG_1730 are reported in their amino-acid form

That being said, not all genes in MEGARes are proteing-coding, as seen with MEG_1 and MEG_3740:
- The 16S and 23S subunits of the ribosome are RNA sequences, and it would therefore be useless to use the amino-acid form when describing mutations
- To identify which mutations belong to rRNA-coding sequences, the SNP_Verification program needs to recognize the *Nuc* notation at the start of a mutation label
- Some rRNA gene have been shown to contain more than one resistance-conferring mutation at the same base
    - Base pair 2071 in MEG_3740 can either be mutated in a C, a G, or a T

## The *Mis* notation, or the original SNP
| Accession Number | Type  | Class            | Mechanism                                              | Group  | Mutations                |                         |                       |                       |                       |                       |
|------------------|-------|------------------|--------------------------------------------------------|--------|--------------------------|-------------------------|-----------------------|-----------------------|-----------------------|-----------------------|
| MEG_3178         | Drugs | Fluoroquinolones | Fluoroquinolone-resistant_DNA_topoisomerases           | GYRA   | Mis:S83L_HPHGD_AVYDT     |                         |                       |                       |                       |                       |
| MEG_3244         | Drugs | Aminocoumarins   | Aminocoumarin-resistant_DNA_topoisomerases             | GYRBA  | Mis:R136CHS_ELVIQ_EGKIH  | Mis:G164V_ETEKT_TMVRF   |                       |                       |                       |                       |
| MEG_6136         | Drugs | Rifampin         | Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB | RPOB   | Mis:H481N_LAELT_KRRLS    | Mis:Q468K_SSQLS_FMDQA   | Mis:A473T_QFMDQ_NPLAE | Mis:Q465R_FFGSS_LSQFM | Mis:L466S_FGSSQ_SQFMD | Mis:A477T_QANPL_ELTHK |

Single-nucleotide polymorphism, or SNPs for short, are exactly what their name implies—a mutation that changed one base pair from the original sequence:
- At times, a SNP can cause a missense mutation which would lead to a change in one amino acid in the protein sequence
    - The S83L missense mutation in MEG_3178 was caused by a SNP in the DNA sequence
- Some of those genes have more than one missense mutation that results in antimicrobial resistance
    - MEG_6136 can have one of six different missense mutations
- Some of those genes also have more than one resistance-conferring missense mutation at the same position
    - R136 in MEG_3244 can either be mutated into a C, an H, or an S

## *Nonsense*, an unexpected stop codon
| Accession Number | Type  | Class      | Mechanism                  | Group | Mutations                 |                            |                            |                       |                       |
|------------------|-------|------------|----------------------------|-------|---------------------------|----------------------------|----------------------------|-----------------------|-----------------------|
| MEG_4094         | Drugs | Fosfomycin | Fosfomycin_target_mutation | MURA  | Nonsense:L42*_DKPSK_VNVPA |                            |                            |                       |                       |
| MEG_7329         | Drugs | Fosfomycin | Fosfomycin_target_mutation | UHPT  | Mis:G112E_TVLIM_FVLSY     | Nonsense:W228*_RAEEI_EEPVD | Nonsense:Y314*_SLLWG_VSDLL | Mis:G358V_SLFAL_ALIFG | Mis:W425R_YTLSG_TDVFI |

Nonsense mutations involves SNP that results in a change from an amino acid to a stop in the protein sequence:
- This is the case for MEG_4094 where L42 can be changed into a stop codon, thus abruptly terminating the sequence too early
- Genes are not limited in the number of resistance-conferring nonsense mutations they can have
    - For example, MEG_7329 can have a nonsense mutation at W228 or Y314
- Genes that can have nonsense mutations can also have other types of mutations
    - For example, MEG_7329 which can also have one of three missense mutations

## *Nonstop* and its different reasons for existing
| Accession Number | Type  | Class                                    | Mechanism                     | Group | Mutations                              |
|------------------|-------|------------------------------------------|-------------------------------|-------|----------------------------------------|
| MEG_3594         | Drugs | Multi-drug_resistance                    | Multi-drug_ABC_efflux_pumps   | LMRA  | Mult:Mis:Q52P_PGGKE_LAIEA;Nonstop:189S |
| MEG_6142         | Drugs | Mycobacterium_tuberculosis-specific_Drug | Pyrazinamide-resistant_mutant | RPSA  | Nonstop:482fs                          |

Nonstop mutations involves a change from the stop codon to an amino acid. There are two ways to cause a nonstop mutation:
- A SNP at one of the stop codon's base
    - The nonstop mutation in MEG_3594 involves changing the stop codon to a serine
- A frameshift caused by an insertion or deletion of a base at the stop codon
    - The nonstop mutation in MEG_6142 involves such a frameshift, although it is currently unclear whether the frameshift is an insertion or a deletion

## *Hyper*: why nonresistant SNPs are included
| Accession Number | Type  | Class            | Mechanism                                    | Group | Mutations            |                       |                      |                           |                      |                                                 |                                                 |
|------------------|-------|------------------|----------------------------------------------|-------|----------------------|-----------------------|----------------------|---------------------------|----------------------|-------------------------------------------------|-------------------------------------------------|
| MEG_3180         | Drugs | Fluoroquinolones | Fluoroquinolone-resistant_DNA_topoisomerases | GYRA  | Mis:G88C_NYHPH_DASIY | Mis:D89GN_YHPHG_ASIYD | Mis:S91P_PHGDA_IYDSL | Mis:D94AYGNTH_DASIY_SLVRM | Mis:A74S_SHAKS_RSVAE | Mult:Mis:A90V_HPHGD_SIYDS;Mis:P102H_VRMAQ_WSLRY | Hyper:Mis:T80A_RSVAE_MGNYH;Mis:A90G_HPHGD_SIYDS |

MEG_3180 has a list of possible resistance-conferring missense mutations; however, it can also have a hypersusceptible double mutation:
- The presence of both T80A and A90G leads to an increase in the gene's susceptibility to fluoroquinolones
- Individually, T80A and A90G are also susceptible mutations, but articles cited in the CARD entry of this gene have shown that the combination of a resistance-conferring mutation and a susceptible mutation still leads to higher MIC levels
- However, because there is no data on the combination of one or more resistance-conferring mutations and the hypersusceptible double mutation, the SNP_Verification program identifies such instances as susceptible

## Straying away from SNPs with *Ins* and *Del*
| Accession Number | Type  | Class         | Mechanism                                              | Group  | Mutations                  |
|------------------|-------|---------------|--------------------------------------------------------|--------|----------------------------|
| MEG_6090         | Drugs | Rifampin      | Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB | RPOB   | Ins:R434_QLSQF_MDQNN       |
| MEG_6094         | Drugs | Rifampin      | Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB | RPOB   | Del:A532_KRRIS_LGPGG       |
| MEG_6963         | Drugs | Tetracyclines | Tetracycline-resistant_16S_ribosomal_subunit_protein   | TET16S | NucDel:C595_CTCCA_CGCTT    |

Despite what its name suggests, the SNP_Verification program is now also able to identify possible resistance-conferring indels:
- Gene MEG_6090 lists the insertion of arginine at position 434 as a resistance-conferring mutation
- Gene MEG_6094 lists the deletion of alanine at position 532 as a resistance-conferring mutation
- rRNA such as MEG_6963 can also have resistance-conferring deletions [^1]

| Accession Number | Type  | Class        | Mechanism                   | Group  | Mutations                  |                                      |                                  |
|------------------|-------|--------------|-----------------------------|--------|----------------------------|--------------------------------------|----------------------------------|
| MEG_1731         | Drugs | Lipopeptides | Daptomycin-resistant_mutant | CLS    | Del:K59/60/61_RGLTD_FYLQQ  |                                      |                                  |
| MEG_1732         | Drugs | Lipopeptides | Daptomycin-resistant_mutant | CLS    | Ins:MPL108/111_SSLNR_TRMNS |                                      |                                  |
| MEG_3588         | Drugs | Lipopeptides | Daptomycin-resistant_mutant | LIAFSR | Mis:T194I_LLEHS_FYGTV      | Ins:I168/169/170/171/172_PKEDN_RKGFG | Del:I168/169/170/171_PKEDN_RKGFG |

In the SNPInfo database, indels have special rules concerning the position listed:
- For MEG_1731 which contains the sequence **RGLTDKKKFYLQQ** from positions 54 to 66, any of the three lysines can be marked as deleted in the SAM file input; therefore all three positions are included in the annotation
- For MEG_1732 which contains the sequence **SSLNRMPLTRMNS** from positions 103 to 115, an extra **MPL** can be inserted either before or after the **MPL** already there; therefore both the position of **M** or the position right after **L** are included in the annotation
    - Additionally, it is possible to have for example an **PLM** or a **LMP** insertion after the original **M** and **P** respectively
        - Although this possibility is not listed in the database, it is accounted for by the SNP_Verification program  

Finally, as seen with MEG_3588, genes are able to have both SNPs and Indels as potential resistance-conferring mutations

[^1]: Currently no instances of resistance-conferring insertions have been found in rRNA

## *Mult* because not all mutations come in ones
Sometimes, genes require a *n*-tuple mutation (i.e. a combination of two or more specific mutations) in order to be resistant to an antimicrobial:
- Such combinations are identified by the SNP_Verification program through the *Mult* notation
- The individual mutations that are part of a combination are divided by semi-colons

| Accession Number | Type           | Class                                    | Mechanism                             | Group  | Mutations                                                              |                                                |                      |
|------------------|----------------|------------------------------------------|---------------------------------------|--------|------------------------------------------------------------------------|------------------------------------------------|----------------------|
| MEG_411          | Multi-compound | Drug_and_biocide_resistance              | Drug_and_biocide_RND_efflux_regulator | ACRR   | Mult:Mis:L204R_LEMYL_CPTLR;Del:C205_EMYLL_PTLRN                        | Mis:R45PC_AAGVT_GAIYW                          | Mis:P85Q_KFPGD_LSVLR |
| MEG_3446         | Drugs          | Mycobacterium_tuberculosis-specific_Drug | Isoniazid-resistant_mutant            | KATG   | Mult:Mis:K433Q_GPLVP_QTLLW;Mis:Q434A_PLVPK_TLLWQ;Mis:T435D_LVPKQ_LLWQD | Mult:Del:W191_GRVDQ_EPDEV;Del:E192_RVDQW_PDEVY |                      |
| MEG_3594         | Drugs          | Multi-drug_resistance                    | Multi-drug_ABC_efflux_pumps           | LMRA   | Mult:Mis:Q52P_PGGKE_LAIEA;Nonstop:189S                                 |                                                |                      |
| MEG_3976         | Drugs          | MLS                                      | Macrolide-resistant_23S_rRNA_mutation | MLS23S | NucMult:Nuc:C2252T_GCCAT_CTTGA;Nuc:G2291A_AACTG_CCTGT                  |                                                |                      |
| MEG_5780         | Drugs          | Lipopeptides                             | Daptomycin-resistant_mutant           | PGSA   | Mult:Mis:K75N_VTNMG_FLDPL;Ins:GE76_TNMGK_FLDPL                         |                                                |                      |

A couple of things to note about the *n*-tuple mutations listed in the SNPInfo database:
- A gene can have both *n*-tuple or single mutations listed as resistance-conferring, as seen with MEG_411
- A gene can have two different *n*-tuple mutations listed as resistance-conferring, as seen with MEG_3446
- *n*-tuple mutations can be made of more than two mutations, as seen with the first *n*-tuple mutation listed under MEG_3446
- *n*-tuple mutations don't have to be made solely of missense mutations[^2]:
    - The *n*-tuple mutation in MEG_411 has a missense and a deletion
    - The *n*-tuple mutation in MEG_3594 has a missense and a nonstop
    - The *n*-tuple mutation in MEG_5780 has a missense and an insertion
    - The second *n*-tuple mutation in MEG_3446 has two deletions
- rRNA genes like MEG_5780 can also have *n*-tuple mutations and are noted with the *NucMult* notation [^3]

[^2]: Out of all of the *n*-tuple mutations currently identified, none contain nonsense mutations, none contain a combination of deletion and nonstop, insertion and nonstop, or insertion and deletion, and none are made of only insertions 
[^3]: Currently only one rRNA gene, MEG_5780, has a resistance-conferring *n*-tuple mutation

## So what about *Must*?
| Accession Number | Type  | Class    | Mechanism                                              | Group  | Mutations                                                                                                                                                                                    |
|------------------|-------|----------|--------------------------------------------------------|--------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| MEG_1628         | Drugs | Phenicol | Chloramphenicol_acetyltransferases                     | CATB   | Must:Amino:H79_AGNQG_KYDWI                                                                                                                                                                   |
| MEG_3983         | Drugs | MLS      | Macrolide-resistant_23S_rRNA_mutation                  | MLS23S | Must:Nuc:T830_TCTTT_CGTCT                                                                                                                                                                    |
| MEG_6092         | Drugs | Rifampin | Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB | RPOB   | Must:Amino:M414_PVVAA_KEFFG;Amino:R430_QFMDQ_NPLAS;Amino:A434_QRNPL_SLTNK;Amino:S435_RNPLA_LTNKR;Amino:N438_LASLT_KRRLS;Amino:Y464_VRDVH_SHYGR;Amino:M484_PNIGL_GYLSV;Amino:Y486_IGLMG_LSVYA |

Although most SNPConfirmation genes require SNPs and other types of mutations in order to be resistant to antimicrobials, a few genes are the opposite—the presence of SNPs or other types of mutations make them susceptible to anitmicrobials
- Those genes have a list of one or more nucleotides/amino acids that *must* remain unchanged in order to keep being resistant
    - This list is identified by the SNP_Verification program through the *Must* notation
    - Because this list doesn't concern mutants, only the wild-types and their respective positions are included
    - Nucleotides/amino acids in the list are seperated by semi-colons
- In the SNPInfo database, both protein-coding genes and rRNA genes can be intrinsically resistant:
    - For intrinsically-resistant protein-coding genes such as MEG_1628 and MEG_6092, each amino acid listed is preceded with the *Amino* notation
    - For intrinsically-resistant rRNA genes such as MEG_3983, each nucleotide listed is preceded with the *Nuc* notation
- The *Must* list can contain either one element (e.g. MEG_1628 and MEG_3983) or multiple elements (e.g. MEG_6092)

| Accession Number | Type  | Class | Mechanism                             | Group  | Mutations             |                       |                                                |
|------------------|-------|-------|---------------------------------------|--------|-----------------------|-----------------------|------------------------------------------------|
| MEG_3979         | Drugs | MLS   | Macrolide-resistant_23S_rRNA_mutation | MLS23S | Nuc:G270A_GATAG_AACCA | Nuc:T823C_GGTCT_TTCGT | Must:Nuc:A271_ATAGG_ACCAA;Nuc:T825_TCTTT_CGTCT |

Currently, the SNPInfo database contains only one gene that is both intrinsically resistant to members of an antimicrobial class and can gain acquired resistant to more members of that class:
- For MEG_3979 to be considered resistant by the SNP_Verification program, either the two nucleotides listed in *Must* must be unchanged, or one of the two SNPs must be found

## *FS*, the new label that requires confirmation
When researching for resistance-conferring SNPs and indels, information about frameshift mutations for 27 genes have been found:
- Frameshift mutations occurs in protein-coding genes when a certain *x* amount of bases, where *x* is defined as a positive integer not divisible by three, are inserted or deleted in the original DNA sequence
- This usually results in a non-functional protein either due to early termination or due to changed conformation/structure in the new protein sequence

| Accession Number | Type           | Class                                    | Mechanism                                  | Group | Mutations                |                            |
|------------------|----------------|------------------------------------------|--------------------------------------------|-------|--------------------------|----------------------------|
| MEG_413          | Multi-compound | Drug_and_biocide_resistance              | Drug_and_biocide_RND_efflux_regulator      | ACRR  | FS-pump_repressor        |                            |
| MEG_2710         | Drugs          | Mycobacterium_tuberculosis-specific_Drug | Ethambutol_resistant_arabinosyltransferase | EMBB  | FS-miscellaneous         | Mis:N13S_RKSTP_RAILG       |
| MEG_2866         | Drugs          | Mycobacterium_tuberculosis-specific_Drug | Ethionamide-resistant_mutant               | ETHA  | FS-drug_activator        | Mis:G43CS_RESMG_TWDLF      |
| MEG_3135         | Drugs          | Aminoglycosides                          | Aminoglycoside-resistant_gidB              | GIDB  | FS-16S_methyltransferase | Mis:R47Q_GRLWD_HLLNC       |
| MEG_3143         | Drugs          | Fosfomycin                               | Fosfomycin_target_mutation                 | GLPT  | FS-transporter           | Mis:W137R_FQGMG_PPSGR      |
| MEG_3626         | Drugs          | Lipopeptides                             | Colistin-resistant_mutant                  | LPXA  | FS-target_synthesis      | Mis:G68D_QFASV_EVCQD       |
| MEG_4296         | Drugs          | betalactams                              | Mutant_porin_proteins                      | OPRD  | FS-porin                 | Nonsense:Q142*_KWGEM_PTAPV |

Those 27 genes are identified with the *FS* label, and 25 of them are considered resistant by the SNP_Verification program if they have a frameshift mutation or a nonsense mutation:
- For those 25 genes, the *FS* label also contains one of seven potential descriptions:
    - "Pump_repressor", "drug_activator", "16S_methyltransferase", "transporter", "target_synthesis", and "porin" describes the mechanism that gets disrupted by frameshift or nonsense mutations
        - Those listed mechanisms have been found to facilitate susceptibility to antimicrobials; rendering one of the involved proteins to a non-functional state would increase antimicrobial resistance
    - The "miscellaneous" description is reserved for genes that are usually considered essential to maintining important cells functions
        - Although frameshift mutitions in those genes have been found to confer resistance, further reasearch as to how such mutations affect overall cell functions must be done
- Those genes can also have resistance-conferring mutations instead of frameshifts, as seen for example with MEG_2710
- A few genes like MEG_413 have only demonstrated resistance through frameshift mutations and therefore don't have any non-frameshift mutations listed in the SNPInfo database

| Accession Number | Type  | Class                                    | Mechanism                                              | Group | Mutations                             |                       |
|------------------|-------|------------------------------------------|--------------------------------------------------------|-------|---------------------------------------|-----------------------|
| MEG_6094         | Drugs | Rifampin                                 | Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB | RPOB  | FS-followed_by_frameshift_suppression | Mis:R687H_GANMQ_QAVPT |
| MEG_6142         | Drugs | Mycobacterium_tuberculosis-specific_Drug | Pyrazinamide-resistant_mutant                          | RPSA  | FS-stop_codon_pos_26*                 | Mis:D342N_VVAVG_DAMVK |

The remaining two genes can only accept frameshift during special cases:
- MEG_6094 can go through a C insertion between nucleotide positions 1592 and 1594
    - The insertion will later get supressed during translation in a process described in this [paper](https://pubmed.ncbi.nlm.nih.gov/31992637/)
    - The insertion/frameshift suppression combination has been shown to result in rifampin resistance
- MEG_6142 can have an A deletion at nucleotide position 76, which results in a early stop codon
    - The frameshift hasn't been shown to confer resistance by itself; however the gene still gets expressed has seen by the fact that one of the resistant strains in this [paper](http://www.ncbi.nlm.nih.gov/pubmed/30402208) contained this mutation
