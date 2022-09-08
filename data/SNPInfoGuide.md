# SNP Info Guide - Understanding the SNP Info annotation scheme
## Table of Contents
- [What does *Nuc* mean?](https://github.com/Isabella136/AmrPlusPlus_SNP/new/main/data#what-does-nuc-mean)
- [The *Mis* notation, or the original SNP](https://github.com/Isabella136/AmrPlusPlus_SNP/new/main/data#the-mis-notation-or-the-original-snp)
- [*Nonsense*, an unexpected stop codon](https://github.com/Isabella136/AmrPlusPlus_SNP/new/main/data#nonsense-an-unexpected-stop-codon)
- [*Nonstop* and the reason why we don't have a *Non* notation](https://github.com/Isabella136/AmrPlusPlus_SNP/new/main/data#nonstop-and-the-reason-why-we-dont-have-a-non-notation)
- [*Hyper*: why nonresistant SNPs are included](https://github.com/Isabella136/AmrPlusPlus_SNP/new/main/data#hyper-why-nonresistant-snps-are-included)
- [Straying away from SNPs with *Ins* and *Del*](https://github.com/Isabella136/AmrPlusPlus_SNP/new/main/data#straying-away-from-snps-with-ins-and-del)
- [*Mult* because not all mutations come in ones](https://github.com/Isabella136/AmrPlusPlus_SNP/new/main/data#mult-because-not-all-mutations-come-in-ones)
- [So what about *Must*?](https://github.com/Isabella136/AmrPlusPlus_SNP/new/main/data#so-what-about-must)
- [*FS*, the new label that requires confirmation](https://github.com/Isabella136/AmrPlusPlus_SNP/new/main/data#so-what-about-must)
## What does *Nuc* mean?
| Accession Number | Type  | Class                 | Mechanism                                              | Group  | Mutations                |                         |                        |                      |
|------------------|-------|-----------------------|--------------------------------------------------------|--------|--------------------------|-------------------------|------------------------|----------------------|
| MEG_1            | Drugs | Aminoglycosides       | Aminoglycoside-resistant_16S_ribosomal_subunit_protein | A16S   | Nuc:A1192G_GGATG_CGTCA   | Nuc:C1193GT_GATGA_GTCAA | Nuc:G1194C_ATGAC_TCAAG |                      |
| MEG_1730         | Drugs | Lipopeptides          | Daptomycin-resistant_mutant                            | CLS    | Mis:L52F_AWLLV_VFLPL     | Mis:T33N_LLFAF_IIFME    | Mis:F60S_LPLFG_ILYLL   | Mis:A23V_SIFIG_FILNL |
| MEG_3740         | Drugs | Multi-drug_resistance | MDR_23S_rRNA_mutation                                  | MDR23S | Nuc:A2071CGT_GACGG_AAGAC | Nuc:A2072GC_ACGGA_AGACC |                        |                      |

*Nuc* stands for nucleic; although the .sam input contains DNA reads, most mutations in protein-coding genes such as MEG_1730 are reported in their amino-acid form.  
That being said, not all genes in MEGARes are proteing-coding, as seen with MEG_1 and MEG_3740.
Indeed, the 16S and 23S subunits of the ribosome are RNA sequences, and it would therefore be useless to use the amino-acid form when describing mutations.  
To identify which mutations belong to rRNA-coding sequences, the SNP_Verification program needs to recognize the *Nuc* notation at the start of a mutation label.
## The *Mis* notation, or the original SNP

## *Nonsense*, an unexpected stop codon
## *Nonstop* and the reason why we don't have a *Non* notation
## *Hyper*: why nonresistant SNPs are included
## Straying away from SNPs with *Ins* and *Del*
## *Mult* because not all mutations come in ones
## So what about *Must*?
## *FS*, the new label that requires confirmation
