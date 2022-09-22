# Disclaimer
This version of the genes-issues.md is outdated. However, because the new version of the genes-issues.md is still being updated, this file is kept as reference only. 

# Genes Issues
## Genes that went through the MEG_6090_SNP_Alignment program:
- MEG_3241: The organism used to determine the SNPs was not the same as the one in CARD. Additionally, the sequence in Megares match the sequence in CARD, and there was too big of a difference between the Megares sequence and the sequence that the SNPs were based of. This program created an alignment between the two sequences to figure out the SNPs' positions based on the CARD/Megares sequence.
  - There were additionally three errors in Kargva, where the single mutations E466V, E466K, and R389P were instead written as E467V, E467K, and R388P respectively
  - There was also an error in Kargva, where the single mutation L444F was instead written as K444F
- MEG_3246: The sequence used to determine the SNPs was an older version of the one in CARD. Additionally, the sequence in Megares match the sequence in CARD, and there was too big of a difference between the Megares sequence and the sequence that the SNPs were based of. This program created an alignment between the two sequences to figure out the SNPs' positions based on the CARD/Megares sequence.
  - There was also an error in Kargva, where the single mutation G124S was instead written as G125S
- MEG_3243: the SNP Verification program only looks up to 30 amino acid away from the position provided by the SNP in order to match it to the Megares sequence. For MEG_3243, where the CARD sequence didn't match the Megares sequence, the distance was too small.
  - This program created an alignment between the two sequences to figure out the affected SNPs’ positions based on the Megares sequence.
  - I replaced the CARD sequence with the Megares sequence in the SNPInfoExtraction program to find the context
-	MEG_6090: an overwhelming majority of SNPs’ positions are based on the E.coli rpoB gene instead of the M.tuberculosis rpoB gene. This program created an alignment between the two genes to figure out the affected SNPs’ positions based on the M.tuberculosis rpoB gene.

## Genes where SNPs had to be changed:

- MEG_2710: I found an error in CARD, where the double mutation A313G,Y319C was instead written as A314G,Y322C
- MEG_3594: The organism used to determine the SNPs was not the same as the one in CARD, so I had to change the sequence used to find the context
- MEG_4057: The sequence used in CARD was outdated, so I had to change the sequence used to find the context
-	MEG_4296: In both Kargva and CARD, the SNP is entered as Q142X; however, when looking at the paper where the SNP was found, I found that Q142X actually meant Q142*
-	MEG_5325: When looking at other parC genes, I realized that SNP S80L had the wrong position and is actually SNP S87L
-	MEG_5406: When looking at the paper where all three SNPs were found, I realized that the positions are based on another version of the sequence that had 6 less amino acids at the beginning, so I simply added 6 to the SNPs’ positions.
- MEG_5803: I found five errors in CARD:
  - D8G was written as A8G
  - W68D was written as Y68D
  - R140S was written as A140S
  - V157W was written as R157W
  - L172P was written as L72P
-	MEG_7301: The strain used to determine the SNPs was not the same as the one in CARD, so I had to change the sequence used to find the context
-	MEG_7333: The organism used to determine the SNPs was not the same as the one in CARD, so I had to change the sequence used to find the context

## Genes where SNPs had to be removed:

- MEG_2712: I found an error in CARD where the single mutation S244T actually comes from a different gene
- MEG_3065: I couldn’t find the true positions for the single mutation F441Y and for the double mutation T387I,E449K
- MEG_3243: I hade to remove three single mutations, given below with the CARD positions:
  - Looking at the article where N510D was found, I realized that the actual position of the SNP was 538. Since the SNP N538D is already accounted for, I decided to just remove N510D
  - I couldn’t find the true positions for the single mutations E501D and N499T.
- MEG_3446: I removed two mutations:
  - A234G: the sequence in CARD has a G at that location; however, when looking at the catalase/peroxidase domain in other organisms, I found that G is the wild-type in other organisms as well, which made me feel unsure of the actual validity of this SNP
  - A431V: the sequence in CARD has a V at that location; additionally, CARD also has V431A as a SNP; to verify which SNP was correct, I looked at the catalase/peroxidase domain in other organisms and found that V is the wild-type in other organisms.
- MEG_4057: this was one of the genes where I had to change the sequence in order to find the context; the problem was that the sequence in Megares matches the sequence in CARD. Both of those sequence had deletions at codons 349 to 352, while the sequence used to find the context of SNPs does not. This led to the SNP at codon 351 being automatically disregarded by the SNP Verification program.
  - The SNP at codon 347 was also disregarded by the program due to the fact that the context was now too different due to the deletions, so I tweaked the code so that it was forced to accept this SNP.
- MEG_6090: the alignment between the E.coli rpoB gene and the M.tuberculosis rpoB gene led to a couple of SNPs with “unaligned” E.coli positions due to the fact that the amino acids at those positions in E.coli were different from the wild-type amino acid in the SNP; however, when using the surrounding amino acids that were aligned, I was able to figure out all M.tuberculosis positions except for one -> S622A. 
  - Since S622A is part of a double mutation, S531L,S622A (using E.coli positions), and the single mutation S531L is a valid resistance-conferring SNP, I decided to simply remove the double mutation
- MEG_7250: I couldn’t find the true position for the single mutation I211V.
- MEG_7309: as a shorter TUFAB gene of the length of 300 aa, it didn't contain SNPs that are normally found in position 320 and after in other TUFAB genes
