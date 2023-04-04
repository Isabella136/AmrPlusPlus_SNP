#   AMRPlusPlus_SNP_Verification
#   Copyright (C) 2022  Nathalie Bonin
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see https://www.gnu.org/licenses/.

from SNP_Verification_Tools import Gene, reverseTranslation
from Bio import SeqIO
import configparser
import pysam
import os

header = """@HD	VN:1.6	SO:coordinate
@SQ	SN:MEG_1|Drugs|Aminoglycosides|Aminoglycoside-resistant_16S_ribosomal_subunit_protein|A16S|RequiresSNPConfirmation	LN:1507
@SQ	SN:MEG_2|Drugs|Aminoglycosides|Aminoglycoside-resistant_16S_ribosomal_subunit_protein|A16S|RequiresSNPConfirmation	LN:1544
@SQ	SN:MEG_3|Drugs|Aminoglycosides|Aminoglycoside-resistant_16S_ribosomal_subunit_protein|A16S|RequiresSNPConfirmation	LN:1537
@SQ	SN:MEG_4|Drugs|Aminoglycosides|Aminoglycoside-resistant_16S_ribosomal_subunit_protein|A16S|RequiresSNPConfirmation	LN:1551
@SQ	SN:MEG_5|Drugs|Aminoglycosides|Aminoglycoside-resistant_16S_ribosomal_subunit_protein|A16S|RequiresSNPConfirmation	LN:1504
@SQ	SN:MEG_6|Drugs|Aminoglycosides|Aminoglycoside-resistant_16S_ribosomal_subunit_protein|A16S|RequiresSNPConfirmation	LN:1551
@SQ	SN:MEG_7|Drugs|Aminoglycosides|Aminoglycoside-resistant_16S_ribosomal_subunit_protein|A16S|RequiresSNPConfirmation	LN:1544
@SQ	SN:MEG_8|Drugs|Aminoglycosides|Aminoglycoside-resistant_16S_ribosomal_subunit_protein|A16S|RequiresSNPConfirmation	LN:1474
@SQ	SN:MEG_9|Drugs|Aminoglycosides|Aminoglycoside-resistant_16S_ribosomal_subunit_protein|A16S|RequiresSNPConfirmation	LN:1528
@SQ	SN:MEG_10|Drugs|Aminoglycosides|Aminoglycoside-resistant_16S_ribosomal_subunit_protein|A16S|RequiresSNPConfirmation	LN:1477
@SQ	SN:MEG_11|Drugs|Aminoglycosides|Aminoglycoside-resistant_16S_ribosomal_subunit_protein|A16S|RequiresSNPConfirmation	LN:1441
@SQ	SN:MEG_12|Drugs|Aminoglycosides|Aminoglycoside-resistant_16S_ribosomal_subunit_protein|A16S|RequiresSNPConfirmation	LN:1454
@SQ	SN:MEG_399|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_RND_efflux_pumps|ACRA	LN:1194
@SQ	SN:MEG_400|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_RND_efflux_pumps|ACRA	LN:1197
@SQ	SN:MEG_401|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_RND_efflux_pumps|ACRB	LN:3150
@SQ	SN:MEG_402|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_RND_efflux_pumps|ACRB	LN:1194
@SQ	SN:MEG_411|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_RND_efflux_regulator|ACRR|RequiresSNPConfirmation	LN:648
@SQ	SN:MEG_412|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_RND_efflux_regulator|ACRR|RequiresSNPConfirmation	LN:540
@SQ	SN:MEG_413|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_RND_efflux_regulator|ACRR|RequiresSNPConfirmation	LN:651
@SQ	SN:MEG_414|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_RND_efflux_regulator|ACRR|RequiresSNPConfirmation	LN:651
@SQ	SN:MEG_1187|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_RND_efflux_pumps|AXYZ|RequiresSNPConfirmation	LN:639
@SQ	SN:MEG_1197|Drugs|Cationic_antimicrobial_peptides|Polymyxin_B_resistance_regulator|BASRS|RequiresSNPConfirmation	LN:666
@SQ	SN:MEG_1490|Drugs|Cationic_antimicrobial_peptides|Cationic_peptide-resistant_16S_ribosomal_subunit_protein|CAP16S|RequiresSNPConfirmation	LN:1542
@SQ	SN:MEG_1594|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1595|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1596|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1597|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1598|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1599|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1600|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1601|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1602|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1603|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1604|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1605|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1606|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1607|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:1260
@SQ	SN:MEG_1608|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1609|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1610|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1611|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1612|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1613|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1614|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1615|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1616|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1617|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1618|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1619|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1620|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:549
@SQ	SN:MEG_1621|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:630
@SQ	SN:MEG_1622|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:630
@SQ	SN:MEG_1623|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1624|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1625|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1626|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1627|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1628|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:630
@SQ	SN:MEG_1629|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1630|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATB|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_1631|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATP|RequiresSNPConfirmation	LN:639
@SQ	SN:MEG_1632|Drugs|Phenicol|Chloramphenicol_acetyltransferases|CATP|RequiresSNPConfirmation	LN:624
@SQ	SN:MEG_1642|Drugs|Lipopeptides|Daptomycin-resistant_mutant|CDSA|RequiresSNPConfirmation	LN:804
@SQ	SN:MEG_1730|Drugs|Lipopeptides|Daptomycin-resistant_mutant|CLS|RequiresSNPConfirmation	LN:1485
@SQ	SN:MEG_1731|Drugs|Lipopeptides|Daptomycin-resistant_mutant|CLS|RequiresSNPConfirmation	LN:1446
@SQ	SN:MEG_1732|Drugs|Lipopeptides|Daptomycin-resistant_mutant|CLS|RequiresSNPConfirmation	LN:1452
@SQ	SN:MEG_2572|Drugs|Trimethoprim|Dihydrofolate_reductase|DFRC|RequiresSNPConfirmation	LN:486
@SQ	SN:MEG_2573|Drugs|Trimethoprim|Dihydrofolate_reductase|DFRC|RequiresSNPConfirmation	LN:486
@SQ	SN:MEG_2615|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:474
@SQ	SN:MEG_2616|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:498
@SQ	SN:MEG_2617|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:561
@SQ	SN:MEG_2618|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:474
@SQ	SN:MEG_2619|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:474
@SQ	SN:MEG_2620|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:459
@SQ	SN:MEG_2621|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:474
@SQ	SN:MEG_2622|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:472
@SQ	SN:MEG_2623|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:474
@SQ	SN:MEG_2624|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:474
@SQ	SN:MEG_2625|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:474
@SQ	SN:MEG_2626|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:510
@SQ	SN:MEG_2627|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:558
@SQ	SN:MEG_2628|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:474
@SQ	SN:MEG_2629|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:459
@SQ	SN:MEG_2630|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:483
@SQ	SN:MEG_2631|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:498
@SQ	SN:MEG_2632|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:474
@SQ	SN:MEG_2633|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:474
@SQ	SN:MEG_2634|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:498
@SQ	SN:MEG_2635|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:474
@SQ	SN:MEG_2636|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:570
@SQ	SN:MEG_2637|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:474
@SQ	SN:MEG_2638|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:606
@SQ	SN:MEG_2639|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:474
@SQ	SN:MEG_2640|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:474
@SQ	SN:MEG_2641|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:522
@SQ	SN:MEG_2642|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:474
@SQ	SN:MEG_2643|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:474
@SQ	SN:MEG_2644|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:483
@SQ	SN:MEG_2645|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:474
@SQ	SN:MEG_2646|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:474
@SQ	SN:MEG_2647|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFRIII|RequiresSNPConfirmation	LN:489
@SQ	SN:MEG_2693|Drugs|Multi-drug_resistance|Multi-drug_ABC_efflux_pumps|EATAV|RequiresSNPConfirmation	LN:1503
@SQ	SN:MEG_2694|Drugs|Multi-drug_resistance|Multi-drug_ABC_efflux_pumps|EATAV|RequiresSNPConfirmation	LN:1503
@SQ	SN:MEG_2709|Drugs|Mycobacterium_tuberculosis-specific_Drug|Ethambutol_resistant_arabinosyltransferase|EMBA|RequiresSNPConfirmation	LN:3285
@SQ	SN:MEG_2710|Drugs|Mycobacterium_tuberculosis-specific_Drug|Ethambutol_resistant_arabinosyltransferase|EMBB|RequiresSNPConfirmation	LN:3297
@SQ	SN:MEG_2711|Drugs|Mycobacterium_tuberculosis-specific_Drug|Ethambutol_resistant_arabinosyltransferase|EMBB|RequiresSNPConfirmation	LN:3297
@SQ	SN:MEG_2712|Drugs|Mycobacterium_tuberculosis-specific_Drug|Ethambutol_resistant_arabinosyltransferase|EMBC|RequiresSNPConfirmation	LN:3285
@SQ	SN:MEG_2713|Drugs|Mycobacterium_tuberculosis-specific_Drug|Ethambutol_resistant_arabinosyltransferase|EMBR|RequiresSNPConfirmation	LN:1167
@SQ	SN:MEG_2866|Drugs|Mycobacterium_tuberculosis-specific_Drug|Ethionamide-resistant_mutant|ETHA|RequiresSNPConfirmation	LN:1470
@SQ	SN:MEG_2873|Biocides|Phenolic_compound_resistance|Triclosan-resistant_mutation|FABG|RequiresSNPConfirmation	LN:735
@SQ	SN:MEG_2933|Drugs|Mycobacterium_tuberculosis-specific_Drug|Para-aminosalicylic_acid_resistant_mutant|FOLC|RequiresSNPConfirmation	LN:1464
@SQ	SN:MEG_2934|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:801
@SQ	SN:MEG_2935|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:855
@SQ	SN:MEG_2936|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:822
@SQ	SN:MEG_2937|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:1143
@SQ	SN:MEG_2938|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:852
@SQ	SN:MEG_2939|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:843
@SQ	SN:MEG_2940|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:1188
@SQ	SN:MEG_2941|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:912
@SQ	SN:MEG_2942|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:792
@SQ	SN:MEG_2943|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:876
@SQ	SN:MEG_2944|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:885
@SQ	SN:MEG_2945|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:804
@SQ	SN:MEG_2946|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:843
@SQ	SN:MEG_2947|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:795
@SQ	SN:MEG_2948|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:849
@SQ	SN:MEG_2949|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:804
@SQ	SN:MEG_2950|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:804
@SQ	SN:MEG_2951|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:804
@SQ	SN:MEG_2952|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:849
@SQ	SN:MEG_2953|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:858
@SQ	SN:MEG_2954|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:840
@SQ	SN:MEG_2955|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:894
@SQ	SN:MEG_2956|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:804
@SQ	SN:MEG_2957|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:978
@SQ	SN:MEG_2958|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:852
@SQ	SN:MEG_2959|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:849
@SQ	SN:MEG_2960|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:1143
@SQ	SN:MEG_2961|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:852
@SQ	SN:MEG_2962|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:804
@SQ	SN:MEG_2963|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:804
@SQ	SN:MEG_2964|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:927
@SQ	SN:MEG_2965|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:843
@SQ	SN:MEG_2966|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:849
@SQ	SN:MEG_2967|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:858
@SQ	SN:MEG_2968|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:804
@SQ	SN:MEG_2969|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:957
@SQ	SN:MEG_2970|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:819
@SQ	SN:MEG_2971|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:843
@SQ	SN:MEG_2972|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:849
@SQ	SN:MEG_3065|Drugs|Fusidic_acid|Fusidic_acid-resistant_mutation|FUSA|RequiresSNPConfirmation	LN:2082
@SQ	SN:MEG_3071|Drugs|Fusidic_acid|Fusidic_acid-resistant_mutation|FUSE|RequiresSNPConfirmation	LN:537
@SQ	SN:MEG_3135|Drugs|Aminoglycosides|Aminoglycoside-resistant_gidB|GIDB|RequiresSNPConfirmation	LN:675
@SQ	SN:MEG_3143|Drugs|Fosfomycin|Fosfomycin_target_mutation|GLPT|RequiresSNPConfirmation	LN:1359
@SQ	SN:MEG_3144|Drugs|Fosfomycin|Fosfomycin_target_mutation|GLPT|RequiresSNPConfirmation	LN:1359
@SQ	SN:MEG_3175|Drugs|Lipopeptides|Daptomycin-resistant_mutant|GSHF|RequiresSNPConfirmation	LN:2271
@SQ	SN:MEG_3176|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2511
@SQ	SN:MEG_3177|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2751
@SQ	SN:MEG_3178|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2628
@SQ	SN:MEG_3179|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2787
@SQ	SN:MEG_3180|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2517
@SQ	SN:MEG_3181|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:3750
@SQ	SN:MEG_3182|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2676
@SQ	SN:MEG_3183|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2715
@SQ	SN:MEG_3184|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2646
@SQ	SN:MEG_3185|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2511
@SQ	SN:MEG_3186|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2661
@SQ	SN:MEG_3187|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2628
@SQ	SN:MEG_3188|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2634
@SQ	SN:MEG_3189|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2634
@SQ	SN:MEG_3190|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2589
@SQ	SN:MEG_3191|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2670
@SQ	SN:MEG_3192|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2715
@SQ	SN:MEG_3193|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2715
@SQ	SN:MEG_3194|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2751
@SQ	SN:MEG_3195|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2535
@SQ	SN:MEG_3196|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2541
@SQ	SN:MEG_3197|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2442
@SQ	SN:MEG_3198|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2469
@SQ	SN:MEG_3199|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2682
@SQ	SN:MEG_3200|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2502
@SQ	SN:MEG_3201|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2637
@SQ	SN:MEG_3202|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2592
@SQ	SN:MEG_3203|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2643
@SQ	SN:MEG_3204|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2487
@SQ	SN:MEG_3205|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2637
@SQ	SN:MEG_3206|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2595
@SQ	SN:MEG_3207|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2664
@SQ	SN:MEG_3208|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2634
@SQ	SN:MEG_3209|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2592
@SQ	SN:MEG_3210|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2424
@SQ	SN:MEG_3211|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2382
@SQ	SN:MEG_3212|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2628
@SQ	SN:MEG_3213|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2493
@SQ	SN:MEG_3214|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2529
@SQ	SN:MEG_3215|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2715
@SQ	SN:MEG_3216|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2685
@SQ	SN:MEG_3217|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2751
@SQ	SN:MEG_3218|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2592
@SQ	SN:MEG_3219|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2427
@SQ	SN:MEG_3220|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2637
@SQ	SN:MEG_3221|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2661
@SQ	SN:MEG_3222|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2664
@SQ	SN:MEG_3223|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2664
@SQ	SN:MEG_3224|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2676
@SQ	SN:MEG_3225|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2469
@SQ	SN:MEG_3226|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2772
@SQ	SN:MEG_3227|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2628
@SQ	SN:MEG_3228|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2715
@SQ	SN:MEG_3229|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2682
@SQ	SN:MEG_3230|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2664
@SQ	SN:MEG_3231|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2682
@SQ	SN:MEG_3232|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2715
@SQ	SN:MEG_3233|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2670
@SQ	SN:MEG_3234|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2637
@SQ	SN:MEG_3235|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2466
@SQ	SN:MEG_3236|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2670
@SQ	SN:MEG_3237|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2517
@SQ	SN:MEG_3238|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRB|RequiresSNPConfirmation	LN:2415
@SQ	SN:MEG_3239|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRB|RequiresSNPConfirmation	LN:2037
@SQ	SN:MEG_3240|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRB|RequiresSNPConfirmation	LN:1953
@SQ	SN:MEG_3241|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRB|RequiresSNPConfirmation	LN:1914
@SQ	SN:MEG_3242|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRB|RequiresSNPConfirmation	LN:2415
@SQ	SN:MEG_3243|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRB|RequiresSNPConfirmation	LN:2028
@SQ	SN:MEG_3244|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|GYRBA|RequiresSNPConfirmation	LN:2415
@SQ	SN:MEG_3245|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|GYRBA|RequiresSNPConfirmation	LN:1932
@SQ	SN:MEG_3246|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|GYRBA|RequiresSNPConfirmation	LN:2430
@SQ	SN:MEG_3247|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|GYRBA|RequiresSNPConfirmation	LN:1718
@SQ	SN:MEG_3248|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|GYRBA|RequiresSNPConfirmation	LN:2034
@SQ	SN:MEG_3249|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|GYRBA|RequiresSNPConfirmation	LN:2034
@SQ	SN:MEG_3250|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|GYRBA|RequiresSNPConfirmation	LN:2034
@SQ	SN:MEG_3251|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRC|RequiresSNPConfirmation	LN:2256
@SQ	SN:MEG_3296|Drugs|Mupirocin|Mupirocin-resistant_isoleucyl-tRNA_synthetase|ILES|RequiresSNPConfirmation	LN:2754
@SQ	SN:MEG_3429|Drugs|Mycobacterium_tuberculosis-specific_Drug|Isoniazid-resistant_mutant|INHA|RequiresSNPConfirmation	LN:810
@SQ	SN:MEG_3430|Drugs|Mycobacterium_tuberculosis-specific_Drug|Ethambutol-resistant_mutant|INIA|RequiresSNPConfirmation	LN:1923
@SQ	SN:MEG_3431|Drugs|Mycobacterium_tuberculosis-specific_Drug|Ethambutol-resistant_mutant|INIB|RequiresSNPConfirmation	LN:1440
@SQ	SN:MEG_3432|Drugs|Mycobacterium_tuberculosis-specific_Drug|Ethambutol-resistant_mutant|INIC|RequiresSNPConfirmation	LN:1482
@SQ	SN:MEG_3445|Drugs|Mycobacterium_tuberculosis-specific_Drug|Isoniazid-resistant_mutant|KASA|RequiresSNPConfirmation	LN:1251
@SQ	SN:MEG_3446|Drugs|Mycobacterium_tuberculosis-specific_Drug|Isoniazid-resistant_mutant|KATG|RequiresSNPConfirmation	LN:2223
@SQ	SN:MEG_3586|Drugs|Lipopeptides|Daptomycin-resistant_mutant|LIAFSR|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_3587|Drugs|Lipopeptides|Daptomycin-resistant_mutant|LIAFSR|RequiresSNPConfirmation	LN:1104
@SQ	SN:MEG_3588|Drugs|Lipopeptides|Daptomycin-resistant_mutant|LIAFSR|RequiresSNPConfirmation	LN:732
@SQ	SN:MEG_3589|Drugs|Lipopeptides|Daptomycin-resistant_mutant|LIAFSR|RequiresSNPConfirmation	LN:1068
@SQ	SN:MEG_3590|Drugs|Lipopeptides|Daptomycin-resistant_mutant|LIAFSR|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_3591|Drugs|Lipopeptides|Daptomycin-resistant_mutant|LIAFSR|RequiresSNPConfirmation	LN:732
@SQ	SN:MEG_3594|Drugs|Multi-drug_resistance|Multi-drug_ABC_efflux_pumps|LMRA|RequiresSNPConfirmation	LN:1446
@SQ	SN:MEG_3626|Drugs|Lipopeptides|Colistin-resistant_mutant|LPXA|RequiresSNPConfirmation	LN:789
@SQ	SN:MEG_3627|Drugs|Lipopeptides|Colistin-resistant_mutant|LPXC|RequiresSNPConfirmation	LN:903
@SQ	SN:MEG_3661|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_RND_efflux_regulator|MARA	LN:384
@SQ	SN:MEG_3662|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_RND_efflux_regulator|MARA	LN:390
@SQ	SN:MEG_3663|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_RND_efflux_regulator|MARR	LN:435
@SQ	SN:MEG_3740|Drugs|Multi-drug_resistance|MDR_23S_rRNA_mutation|MDR23S|RequiresSNPConfirmation	LN:2734
@SQ	SN:MEG_3831|Drugs|Lipopeptides|Lysocin-resistant_mutant|MENA|RequiresSNPConfirmation	LN:939
@SQ	SN:MEG_3931|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_RND_efflux_regulator|MEXS|RequiresSNPConfirmation	LN:1020
@SQ	SN:MEG_3932|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_RND_efflux_regulator|MEXT|RequiresSNPConfirmation	LN:915
@SQ	SN:MEG_3933|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_RND_efflux_regulator|MEXT|RequiresSNPConfirmation	LN:1044
@SQ	SN:MEG_3941|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_RND_efflux_regulator|MEXZ|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_3975|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation	LN:3162
@SQ	SN:MEG_3976|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation	LN:2975
@SQ	SN:MEG_3977|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation	LN:2904
@SQ	SN:MEG_3978|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation	LN:3156
@SQ	SN:MEG_3979|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation	LN:2884
@SQ	SN:MEG_3980|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation	LN:3113
@SQ	SN:MEG_3981|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation	LN:3403
@SQ	SN:MEG_3982|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation	LN:3112
@SQ	SN:MEG_3983|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation	LN:2886
@SQ	SN:MEG_3984|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation	LN:3112
@SQ	SN:MEG_3985|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation	LN:3103
@SQ	SN:MEG_3986|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation	LN:2940
@SQ	SN:MEG_3987|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation	LN:2912
@SQ	SN:MEG_3988|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation	LN:2900
@SQ	SN:MEG_3989|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation	LN:2905
@SQ	SN:MEG_3990|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation	LN:3149
@SQ	SN:MEG_3991|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation	LN:2996
@SQ	SN:MEG_3992|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation	LN:2904
@SQ	SN:MEG_3993|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation	LN:3118
@SQ	SN:MEG_3994|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation	LN:3135
@SQ	SN:MEG_4057|Drugs|Lipopeptides|Daptomycin-resistant_mutant|MPRFD|RequiresSNPConfirmation	LN:2511
@SQ	SN:MEG_4087|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_MFS_efflux_regulator|MTRR|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_4088|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_MFS_efflux_regulator|MTRR|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_4092|Drugs|Fosfomycin|Fosfomycin_target_mutation|MURA|RequiresSNPConfirmation	LN:1284
@SQ	SN:MEG_4093|Drugs|Fosfomycin|Fosfomycin_target_mutation|MURA|RequiresSNPConfirmation	LN:1335
@SQ	SN:MEG_4094|Drugs|Fosfomycin|Fosfomycin_target_mutation|MURA|RequiresSNPConfirmation	LN:1266
@SQ	SN:MEG_4095|Drugs|Fosfomycin|Fosfomycin_target_mutation|MURA|RequiresSNPConfirmation	LN:1260
@SQ	SN:MEG_4096|Drugs|Fosfomycin|Fosfomycin_target_mutation|MURA|RequiresSNPConfirmation	LN:1257
@SQ	SN:MEG_4097|Drugs|Glycopeptides|Vancomycin-resistant_mutation|MURG|RequiresSNPConfirmation	LN:1080
@SQ	SN:MEG_4109|Drugs|Multi-drug_resistance|Multi-drug_RND_efflux_pumps|NALC|RequiresSNPConfirmation	LN:642
@SQ	SN:MEG_4110|Drugs|Multi-drug_resistance|Multi-drug_RND_efflux_regulator|NALD|RequiresSNPConfirmation	LN:639
@SQ	SN:MEG_4130|Drugs|Mycobacterium_tuberculosis-specific_Drug|Isoniazid-resistant_mutant|NDH|RequiresSNPConfirmation	LN:1392
@SQ	SN:MEG_4131|Drugs|Mycobacterium_tuberculosis-specific_Drug|Isoniazid-resistant_mutant|NDH|RequiresSNPConfirmation	LN:1392
@SQ	SN:MEG_4132|Drugs|Mycobacterium_tuberculosis-specific_Drug|Isoniazid-resistant_mutant|NDH|RequiresSNPConfirmation	LN:1374
@SQ	SN:MEG_4161|Drugs|Multi-drug_resistance|Multi-drug_RND_efflux_pumps|NFXB|RequiresSNPConfirmation	LN:564
@SQ	SN:MEG_4223|Drugs|Oxazolidinone|Oxazolidinone-resistant_23S_rRNA_mutation|O23S|RequiresSNPConfirmation	LN:2926
@SQ	SN:MEG_4279|Drugs|betalactams|Mutant_porin_proteins|OMP36|RequiresSNPConfirmation	LN:1128
@SQ	SN:MEG_4280|Drugs|betalactams|Mutant_porin_proteins|OMP36|RequiresSNPConfirmation	LN:1131
@SQ	SN:MEG_4281|Drugs|betalactams|Mutant_porin_proteins|OMP36|RequiresSNPConfirmation	LN:1128
@SQ	SN:MEG_4282|Drugs|betalactams|Mutant_porin_proteins|OMP36|RequiresSNPConfirmation	LN:1128
@SQ	SN:MEG_4286|Drugs|Multi-drug_resistance|MDR_mutant_porin_proteins|OMPF|RequiresSNPConfirmation	LN:927
@SQ	SN:MEG_4287|Drugs|Multi-drug_resistance|MDR_mutant_porin_proteins|OMPF|RequiresSNPConfirmation	LN:1098
@SQ	SN:MEG_4288|Drugs|Multi-drug_resistance|MDR_mutant_porin_proteins|OMPF|RequiresSNPConfirmation	LN:1092
@SQ	SN:MEG_4289|Drugs|betalactams|Mutant_porin_proteins|OMPFB|RequiresSNPConfirmation	LN:1089
@SQ	SN:MEG_4296|Drugs|betalactams|Mutant_porin_proteins|OPRD|RequiresSNPConfirmation	LN:1332
@SQ	SN:MEG_5322|Drugs|Pactamycin|Pactamycin-resistant_16S_ribosomal_subunit_protein|P16S|RequiresSNPConfirmation	LN:1473
@SQ	SN:MEG_5323|Drugs|Pleuromutilin|Pleuromutilin-resistant_23S_rRNA_mutation|P23S|RequiresSNPConfirmation	LN:2913
@SQ	SN:MEG_5325|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2265
@SQ	SN:MEG_5326|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2304
@SQ	SN:MEG_5327|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2472
@SQ	SN:MEG_5328|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2553
@SQ	SN:MEG_5329|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2403
@SQ	SN:MEG_5330|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2802
@SQ	SN:MEG_5331|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2346
@SQ	SN:MEG_5332|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2259
@SQ	SN:MEG_5333|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2259
@SQ	SN:MEG_5334|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2259
@SQ	SN:MEG_5335|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2220
@SQ	SN:MEG_5336|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2220
@SQ	SN:MEG_5337|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2256
@SQ	SN:MEG_5338|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2544
@SQ	SN:MEG_5339|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2871
@SQ	SN:MEG_5340|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2481
@SQ	SN:MEG_5341|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2463
@SQ	SN:MEG_5342|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2271
@SQ	SN:MEG_5343|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2451
@SQ	SN:MEG_5344|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2403
@SQ	SN:MEG_5345|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2403
@SQ	SN:MEG_5346|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2958
@SQ	SN:MEG_5347|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2256
@SQ	SN:MEG_5348|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2259
@SQ	SN:MEG_5349|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2424
@SQ	SN:MEG_5350|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2220
@SQ	SN:MEG_5351|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2304
@SQ	SN:MEG_5352|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2220
@SQ	SN:MEG_5353|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2259
@SQ	SN:MEG_5354|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2403
@SQ	SN:MEG_5355|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2403
@SQ	SN:MEG_5356|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2481
@SQ	SN:MEG_5357|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2265
@SQ	SN:MEG_5358|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2220
@SQ	SN:MEG_5359|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2403
@SQ	SN:MEG_5360|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2316
@SQ	SN:MEG_5361|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2259
@SQ	SN:MEG_5362|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2403
@SQ	SN:MEG_5363|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2193
@SQ	SN:MEG_5364|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1974
@SQ	SN:MEG_5365|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1896
@SQ	SN:MEG_5366|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1998
@SQ	SN:MEG_5367|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1884
@SQ	SN:MEG_5368|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1887
@SQ	SN:MEG_5369|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1881
@SQ	SN:MEG_5370|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1971
@SQ	SN:MEG_5371|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1956
@SQ	SN:MEG_5372|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1944
@SQ	SN:MEG_5373|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:2058
@SQ	SN:MEG_5374|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1893
@SQ	SN:MEG_5375|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1896
@SQ	SN:MEG_5376|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1881
@SQ	SN:MEG_5377|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1968
@SQ	SN:MEG_5378|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1971
@SQ	SN:MEG_5379|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1998
@SQ	SN:MEG_5380|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1986
@SQ	SN:MEG_5381|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1893
@SQ	SN:MEG_5382|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1893
@SQ	SN:MEG_5383|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1998
@SQ	SN:MEG_5384|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1998
@SQ	SN:MEG_5385|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1944
@SQ	SN:MEG_5386|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1890
@SQ	SN:MEG_5387|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1998
@SQ	SN:MEG_5388|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1884
@SQ	SN:MEG_5389|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1962
@SQ	SN:MEG_5390|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1998
@SQ	SN:MEG_5391|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1992
@SQ	SN:MEG_5392|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:2079
@SQ	SN:MEG_5393|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1893
@SQ	SN:MEG_5394|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PAREF|RequiresSNPConfirmation	LN:1893
@SQ	SN:MEG_5395|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PAREF|RequiresSNPConfirmation	LN:1998
@SQ	SN:MEG_5399|Drugs|betalactams|Penicillin_binding_protein|PBP1A|RequiresSNPConfirmation	LN:2160
@SQ	SN:MEG_5401|Drugs|betalactams|Penicillin_binding_protein|PBP2|RequiresSNPConfirmation	LN:1746
@SQ	SN:MEG_5405|Drugs|betalactams|Penicillin_binding_protein|PBP2|RequiresSNPConfirmation	LN:1746
@SQ	SN:MEG_5406|Drugs|betalactams|Penicillin_binding_protein|PBP2B|RequiresSNPConfirmation	LN:2058
@SQ	SN:MEG_5407|Drugs|betalactams|Penicillin_binding_protein|PBP2X|RequiresSNPConfirmation	LN:2253
@SQ	SN:MEG_5408|Drugs|betalactams|Penicillin_binding_protein|PBP3|RequiresSNPConfirmation	LN:1833
@SQ	SN:MEG_5779|Drugs|Lipopeptides|Daptomycin-resistant_mutant|PGSA|RequiresSNPConfirmation	LN:1143
@SQ	SN:MEG_5780|Drugs|Lipopeptides|Daptomycin-resistant_mutant|PGSA|RequiresSNPConfirmation	LN:579
@SQ	SN:MEG_5781|Drugs|Phenicol|Phenicol-resistant_23S_rRNA_mutation|PH23S|RequiresSNPConfirmation	LN:2905
@SQ	SN:MEG_5782|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_ABC_efflux_regulator|PHOB|RequiresSNPConfirmation	LN:714
@SQ	SN:MEG_5783|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_ABC_efflux_regulator|PHOB|RequiresSNPConfirmation	LN:672
@SQ	SN:MEG_5784|Drugs|Lipopeptides|Colistin-resistant_mutant|PHOP|RequiresSNPConfirmation	LN:678
@SQ	SN:MEG_5785|Drugs|Lipopeptides|Colistin-resistant_mutant|PHOQ|RequiresSNPConfirmation	LN:1347
@SQ	SN:MEG_5803|Drugs|Mycobacterium_tuberculosis-specific_Drug|Pyrazinamide-resistant_mutant|PNCA|RequiresSNPConfirmation	LN:561
@SQ	SN:MEG_5808|Drugs|Multi-drug_resistance|MDR_mutant_porin_proteins|POR|RequiresSNPConfirmation	LN:1047
@SQ	SN:MEG_5819|Drugs|Fosfomycin|Fosfomycin_target_mutation|PTSL|RequiresSNPConfirmation	LN:1728
@SQ	SN:MEG_6044|Drugs|Multi-drug_resistance|Multi-drug_RND_efflux_regulator|RAMA	LN:375
@SQ	SN:MEG_6045|Drugs|Multi-drug_resistance|Multi-drug_RND_efflux_regulator|RAMR|RequiresSNPConfirmation	LN:582
@SQ	SN:MEG_6046|Drugs|Multi-drug_resistance|Multi-drug_RND_efflux_regulator|RAMR|RequiresSNPConfirmation	LN:582
@SQ	SN:MEG_6055|Drugs|Mycobacterium_tuberculosis-specific_Drug|Para-aminosalicylic_acid_resistant_mutant|RIBD|RequiresSNPConfirmation	LN:777
@SQ	SN:MEG_6082|Multi-compound|Drug_and_biocide_and_metal_resistance|Drug_and_biocide_and_metal_RND_efflux_pumps|ROBA	LN:870
@SQ	SN:MEG_6090|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:3519
@SQ	SN:MEG_6091|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:3537
@SQ	SN:MEG_6092|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:3489
@SQ	SN:MEG_6093|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:3498
@SQ	SN:MEG_6094|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:4029
@SQ	SN:MEG_6095|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:3729
@SQ	SN:MEG_6096|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:4134
@SQ	SN:MEG_6097|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:4011
@SQ	SN:MEG_6098|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:3552
@SQ	SN:MEG_6099|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:4062
@SQ	SN:MEG_6100|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:4089
@SQ	SN:MEG_6101|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:4107
@SQ	SN:MEG_6102|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:3540
@SQ	SN:MEG_6103|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:3699
@SQ	SN:MEG_6104|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:4110
@SQ	SN:MEG_6105|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:3651
@SQ	SN:MEG_6106|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:3624
@SQ	SN:MEG_6107|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:4029
@SQ	SN:MEG_6108|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:4071
@SQ	SN:MEG_6109|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:4185
@SQ	SN:MEG_6110|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:3552
@SQ	SN:MEG_6111|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:3552
@SQ	SN:MEG_6112|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:3552
@SQ	SN:MEG_6113|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:4029
@SQ	SN:MEG_6114|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:4107
@SQ	SN:MEG_6115|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:3591
@SQ	SN:MEG_6116|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:3543
@SQ	SN:MEG_6117|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:3552
@SQ	SN:MEG_6118|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:3627
@SQ	SN:MEG_6119|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:4137
@SQ	SN:MEG_6120|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:3717
@SQ	SN:MEG_6121|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:4089
@SQ	SN:MEG_6122|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:4029
@SQ	SN:MEG_6123|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:4029
@SQ	SN:MEG_6124|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:3552
@SQ	SN:MEG_6125|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:4179
@SQ	SN:MEG_6126|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:3612
@SQ	SN:MEG_6127|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:4074
@SQ	SN:MEG_6128|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:3552
@SQ	SN:MEG_6129|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:4113
@SQ	SN:MEG_6130|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:4029
@SQ	SN:MEG_6131|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:4089
@SQ	SN:MEG_6132|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:3582
@SQ	SN:MEG_6133|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:3552
@SQ	SN:MEG_6134|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:3537
@SQ	SN:MEG_6135|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:4029
@SQ	SN:MEG_6136|Drugs|Multi-drug_resistance|MDR_mutation|MRPOB|RequiresSNPConfirmation	LN:3552
@SQ	SN:MEG_6138|Drugs|Glycopeptides|Vancomycin-resistant_mutation|RPOC|RequiresSNPConfirmation	LN:3486
@SQ	SN:MEG_6139|Drugs|Lipopeptides|Daptomycin-resistant_beta-subunit_of_RNA_polymerase_RpoC|RPOCL|RequiresSNPConfirmation	LN:3624
@SQ	SN:MEG_6142|Drugs|Mycobacterium_tuberculosis-specific_Drug|Pyrazinamide-resistant_mutant|RPSA|RequiresSNPConfirmation	LN:1446
@SQ	SN:MEG_6143|Drugs|Tetracyclines|Tetracycline_resistance_ribosomal_protection_proteins|RPSJ|RequiresSNPConfirmation	LN:312
@SQ	SN:MEG_6144|Drugs|Aminoglycosides|Aminoglycoside-resistant_30S_ribosomal_protein_S12|RPSL|RequiresSNPConfirmation	LN:375
@SQ	SN:MEG_6145|Drugs|Aminoglycosides|Aminoglycoside-resistant_16S_ribosomal_subunit_protein|RRSA|RequiresSNPConfirmation	LN:1529
@SQ	SN:MEG_6146|Drugs|Aminoglycosides|Aminoglycoside-resistant_16S_ribosomal_subunit_protein|RRSC|RequiresSNPConfirmation	LN:1542
@SQ	SN:MEG_6147|Drugs|Aminoglycosides|Aminoglycoside-resistant_16S_ribosomal_subunit_protein|RRSH|RequiresSNPConfirmation	LN:1542
@SQ	SN:MEG_6548|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_ABC_efflux_regulator|SOXR|RequiresSNPConfirmation	LN:459
@SQ	SN:MEG_6549|Biocides|Multi-biocide_resistance|Multi-biocide_resistance_regulator|SOXRB	LN:465
@SQ	SN:MEG_6550|Biocides|Multi-biocide_resistance|Multi-biocide_resistance_regulator|SOXRB	LN:465
@SQ	SN:MEG_6551|Multi-compound|Drug_and_biocide_and_metal_resistance|Drug_and_biocide_and_metal_resistance_regulator|SOXS	LN:324
@SQ	SN:MEG_6552|Multi-compound|Drug_and_biocide_and_metal_resistance|Drug_and_biocide_and_metal_resistance_regulator|SOXS|RequiresSNPConfirmation	LN:324
@SQ	SN:MEG_6963|Drugs|Tetracyclines|Tetracycline-resistant_16S_ribosomal_subunit_protein|TET16S|RequiresSNPConfirmation	LN:1501
@SQ	SN:MEG_6964|Drugs|Tetracyclines|Tetracycline-resistant_16S_ribosomal_subunit_protein|TET16S|RequiresSNPConfirmation	LN:1486
@SQ	SN:MEG_7185|Drugs|Tetracyclines|Tetracycline_transcriptional_repressor|TETR|RequiresSNPConfirmation	LN:651
@SQ	SN:MEG_7186|Drugs|Tetracyclines|Tetracycline_transcriptional_repressor|TETR|RequiresSNPConfirmation	LN:633
@SQ	SN:MEG_7187|Drugs|Tetracyclines|Tetracycline_transcriptional_repressor|TETR|RequiresSNPConfirmation	LN:642
@SQ	SN:MEG_7188|Drugs|Tetracyclines|Tetracycline_transcriptional_repressor|TETR|RequiresSNPConfirmation	LN:648
@SQ	SN:MEG_7189|Drugs|Tetracyclines|Tetracycline_transcriptional_repressor|TETR|RequiresSNPConfirmation	LN:612
@SQ	SN:MEG_7190|Drugs|Tetracyclines|Tetracycline_transcriptional_repressor|TETR|RequiresSNPConfirmation	LN:615
@SQ	SN:MEG_7191|Drugs|Tetracyclines|Tetracycline_transcriptional_repressor|TETR|RequiresSNPConfirmation	LN:600
@SQ	SN:MEG_7192|Drugs|Tetracyclines|Tetracycline_transcriptional_repressor|TETR|RequiresSNPConfirmation	LN:732
@SQ	SN:MEG_7193|Drugs|Tetracyclines|Tetracycline_transcriptional_repressor|TETR|RequiresSNPConfirmation	LN:627
@SQ	SN:MEG_7194|Drugs|Tetracyclines|Tetracycline_transcriptional_repressor|TETR|RequiresSNPConfirmation	LN:624
@SQ	SN:MEG_7195|Drugs|Tetracyclines|Tetracycline_transcriptional_repressor|TETR|RequiresSNPConfirmation	LN:642
@SQ	SN:MEG_7196|Drugs|Tetracyclines|Tetracycline_resistance_MFS_efflux_regulator|TETRM|RequiresSNPConfirmation	LN:627
@SQ	SN:MEG_7197|Drugs|Tetracyclines|Tetracycline_resistance_MFS_efflux_regulator|TETRM|RequiresSNPConfirmation	LN:567
@SQ	SN:MEG_7250|Drugs|Mycobacterium_tuberculosis-specific_Drug|Para-aminosalicylic_acid_resistant_mutant|THYA|RequiresSNPConfirmation	LN:792
@SQ	SN:MEG_7258|Drugs|Aminoglycosides|Aminoglycoside-resistant_arabinosyltransferase|TLYA|RequiresSNPConfirmation	LN:807
@SQ	SN:MEG_7259|Drugs|Rifampin|Rifamycin-resistant_arabinosyltransferase|REMB|RequiresSNPConfirmation	LN:3297
@SQ	SN:MEG_7301|Drugs|Elfamycins|EF-Tu_inhibition|TUFAB|RequiresSNPConfirmation	LN:1230
@SQ	SN:MEG_7302|Drugs|Elfamycins|EF-Tu_inhibition|TUFAB|RequiresSNPConfirmation	LN:1185
@SQ	SN:MEG_7303|Drugs|Elfamycins|EF-Tu_inhibition|TUFAB|RequiresSNPConfirmation	LN:1272
@SQ	SN:MEG_7304|Drugs|Elfamycins|EF-Tu_inhibition|TUFAB|RequiresSNPConfirmation	LN:1194
@SQ	SN:MEG_7305|Drugs|Elfamycins|EF-Tu_inhibition|TUFAB|RequiresSNPConfirmation	LN:1194
@SQ	SN:MEG_7306|Drugs|Elfamycins|EF-Tu_inhibition|TUFAB|RequiresSNPConfirmation	LN:1194
@SQ	SN:MEG_7307|Drugs|Elfamycins|EF-Tu_inhibition|TUFAB|RequiresSNPConfirmation	LN:705
@SQ	SN:MEG_7308|Drugs|Elfamycins|EF-Tu_inhibition|TUFAB|RequiresSNPConfirmation	LN:1191
@SQ	SN:MEG_7309|Drugs|Elfamycins|EF-Tu_inhibition|TUFAB|RequiresSNPConfirmation	LN:924
@SQ	SN:MEG_7310|Drugs|Elfamycins|EF-Tu_inhibition|TUFAB|RequiresSNPConfirmation	LN:1185
@SQ	SN:MEG_7311|Drugs|Elfamycins|EF-Tu_inhibition|TUFAB|RequiresSNPConfirmation	LN:1200
@SQ	SN:MEG_7312|Drugs|Elfamycins|EF-Tu_inhibition|TUFAB|RequiresSNPConfirmation	LN:1185
@SQ	SN:MEG_7313|Drugs|Elfamycins|EF-Tu_inhibition|TUFAB|RequiresSNPConfirmation	LN:1185
@SQ	SN:MEG_7314|Drugs|Elfamycins|EF-Tu_inhibition|TUFAB|RequiresSNPConfirmation	LN:1194
@SQ	SN:MEG_7315|Drugs|Elfamycins|EF-Tu_inhibition|TUFAB|RequiresSNPConfirmation	LN:1185
@SQ	SN:MEG_7316|Drugs|Elfamycins|EF-Tu_inhibition|TUFAB|RequiresSNPConfirmation	LN:1236
@SQ	SN:MEG_7317|Drugs|Elfamycins|EF-Tu_inhibition|TUFAB|RequiresSNPConfirmation	LN:1185
@SQ	SN:MEG_7318|Drugs|Elfamycins|EF-Tu_inhibition|TUFAB|RequiresSNPConfirmation	LN:1188
@SQ	SN:MEG_7319|Drugs|Elfamycins|EF-Tu_inhibition|TUFAB|RequiresSNPConfirmation	LN:1185
@SQ	SN:MEG_7320|Drugs|Elfamycins|EF-Tu_inhibition|TUFAB|RequiresSNPConfirmation	LN:1194
@SQ	SN:MEG_7321|Drugs|Elfamycins|EF-Tu_inhibition|TUFAB|RequiresSNPConfirmation	LN:1185
@SQ	SN:MEG_7329|Drugs|Fosfomycin|Fosfomycin_target_mutation|UHPT|RequiresSNPConfirmation	LN:1380
@SQ	SN:MEG_7330|Drugs|Fosfomycin|Fosfomycin_target_mutation|UHPT|RequiresSNPConfirmation	LN:1392
@SQ	SN:MEG_7333|Drugs|Glycopeptides|VanA-type_resistance_protein|VAN|RequiresSNPConfirmation	LN:1056
@SQ	SN:MEG_7806|Drugs|Lipopeptides|Daptomycin-resistant_mutant|WALK|RequiresSNPConfirmation	LN:1827
@SQ	SN:MEG_7840|Drugs|Lipopeptides|Defensin-resistant_mutant|YKKCL|RequiresSNPConfirmation	LN:2571
@SQ	SN:MEG_7843|Drugs|Lipopeptides|Daptomycin-resistant_mutant|YYBT|RequiresSNPConfirmation	LN:1977
@SQ	SN:MEG_8147|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2568
@SQ	SN:MEG_8677|Drugs|Aminoglycosides|Aminoglycoside-resistant_16S_ribosomal_subunit_protein|RPSL|RequiresSNPConfirmation	LN:387
@SQ	SN:MEG_8145|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2469
@SQ	SN:MEG_8461|Drugs|betalactams|Penicillin_binding_protein|PBP5|RequiresSNPConfirmation	LN:2037
@SQ	SN:MEG_8249|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation	LN:2904
@SQ	SN:MEG_7974|Drugs|betalactams|Class_C_betalactamase_regulator|AMPCR|RequiresSNPConfirmation	LN:145
@SQ	SN:MEG_8445|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2259
@SQ	SN:MEG_8634|Drugs|Lipopeptides|Colistin-resistant_mutant|PMRAR|RequiresSNPConfirmation	LN:669
@SQ	SN:MEG_8672|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:8673
@SQ	SN:MEG_7884|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_RND_efflux_regulator|ACRR|RequiresSNPConfirmation	LN:651
@SQ	SN:MEG_8144|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:1005
@SQ	SN:MEG_8151|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRB|RequiresSNPConfirmation	LN:2418
@SQ	SN:MEG_8289|Drugs|Multi-drug_resistance|MDR_mutant_porin_proteins|OMPK|RequiresSNPConfirmation	LN:1080
@SQ	SN:MEG_8290|Drugs|betalactams|Mutant_porin_proteins|OMPK36|RequiresSNPConfirmation	LN:1092
@SQ	SN:MEG_8447|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:522
@SQ	SN:MEG_8659|Drugs|Multi-drug_resistance|Multi-drug_RND_efflux_regulator|RAMR|RequiresSNPConfirmation	LN:585
@SQ	SN:MEG_8678|Drugs|Aminoglycosides|Aminoglycoside-resistant_16S_ribosomal_subunit_protein|RPSL|RequiresSNPConfirmation	LN:387
@SQ	SN:MEG_8687|Drugs|Mycobacterium_tuberculosis-specific_Drug|MDR_mutant_promoter|RV0678|RequiresSNPConfirmation	LN:498
@SQ	SN:MEG_7971|Drugs|Mycobacterium_tuberculosis-specific_Drug|Isoniazid-resistant_mutant|AHPC|RequiresSNPConfirmation	LN:768
@SQ	SN:MEG_7973|Drugs|Mycobacterium_tuberculosis-specific_Drug|Cycloserine-resistant_mutant|ALR|RequiresSNPConfirmation	LN:1227
@SQ	SN:MEG_8078|Drugs|Aminoglycosides|Aminoglycoside-resistant_eis_promoter|EIS|RequiresSNPConfirmation	LN:1328
@SQ	SN:MEG_8079|Drugs|Mycobacterium_tuberculosis-specific_Drug|Ethambutol_resistant_arabinosyltransferase|EMBA|RequiresSNPConfirmation	LN:3400
@SQ	SN:MEG_8084|Drugs|Mycobacterium_tuberculosis-specific_Drug|Ethionamide-resistant_mutant|ETHA|RequiresSNPConfirmation	LN:1544
@SQ	SN:MEG_8085|Drugs|Mycobacterium_tuberculosis-specific_Drug|Ethionamide-resistant_mutant|ETHR|RequiresSNPConfirmation	LN:651
@SQ	SN:MEG_8157|Drugs|Aminoglycosides|Aminoglycoside-resistant_idsA|IDSA|RequiresSNPConfirmation	LN:1240
@SQ	SN:MEG_8171|Drugs|Mycobacterium_tuberculosis-specific_Drug|Isoniazid-resistant_mutant|KATG|RequiresSNPConfirmation	LN:2330
@SQ	SN:MEG_8275|Drugs|Mycobacterium_tuberculosis-specific_Drug|MDR_mutant_promoter|NUOA|RequiresSNPConfirmation	LN:752
@SQ	SN:MEG_8444|Drugs|Mycobacterium_tuberculosis-specific_Drug|Pyrazinamide-resistant_mutant|PAND|RequiresSNPConfirmation	LN:420
@SQ	SN:MEG_8636|Drugs|Mycobacterium_tuberculosis-specific_Drug|Pyrazinamide-resistant_mutant|PNCA|RequiresSNPConfirmation	LN:668
@SQ	SN:MEG_8670|Drugs|Oxazolidinone|Oxazolidinone-resistant_23S_rRNA_mutation|RPLC|RequiresSNPConfirmation	LN:654
@SQ	SN:MEG_8675|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOCR|RequiresSNPConfirmation	LN:3951
@SQ	SN:MEG_8679|Drugs|Oxazolidinone|Oxazolidinone-resistant_23S_rRNA_mutation|RRL|RequiresSNPConfirmation	LN:3138
@SQ	SN:MEG_8680|Drugs|Aminoglycosides|Aminoglycoside-resistant_16S_ribosomal_subunit_protein|RRS|RequiresSNPConfirmation	LN:1584
@SQ	SN:MEG_8720|Drugs|Aminoglycosides|Aminoglycoside-resistant_arabinosyltransferase|TLYA|RequiresSNPConfirmation	LN:807
@SQ	SN:MEG_8250|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation	LN:2899
@SQ	SN:MEG_8251|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation	LN:2898
@SQ	SN:MEG_8252|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation	LN:2898
@SQ	SN:MEG_8253|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation	LN:2905
@SQ	SN:MEG_8234|Drugs|MLS|MLS_resistance_ABC_efflux_pumps|MACAB|RequiresSNPConfirmation	LN:228
@SQ	SN:MEG_8263|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_RND_efflux_pumps|MTRC|RequiresSNPConfirmation	LN:1413
@SQ	SN:MEG_8264|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_MFS_efflux_regulator|MTRR|RequiresSNPConfirmation	LN:699
@SQ	SN:MEG_8624|Drugs|betalactams|Class_A_betalactamases|PENA|RequiresSNPConfirmation	LN:1746
@SQ	SN:MEG_8632|Drugs|betalactams|Penicillin_binding_protein|PILQ|RequiresSNPConfirmation	LN:2196
@SQ	SN:MEG_8453|Drugs|betalactams|Penicillin_binding_protein|PBP1|RequiresSNPConfirmation	LN:2397
@SQ	SN:MEG_8637|Drugs|Multi-drug_resistance|MDR_mutant_porin_proteins|POR|RequiresSNPConfirmation	LN:1047
@SQ	SN:MEG_8676|Drugs|Aminoglycosides|Aminoglycoside-resistant_16S_ribosomal_subunit_protein|RPSE|RequiresSNPConfirmation	LN:519
@SQ	SN:MEG_8073|Drugs|Trimethoprim|Dihydrofolate_reductase|DHFR|RequiresSNPConfirmation	LN:1827
@SQ	SN:MEG_8074|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|DHPS|RequiresSNPConfirmation	LN:2121
@SQ	SN:MEG_8242|Drugs|Multi-drug_resistance|MDR_mutation|MDR|RequiresSNPConfirmation	LN:4260
@SQ	SN:MEG_7883|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_RND_efflux_pumps|ACRB|RequiresSNPConfirmation	LN:3150
@SQ	SN:MEG_8146|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2637
@SQ	SN:MEG_8150|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRB|RequiresSNPConfirmation	LN:2415
@SQ	SN:MEG_8446|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2259
@SQ	SN:MEG_8449|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1893
@SQ	SN:MEG_8635|Drugs|Lipopeptides|Colistin-resistant_mutant|PMRAR|RequiresSNPConfirmation	LN:669
@SQ	SN:MEG_8070|Drugs|Trimethoprim|Dihydrofolate_reductase|DFRB|RequiresSNPConfirmation	LN:480
@SQ	SN:MEG_8100|Drugs|Fusidic_acid|Fusidic_acid-resistant_mutation|FUSA|RequiresSNPConfirmation	LN:2082
@SQ	SN:MEG_8142|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GRLA|RequiresSNPConfirmation	LN:2403
@SQ	SN:MEG_8143|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GRLB|RequiresSNPConfirmation	LN:1992
@SQ	SN:MEG_8456|Drugs|betalactams|Penicillin_binding_protein|PBP2|RequiresSNPConfirmation	LN:2184
@SQ	SN:MEG_8459|Drugs|betalactams|Penicillin_binding_protein|PBP4|RequiresSNPConfirmation	LN:1296
@SQ	SN:MEG_8460|Drugs|betalactams|Penicillin_binding_protein|PBP4|RequiresSNPConfirmation	LN:394
@SQ	SN:MEG_8450|Drugs|Aminocoumarins|Aminocoumarin-resistant_DNA_topoisomerases|PARE|RequiresSNPConfirmation	LN:1992
@SQ	SN:MEG_8153|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRB|RequiresSNPConfirmation	LN:2145
@SQ	SN:MEG_8287|Drugs|betalactams|Mutant_porin_proteins|OMP36|RequiresSNPConfirmation	LN:1128
@SQ	SN:MEG_8448|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation	LN:2220
@SQ	SN:MEG_8257|Drugs|Mycobacterium_tuberculosis-specific_Drug|Isoniazid-resistant_mutant|MSHA|RequiresSNPConfirmation	LN:1443
@SQ	SN:MEG_8086|Drugs|Sulfonamides|Sulfonamide-resistant_dihydropteroate_synthases|FOLP|RequiresSNPConfirmation	LN:855
@SQ	SN:MEG_8154|Drugs|Spiropyrimidinetriones|Zoliflodacin-resistant_DNA_topoisomerases|GYRBZ|RequiresSNPConfirmation	LN:2391
@SQ	SN:MEG_8148|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2604
@SQ	SN:MEG_8455|Drugs|betalactams|Penicillin_binding_protein|PBP1|RequiresSNPConfirmation	LN:2397
@SQ	SN:MEG_8633|Drugs|betalactams|Penicillin_binding_protein|PILQ|RequiresSNPConfirmation	LN:2196
@SQ	SN:MEG_8723|Drugs|Pleuromutilin|Pleuromutilin-resistant_23S_rRNA_mutation|UL3|RequiresSNPConfirmation	LN:621
@SQ	SN:MEG_8149|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation	LN:2484
@SQ	SN:MEG_8152|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRB|RequiresSNPConfirmation	LN:2322
@SQ	SN:MEG_8458|Drugs|betalactams|Penicillin_binding_protein|PBP3|RequiresSNPConfirmation	LN:1848
@SQ	SN:MEG_8457|Drugs|betalactams|Penicillin_binding_protein|PBP2|RequiresSNPConfirmation	LN:1767
@SQ	SN:MEG_8454|Drugs|betalactams|Penicillin_binding_protein|PBP1|RequiresSNPConfirmation	LN:1980
@SQ	SN:MEG_8101|Drugs|Metronidazole|fxr_nitroimidazole_reductase|FXRA|RequiresSNPConfirmation	LN:654
@SQ	SN:MEG_8673|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation	LN:8673
@SQ	SN:MEG_8671|Drugs|MLS|Macrolide-resistant_50S_rRNA_mutation|RPLD|RequiresSNPConfirmation	LN:621
@SQ	SN:MEG_8254|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation	LN:2910
@SQ	SN:MEG_8443|Drugs|Pleuromutilin|Pleuromutilin-resistant_23S_rRNA_mutation|P23S|RequiresSNPConfirmation	LN:2893
"""

gene_dict = {}

full_gene_name = []
full_gene_name.append("MEG_4095|Drugs|Fosfomycin|Fosfomycin_target_mutation|MURA|RequiresSNPConfirmation")
full_gene_name.append("MEG_4130|Drugs|Mycobacterium_tuberculosis-specific_Drug|Isoniazid-resistant_mutant|NDH|RequiresSNPConfirmation")
full_gene_name.append("MEG_5328|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation")
full_gene_name.append("MEG_5331|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|PARC|RequiresSNPConfirmation")
full_gene_name.append("MEG_5401|Drugs|betalactams|Penicillin_binding_protein|PBP2|RequiresSNPConfirmation")
full_gene_name.append("MEG_411|Multi-compound|Drug_and_biocide_resistance|Drug_and_biocide_RND_efflux_regulator|ACRR|RequiresSNPConfirmation")
full_gene_name.append("MEG_6090|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation")
full_gene_name.append("MEG_6094|Drugs|Rifampin|Rifampin-resistant_beta-subunit_of_RNA_polymerase_RpoB|RPOB|RequiresSNPConfirmation")
full_gene_name.append("MEG_6142|Drugs|Mycobacterium_tuberculosis-specific_Drug|Pyrazinamide-resistant_mutant|RPSA|RequiresSNPConfirmation")

# Example: transforms 3M2I3M to MMMIIMMM
def extend_cigar(cigar):
    extended_cigar = ""
    count = 0
    for c in cigar:
        if c.isdigit():
            count = count * 10 + int(c)
        else:
            while count > 0:
                extended_cigar = extended_cigar + c
                count -= 1
    return extended_cigar

def snp_test(snp_num, full_name, start, end, cigar, SEQ):
    line1 = []
    line2 = []
    line1.append("test_" + str(snp_num))
    line2.append(line1[0])
    line1.append("99")
    line2.append("147")
    line1.append(full_name)
    line2.append(line1[2])
    line1.append(str(start))
    line2.append(line1[3])
    line1.append("255")
    line2.append(line1[4])
    line1.append(cigar)
    line2.append(cigar)
    line1.append("=")
    line2.append(line1[6])
    line1.append(line1[3])
    line2.append(line1[3])
    line1.append(str(end-(start-1)))
    line2.append("-" + line1[8])
    line1.append(SEQ)
    line2.append(SEQ)
    line1.append("*")
    line2.append(line1[10])
    to_return = ""
    for tab in line1:   to_return = to_return + tab + "\t"
    to_return = to_return + "\n"
    for tab in line2:   to_return = to_return + tab + "\t"
    to_return = to_return + "\n"
    return to_return

def make_test(snp_num, full_name, seq, cigar, start, end, seq_is_nt=False):
    nt_seq = seq
    if not seq_is_nt:
        nt_seq = ""
        for aa in seq:   nt_seq = nt_seq + reverseTranslation(aa)
        index_start = (start-1) % 3
        extended_cigar = extend_cigar(cigar)
        index_end = index_start + extended_cigar.count('M') + extended_cigar.count('I')
        nt_seq = nt_seq[index_start:index_end]
    return snp_test(snp_num, full_name, start, end, cigar, nt_seq)

def sort_snp_info(snp_info):
    snp_first = snp_info[0:int(len(snp_info)/2)]
    if len(snp_first) > 1:      snp_first = sort_snp_info(snp_first)

    snp_last = snp_info[int(len(snp_info)/2):len(snp_info)]
    if len(snp_last) > 1:       snp_last = sort_snp_info(snp_last)

    to_return = []
    for i in range(0, len(snp_info)):
        if len(snp_first) == 0:         to_return.append(snp_last.pop(0))
        elif len(snp_last) == 0:        to_return.append(snp_first.pop(0))

        else:
            if ((snp_first[0][1][0] if type(snp_first[0][1]) == list else snp_first[0][1]) < 
                    (snp_last[0][1][0] if type(snp_last[0][1]) == list else snp_last[0][1])):
                to_return.append(snp_first.pop(0))
            else:
                to_return.append(snp_last.pop(0))
    return to_return

def nucleic_must_test(must_info, gene, snp_num):
    start = must_info.getPos()-26
    if start < 1: start = 1

    gene.seq = gene.ntSequence()
    SEQ = gene.seq[start-1:(must_info.getPos()-1)]
    
    cigar = "20H"
    last_nt = must_info.getPos()-1

    while must_info != None:
        if last_nt < must_info.getPos()-1:  
            SEQ = SEQ + gene.seq[last_nt:must_info.getPos()-1]

        SEQ = SEQ + must_info.getWt()
        last_nt = must_info.getPos()
        end = must_info.getPos()+ 25
        must_info = must_info.getNext()

    if end > gene.ntSeqLength(): end = gene.ntSeqLength()
    cigar = cigar + str(end - (start - 1)) + "M20H"
    SEQ = SEQ + gene.seq[last_nt:end]
    return snp_test(snp_num, gene.getFullName(), start, end, cigar, str(SEQ))

def must_test(must_info, gene, snp_num):
    if gene.rRna(): return nucleic_must_test(must_info, gene, snp_num)
    start = must_info.getPos()*3-23
    if start < 1: start = 1
    gene.seq = gene.ntSequence()
    SEQ = gene.seq[start-1:(must_info.getPos()*3-3)]
    cigar = "20H"
    last_nt = must_info.getPos()*3-3
    while must_info != None:
        if last_nt < must_info.getPos()*3-3:
            SEQ = SEQ + gene.seq[last_nt:must_info.getPos()*3-3]
        SEQ = SEQ + reverseTranslation(must_info.getWt())
        last_nt = must_info.getPos()*3
        end = must_info.getPos()*3 + 21
        must_info = must_info.getNext()
    if end > gene.ntSeqLength(): end = gene.ntSeqLength()
    cigar = cigar + str(end - (start - 1)) + "M20H"
    SEQ = SEQ + gene.seq[last_nt:end]
    return snp_test(snp_num, gene.getFullName(), start, end, cigar, str(SEQ))

def multiple_snp_test(snp_info, gene, snp_num):
    if gene.rRna(): return multiple_nucleic_snp_test(snp_info, gene, snp_num)
    snp_info = sort_snp_info(snp_info)
    start = (snp_info[0][1][0] if type(snp_info[0][1]) == list else snp_info[0][1])*3-23
    end = (snp_info[-1][1][0] if type(snp_info[-1][1]) == list else snp_info[-1][1])*3+21
    if (start < 1) : start = 1
    if end > gene.ntSeqLength(): end = gene.ntSeqLength()
    in_between_H = []
    last_nt = start - 1
    last_is_M = False
    last_is_I = False
    last_is_D = False
    for snp in snp_info:
        if (last_nt < (snp[1][0] if type(snp[1]) == list else snp[1])*3 - 3):
            if not(last_is_M): #wt435-, wt437-/mt
                in_between_H.append(((((snp[1][0] if type(snp[1]) == list else snp[1])*3-3) - last_nt),'M'))
                last_is_M = True
                last_is_I = False
                last_is_D = False
            else: #wt435mt, wt437-/mt
                in_between_H[-1] = (in_between_H[-1][0] + (((snp[1][0] if type(snp[1]) == list else snp[1])*3-3) - last_nt), 'M')
        if snp[2][0] == "-":
            if not(last_is_D): #wt436mt, wt437-; wt435-/mt, wt437-
                last_is_D = True
                last_is_I = False
                last_is_M = False
                in_between_H.append((3, 'D'))
            else: #wt436-,wt437-
                in_between_H[-1] = (in_between_H[-1][0] + 3, 'D')
        elif snp[2][0] == "+":
            if not(last_is_I): #wt436mt, wt437-; wt435-/mt, wt437-
                last_is_M = False
                last_is_D = False
                last_is_I = True
                in_between_H.append((3 * len(snp[0]), 'I'))
            else:
                in_between_H[-1] = (in_between_H[-1][0] + 3 * len(snp[0]), 'I')

        else:
            if (last_is_M): #wt436mt, wt437mt; #wt435-/mt, wt437mt:
                in_between_H[-1] = (in_between_H[-1][0] + 3, 'M')
            else: #wt436-,wt437mt
                last_is_M = True
                last_is_D = False
                last_is_I = False
                in_between_H.append((3, 'M'))
        if (snp[2]) == "-":
            last_nt = snp[1][0]*3
        elif (snp[2]) != "+":
            last_nt = snp[1]*3
    if last_is_M:
        in_between_H[-1] = (in_between_H[-1][0] + end - (snp_info[-1][1][0] if type(snp_info[-1][1]) == list else snp_info[-1][1])*3, 'M')
    elif end > (snp_info[-1][1][0] if type(snp_info[-1][1]) == list else snp_info[-1][1])*3:
        in_between_H.append((end - last_nt, 'M'))
    cigar = "20H"
    for next in in_between_H:
        cigar = cigar + str(next[0]) + next[1]
    cigar = cigar + "20H"
    gene.seq = gene.ntSequence()
    SEQ = gene.seq[start-1:((snp_info[0][1][0] if type(snp_info[0][1]) == list else snp_info[0][1])*3-3)]
    for i in range(0, len(snp_info)):
        if snp_info[i][2][-1] != '-':
            if snp_info[i][2][-1] == '+':
                for aa in snp_info[i][0]:
                    SEQ = SEQ + reverseTranslation(aa)
            else:
                SEQ = SEQ + reverseTranslation(snp_info[i][2][-1])
        if i < (len(snp_info) - 1):
            SEQ = SEQ + gene.seq[(snp_info[i][1][0] if type(snp_info[i][1]) == list else snp_info[i][1])*3:(snp_info[i+1][1][0] if type(snp_info[i+1][1]) == list else snp_info[i+1][1])*3-3]
        else:
            SEQ = SEQ + gene.seq[last_nt:end]
    return snp_test(snp_num, gene.getFullName(), start, end, cigar, str(SEQ))

def multiple_nucleic_snp_test(snp_info, gene, snp_num):
    snp_info = sort_snp_info(snp_info)
    start = snp_info[0][1]-26
    end = snp_info[-1][1]+25
    if (start < 1) : start = 1
    if end > gene.ntSeqLength(): end = gene.ntSeqLength()
    in_between_H = []
    last_nt = start - 1
    last_is_M = False
    for snp in snp_info:
        if (last_nt < (snp[1][0] if type(snp[1]) == list else snp[1])-1):
            if not(last_is_M): #wt435-, wt437-/mt
                in_between_H.append((((snp[1][0] if type(snp[1]) == list else snp[1])-1 - last_nt),'M'))
                last_is_M = True
            else: #wt435mt, wt437-/mt
                in_between_H[-1] = (in_between_H[-1][0] + ((snp[1][0] if type(snp[1]) == list else snp[1])-1 - last_nt), 'M')
        if snp[2][0] == "-":
            if (last_is_M): #wt436mt, wt437-; wt435-/mt, wt437-
                last_is_M = False
                in_between_H.append((1, 'D'))
            else: #wt436-,wt437-
                in_between_H[-1] = (in_between_H[-1][0] + 1, 'D')
        else:
            if (last_is_M): #wt436mt, wt437mt; #wt435-/mt, wt437mt:
                in_between_H[-1] = (in_between_H[-1][0] + 1, 'M')
            else: #wt436-,wt437mt
                last_is_M = True
                in_between_H.append((1, 'M'))
        if (type(snp[1]) == list):
            last_nt = snp[1][0]
        else:
            last_nt = snp[1]
    if last_is_M:
        in_between_H[-1] = (in_between_H[-1][0] + end - (snp_info[-1][1][0] if type(snp_info[-1][1]) == list else snp_info[-1][1]), 'M')
    elif end > (snp_info[-1][1][0] if type(snp_info[-1][1]) == list else snp_info[-1][1]):
        in_between_H.append((end - (snp_info[-1][1][0] if type(snp_info[-1][1]) == list else snp_info[-1][1]), 'M'))
    cigar = "20H"
    for next in in_between_H:
        cigar = cigar + str(next[0]) + next[1]
    cigar = cigar + "20H"
    gene.seq = gene.ntSequence()
    SEQ = gene.seq[start-1:((snp_info[0][1][0] if type(snp_info[0][1]) == list else snp_info[0][1])-1)]
    for i in range(0, len(snp_info)):
        if snp_info[i][2][-1] != '-':
            SEQ = SEQ + snp_info[i][2][-1]
        if i < (len(snp_info) - 1):
            SEQ = SEQ + gene.seq[(snp_info[i][1][0] if type(snp_info[i][1]) == list else snp_info[i][1]):(snp_info[i+1][1][0] if type(snp_info[i+1][1]) == list else snp_info[i+1][1])-1]
        else:
            SEQ = SEQ + gene.seq[(snp_info[i][1][0] if type(snp_info[i][1]) == list else snp_info[i][1]):end]
    return snp_test(snp_num, gene.getFullName(), start, end, cigar, str(SEQ))

def deletion_test(snp_info, gene, snp_num):
    to_return = ""
    for i in snp_info[1]:
        start = i*3-23
        end = i*3+21
        if (start < 1) : start = 1
        if end > gene.ntSeqLength(): end = gene.ntSeqLength()
        cigar = "20H"
        if i != 1:
            cigar = cigar + str((i*3-3)-(start-1))+ "M3D"
        if i != gene.ntSeqLength():
            cigar = cigar + str(end-i*3) + "M20H"
        gene.seq = gene.ntSequence()
        SEQ = gene.seq[start-1:(i*3-3)] + gene.seq[(i*3):end]
        to_return = to_return + snp_test(str(snp_num) +"a"+str(i), gene.getFullName(), start, end, cigar, str(SEQ))
        cigar = "20H"
        if i != 1:
            cigar = cigar + str((i*3-4)-(start-1)) + "M3D"
        if i != gene.ntSeqLength():
            cigar = cigar + str(end-(i*3)+1) + "M20H"
        gene.seq = gene.ntSequence()
        SEQ = gene.seq[start-1:(i*3-4)] + gene.seq[(i*3)-1:end]
        to_return = to_return + snp_test(str(snp_num) +"b"+str(i), gene.getFullName(), start, end, cigar, str(SEQ))
        cigar = "20H"
        if i != 1:
            cigar = cigar + str((i*3-2)-(start-1)) + "M3D"
        if i != gene.ntSeqLength():
            cigar = cigar + str(end-(i*3)-1) + "M20H"
        gene.seq = gene.ntSequence()
        SEQ = gene.seq[start-1:(i*3-2)] + gene.seq[(i*3)+1:end]
        to_return = to_return + snp_test(str(snp_num) +"c"+str(i), gene.getFullName(), start, end, cigar, str(SEQ))
    return to_return

def insertion_test(snp_info, gene, snp_num):
    to_return = ""
    for i in snp_info[1]:
        start = i*3-23
        end = i*3+21
        if (start < 1) : start = 1
        if end > gene.ntSeqLength(): end = gene.ntSeqLength()
        cigar = "20H"
        if i != 1:
            cigar = cigar + str((i*3-3)-(start-1))+ "M" + str(3 * len(snp_info[0])) + "I"
        if i != gene.ntSeqLength():
            cigar = cigar + str(end-i*3+3) + "M20H"
        gene.seq = gene.ntSequence()
        inserted = ""
        for aa in snp_info[0]:
            inserted = inserted + reverseTranslation(aa)
        SEQ = gene.seq[start-1:(i*3-3)] + inserted + gene.seq[(i*3)-3:end]
        to_return = to_return + snp_test(str(snp_num) +"a"+str(i), gene.getFullName(), start, end, cigar, str(SEQ))
    if (len(snp_info[1]) > 1) and (len(snp_info[0]) > 1):
        start = snp_info[1][0]*3-20
        end = snp_info[1][0]*3+24
        if (start < 1) : start = 1
        if end > gene.ntSeqLength(): end = gene.ntSeqLength()
        cigar = "20H"
        if snp_info[1][0] != 1:
            cigar = cigar + str((snp_info[1][0]*3)-(start-1))+ "M" + str(3 * len(snp_info[0])) + "I"
        if snp_info[1][0] != gene.ntSeqLength():
            cigar = cigar + str(end-snp_info[1][0]*3) + "M20H"
        gene.seq = gene.ntSequence()
        inserted = ""
        for aa in snp_info[0][1:]:
            inserted = inserted + reverseTranslation(aa)
        inserted = inserted + reverseTranslation(snp_info[0][0])
        SEQ = gene.seq[start-1:(snp_info[1][0]*3)] + inserted + gene.seq[(snp_info[1][0]*3):end]
        to_return = to_return + snp_test(str(snp_num) +"a"+str(snp_info[1]), gene.getFullName(), start, end, cigar, str(SEQ))
    return to_return

def single_snp_test(snp_info, gene, snp_num):
    if gene.rRna(): return single_nucleic_snp_test(snp_info, gene, snp_num)
    if snp_info[2] == "-": return deletion_test(snp_info, gene, snp_num)
    elif snp_info[2] == "+": return insertion_test(snp_info, gene, snp_num)
    start = snp_info[1]*3-23
    end = snp_info[1]*3+21
    if (start < 1) : start = 1
    if end > gene.ntSeqLength(): end = gene.ntSeqLength()
    cigar = "20H" + str(end - (start -1)) + "M20H"
    if snp_info[2][0] == "-":
        cigar = "20H"
        if snp_info[1] != 1:
            cigar = cigar + str((snp_info[1]*3-3)-(start-1)) + "M3D"
        if snp_info[1] != gene.ntSeqLength():
            cigar = cigar + str(end-snp_info[1]*3)+ "M20H"
    gene.seq = gene.ntSequence()
    SEQ = gene.seq[start-1:(snp_info[1]*3-3)] + reverseTranslation(snp_info[2][-1]) + gene.seq[(snp_info[1]*3):end]
    return snp_test(snp_num, gene.getFullName(), start, end, cigar, str(SEQ))

def single_nucleic_snp_test(snp_info, gene, snp_num):
    snpPos = (snp_info[1][0] if type(snp_info[1]) == list else snp_info[1])
    start = snpPos - 26
    end = snpPos + 25
    if start < 1: start = 1
    if end > gene.ntSeqLength(): end = gene.ntSeqLength()
    cigar = "20H" + str(end - (start -1)) + "M20H"
    if snp_info[2][0] == "-":
        cigar = "20H"
        if snpPos != 1:
            cigar = cigar + str((snpPos-1)-(start-1)) + "M1D"
        if snpPos != gene.ntSeqLength():
            cigar = cigar + str(end-snpPos) + "M20H"
    gene.seq = gene.ntSequence()
    if snp_info[2][0] == "-":
        SEQ = gene.seq[start-1:(snpPos-1)] + gene.seq[(snpPos):end]
    else:
        SEQ = gene.seq[start-1:(snpPos-1)] + snp_info[2][-1] + gene.seq[(snpPos):end]
    return snp_test(snp_num, gene.getFullName(), start, end, cigar, str(SEQ))


# Standard test of all variants found so far
def test1():
    SAM_file = open('test/test1.sam', 'w')
    snp_num = 1
    SAM_file.write(header)
    for gene in gene_dict.values():
        for snp in gene.condensedMisInDelInfo():
            SAM_file.write(single_snp_test(snp, gene, snp_num))
            snp_num+=1
        for snp in gene.condensedMultInfo():
            SAM_file.write(multiple_snp_test(snp[0], gene, snp_num))
            snp_num+=1
        if not gene.rRna():
            for snp in gene.condensedNonInfo():
                SAM_file.write(single_snp_test(snp, gene, snp_num))
                snp_num+=1
        if gene.getGeneTag() == 'I':
            must = gene.getFirstMustBetweenParams(1, gene.ntSeqLength())
            if must == None:
                continue
            SAM_file.write(must_test(must[0], gene, snp_num))
            snp_num+=1
            if must[0].getNext() != None:
                SAM_file.write(must_test(must[0].getNext(), gene, snp_num))
                snp_num+=1
    SAM_file.close()
    pysam.sort('-o', 'test/test1.bam', 'test/test1.sam')

def test2():
    SAM_file = open('test/test2.sam', 'w')
    SAM_file.write(header)
    snp_num = 1
    SAM_file.write(multiple_snp_test([('Y', 147, ('*',)), ('G', 146, ('*',))], gene_dict.get("MEG_2866|Drugs|Mycobacterium_tuberculosis-specific_Drug|Ethionamide-resistant_mutant|ETHA|RequiresSNPConfirmation"), snp_num))
    snp_num+=1
    SAM_file.write(multiple_snp_test([('Y', 147, ('*',)), ('S', 149, ('*',))], gene_dict.get("MEG_2866|Drugs|Mycobacterium_tuberculosis-specific_Drug|Ethionamide-resistant_mutant|ETHA|RequiresSNPConfirmation"), snp_num))
    snp_num+=1
    SAM_file.write(multiple_snp_test([('Y', 147, ('*',)), ('C', 131, ('*',))], gene_dict.get("MEG_2866|Drugs|Mycobacterium_tuberculosis-specific_Drug|Ethionamide-resistant_mutant|ETHA|RequiresSNPConfirmation"), snp_num))
    snp_num+=1
    SAM_file.write(multiple_snp_test([('W', 98, ('*',)), ('R', 99, ('*',))], gene_dict.get("MEG_7250|Drugs|Mycobacterium_tuberculosis-specific_Drug|Para-aminosalicylic_acid_resistant_mutant|THYA|RequiresSNPConfirmation"), snp_num))
    snp_num+=1
    SAM_file.write(multiple_snp_test([('W', 98, ('*',)), ('S', 105, ('*',))], gene_dict.get("MEG_7250|Drugs|Mycobacterium_tuberculosis-specific_Drug|Para-aminosalicylic_acid_resistant_mutant|THYA|RequiresSNPConfirmation"), snp_num))
    snp_num+=1
    SAM_file.write(multiple_snp_test([('W', 98, ('*',)), ('V', 96, ('*',))], gene_dict.get("MEG_7250|Drugs|Mycobacterium_tuberculosis-specific_Drug|Para-aminosalicylic_acid_resistant_mutant|THYA|RequiresSNPConfirmation"), snp_num))
    snp_num+=1
    SAM_file.write(multiple_snp_test([('G', 84, ('D',)), ('L', 100, ('P',)), ('A', 115, ('T',)), ('Y', 122, ('N',))], gene_dict.get("MEG_4132|Drugs|Mycobacterium_tuberculosis-specific_Drug|Isoniazid-resistant_mutant|NDH|RequiresSNPConfirmation"), snp_num))
    snp_num+=1
    SAM_file.write(multiple_snp_test([('Q', 10, ('*',)), ('D', 12, ('A',)), ('L', 19, ('P',)), ('C', 14, ('Y',)), ('G', 23, ('V',))], gene_dict.get("MEG_5803|Drugs|Mycobacterium_tuberculosis-specific_Drug|Pyrazinamide-resistant_mutant|PNCA|RequiresSNPConfirmation"), snp_num))
    snp_num+=1
    SAM_file.write(multiple_snp_test([('G', 88, ('C',)), ('S', 91, ('P',)), ('T', 80, ('A',)), ('A', 90, ('G',))], gene_dict.get("MEG_3180|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation"), snp_num))
    snp_num+=1
    SAM_file.write(multiple_snp_test([('G', 88, ('C',)), ('S', 91, ('P',)), ('T', 80, ('A',)), ('A', 90, ('V',))], gene_dict.get("MEG_3180|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation"), snp_num))
    snp_num+=1
    SAM_file.write(multiple_snp_test([('G', 88, ('C',)), ('S', 91, ('P',)), ('T', 80, ('A',))], gene_dict.get("MEG_3180|Drugs|Fluoroquinolones|Fluoroquinolone-resistant_DNA_topoisomerases|GYRA|RequiresSNPConfirmation"), snp_num))
    SAM_file.close()
    pysam.sort('-o', 'test/test2.bam', 'test/test2.sam')

def insertion_test1():
    SAM_file = open("test/insertion1.sam", "w")
    SAM_file.write(header)
    snp_num = 1
    SAM_file.write(make_test(snp_num, full_gene_name[0], "TNLR", "20H1M3I8M20H", 369 * 3 - 2, 371 * 3)) #test with mt second
    snp_num  += 1
    SAM_file.write(make_test(snp_num, full_gene_name[1], "CFRH", "20H1M3I8M20H", 13 * 3 - 2, 15 * 3)) #test with mt first
    snp_num  += 1
    SAM_file.write(make_test(snp_num, full_gene_name[2], "WWSI", "20H1M3I8M20H", 83 * 3 - 2, 85 * 3)) #test with mt in both
    snp_num  += 1
    SAM_file.write(make_test(snp_num, full_gene_name[3], "VINN", "20H1M3I8M20H", 103 * 3 - 2, 105 * 3)) #test with wt first, mt second
    snp_num  += 1
    SAM_file.write(make_test(snp_num, full_gene_name[4], "VARK", "20H1M3I8M20H", 501 * 3 - 2, 503 * 3)) #test with mt first, wt second
    SAM_file.close()
    pysam.sort('-o', 'test/insertion1.bam', 'test/insertion1.sam')

def insertion_test2():
    SAM_file = open("test/insertion2_1.sam", "w")
    SAM_file.write(header)
    snp_num = 1
    SAM_file.write(make_test(snp_num, full_gene_name[5], "TTCGA", "20H3M3I9M20H", 44 * 3 - 2, 47 * 3)) #test with mt after insertion
    SAM_file.close()
    pysam.sort('-o', 'test/insertion2_1.bam', 'test/insertion2_1.sam')

    SAM_file = open("test/insertion2_2.sam", "w")
    SAM_file.write(header)
    SAM_file.write(make_test(snp_num, full_gene_name[5], "TRCGA", "20H3M3I9M20H", 44 * 3 - 2, 47 * 3)) #test with wt in insertion, mt after insertion
    SAM_file.close()
    pysam.sort('-o', 'test/insertion2_2.bam', 'test/insertion2_2.sam')

    SAM_file = open("test/insertion2_3.sam", "w")
    SAM_file.write(header)
    SAM_file.write(make_test(snp_num, full_gene_name[5], "TCRGA", "20H3M3I9M20H", 44 * 3 - 2, 47 * 3)) #test with mt in insertion, wt after insertion
    SAM_file.close()
    pysam.sort('-o', 'test/insertion2_3.bam', 'test/insertion2_3.sam')

def insertion_test3():
    SAM_file = open("test/insertion3.sam", "w")
    SAM_file.write(header)
    snp_num = 1
    SAM_file.write(make_test(snp_num, full_gene_name[0], "ATNLR", "20H3M1I3M2I6M20H", 368 * 3 - 2, 371 * 3)) #test with mt second
    snp_num  += 1
    SAM_file.write(make_test(snp_num, full_gene_name[1], "PCFRH", "20H3M1I3M2I6M20H", 12 * 3 - 2, 15 * 3)) #test with mt first
    snp_num  += 1
    SAM_file.write(make_test(snp_num, full_gene_name[2], "DWWSI", "20H3M1I3M2I6M20H", 82 * 3 - 2, 85 * 3)) #test with mt in both
    snp_num  += 1
    SAM_file.write(make_test(snp_num, full_gene_name[3], "TVINN", "20H3M1I3M2I6M20H", 102 * 3 - 2, 105 * 3)) #test with wt first, mt second
    snp_num  += 1
    SAM_file.write(make_test(snp_num, full_gene_name[4], "TVARK", "20H3M1I3M2I6M20H", 500 * 3 - 2, 503 * 3)) #test with mt first, wt second
    SAM_file.close()
    pysam.sort('-o', 'test/insertion3.bam', 'test/insertion3.sam')

def insertion_test4():
    SAM_file = open("test/insertion4_1.sam", "w")
    SAM_file.write(header)
    snp_num = 1
    SAM_file.write(make_test(snp_num, full_gene_name[0], "MATNL", "20H3M1I6M2I3M20H", 367 * 3 - 2, 370 * 3)) #test with mt second; second insertion chunk
    snp_num  += 1
    SAM_file.write(make_test(snp_num, full_gene_name[1], "PPCFR", "20H3M1I6M2I3M20H", 11 * 3 - 2, 14 * 3)) #test with mt first; second insertion chunk
    snp_num  += 1
    SAM_file.write(make_test(snp_num, full_gene_name[2], "GDWWS", "20H3M1I6M2I3M20H", 81 * 3 - 2, 84 * 3)) #test with mt in both; second insertion chunk
    snp_num  += 1
    SAM_file.write(make_test(snp_num, full_gene_name[3], "RTVIN", "20H3M1I6M2I3M20H", 101 * 3 - 2, 104 * 3)) #test with wt first, mt second; second insertion chunk
    snp_num  += 1
    SAM_file.write(make_test(snp_num, full_gene_name[4], "GTVAR", "20H3M1I6M2I3M20H", 499 * 3 - 2, 502 * 3)) #test with mt first, wt second; second insertion chunk
    SAM_file.close()
    pysam.sort('-o', 'test/insertion4_1.bam', 'test/insertion4_1.sam')

    SAM_file = open("test/insertion4_2.sam", "w")
    SAM_file.write(header)
    snp_num = 1
    SAM_file.write(make_test(snp_num, full_gene_name[0], "ATNLR", "20H3M1I6M2I3M20H", 368 * 3 - 2, 371 * 3)) #test with mt second; first insertion chunk
    snp_num  += 1
    SAM_file.write(make_test(snp_num, full_gene_name[1], "PCFRH", "20H3M1I6M2I3M20H", 12 * 3 - 2, 15 * 3)) #test with mt first; first insertion chunk
    snp_num  += 1
    SAM_file.write(make_test(snp_num, full_gene_name[2], "DWWSI", "20H3M1I6M2I3M20H", 82 * 3 - 2, 85 * 3)) #test with mt in both; first insertion chunk
    snp_num  += 1
    SAM_file.write(make_test(snp_num, full_gene_name[3], "TVINN", "20H3M1I6M2I3M20H", 102 * 3 - 2, 105 * 3)) #test with wt first, mt second; first insertion chunk
    snp_num  += 1
    SAM_file.write(make_test(snp_num, full_gene_name[4], "TVARK", "20H3M1I6M2I3M20H", 500 * 3 - 2, 503 * 3)) #test with mt first, wt second; first insertion chunk
    SAM_file.close()
    pysam.sort('-o', 'test/insertion4_2.bam', 'test/insertion4_2.sam')

def insertion_test5():
    SAM_file = open("test/insertion5.sam", "w")
    SAM_file.write(header)
    snp_num = 1
    SAM_file.write(make_test(snp_num, full_gene_name[5], "RCGA", "20H1M3I7M20H", 45 * 3 - 1, 47 * 3)) #test with mt in insertion
    SAM_file.close()
    pysam.sort('-o', 'test/insertion5.bam', 'test/insertion5.sam')

def deletion_test1():
    SAM_file = open("test/deletion1.sam", "w")
    SAM_file.write(header)
    snp_num = 1
    SAM_file.write(make_test(snp_num, full_gene_name[6], "HQN", "20H1M9D8M20H", 432 * 3 - 2, 437 * 3)) #test with deletion mutations
    snp_num  += 1
    SAM_file.write(make_test(snp_num, full_gene_name[5], "CGA", "20H1M3D8M20H", 44 * 3 - 2, 47 * 3)) #test with deletion making mt
    SAM_file.close()
    pysam.sort('-o', 'test/deletion1.bam', 'test/deletion1.sam')

def deletion_test2():
    SAM_file = open("test/deletion2.sam", "w")
    SAM_file.write(header)
    snp_num = 1
    SAM_file.write(make_test(snp_num, full_gene_name[5], "GVTC", "20H1M2D8M1D2M20H", 41 * 3 - 1, 45 * 3)) #test with mt at last deletion chunk
    SAM_file.close()
    pysam.sort('-o', 'test/deletion2.bam', 'test/deletion2.sam')

def insertion_test6():
    SAM_file = open("test/insertion6.sam", "w")
    SAM_file.write(header)
    snp_num = 1
    SAM_file.write(make_test(snp_num, full_gene_name[5], "GVTCG", "20H1M1I6M2I4M20H", 43 * 3 - 1, 46 * 3)) #test with mt last insertion chunk
    SAM_file.close()
    pysam.sort('-o', 'test/insertion6.bam', 'test/insertion6.sam')

def deletion_test3():
    SAM_file = open("test/deletion3_1.sam", "w")
    SAM_file.write(header)
    snp_num = 1
    SAM_file.write(make_test(snp_num, full_gene_name[5], "CGAI", "20H1M1D8M2D3M20H", 45 * 3 - 2, 49 * 3)) #test with mt at first deletion chunk
    SAM_file.close()
    pysam.sort('-o', 'test/deletion3_1.bam', 'test/deletion3_1.sam')

    SAM_file = open("test/deletion3_2.sam", "w")
    SAM_file.write(header)
    snp_num = 1
    SAM_file.write(make_test(snp_num, full_gene_name[5], "GVTC", "20H1M1D8M2D3M20H", 41 * 3 - 2, 45 * 3)) #test with mt at last deletion chunk
    SAM_file.close()
    pysam.sort('-o', 'test/deletion3_2.bam', 'test/deletion3_2.sam')

def insertion_deletion_test1():
    SAM_file = open("test/insertion_deletion1_1.sam", "w")
    SAM_file.write(header)
    snp_num = 1
    SAM_file.write(make_test(snp_num, full_gene_name[5], "TCGAI", "20H3M1D9M1I2M20H", 44 * 3 - 2, 48 * 3)) #test with mt at deletion chunk
    SAM_file.close()
    pysam.sort('-o', 'test/insertion_deletion1_1.bam', 'test/insertion_deletion1_1.sam')

    SAM_file = open("test/insertion_deletion1_2.sam", "w")
    SAM_file.write(header)
    snp_num = 1
    SAM_file.write(make_test(snp_num, full_gene_name[5], "AGVTC", "20H3M1D9M1I2M20H", 41 * 3 - 2, 45 * 3)) #test with mt at insertion chunk
    SAM_file.close()
    pysam.sort('-o', 'test/insertion_deletion1_2.bam', 'test/insertion_deletion1_2.sam')

def insertion_deletion_test2():
    SAM_file = open("test/insertion_deletion2_1.sam", "w")
    SAM_file.write(header)
    snp_num = 1
    SAM_file.write(make_test(snp_num, full_gene_name[5], "CGAI", "20H1M1I8M1D2M20H", 45 * 3 - 2, 48 * 3)) #test with mt at insertion chunk
    SAM_file.close()
    pysam.sort('-o', 'test/insertion_deletion2_1.bam', 'test/insertion_deletion2_1.sam')

    SAM_file = open("test/insertion_deletion2_2.sam", "w")
    SAM_file.write(header)
    snp_num = 1
    SAM_file.write(make_test(snp_num, full_gene_name[5], "GVTC", "20H1M1I8M1D2M20H", 42 * 3 - 2, 45 * 3)) #test with mt at deletion chunk
    SAM_file.close()
    pysam.sort('-o', 'test/insertion_deletion2_2.bam', 'test/insertion_deletion2_2.sam')

def insertion_deletion_test3():
    SAM_file = open("test/insertion_deletion3.sam", "w")
    SAM_file.write(header)
    snp_num = 1
    SAM_file.write(make_test(snp_num, full_gene_name[5], "TCG", "20H4M1D1M1I3M20H", 44 * 3 - 2, 46 * 3)) 
    SAM_file.close()
    pysam.sort('-o', 'test/insertion_deletion3.bam', 'test/insertion_deletion3.sam')

def insertion_deletion_test4():
    SAM_file = open("test/insertion_deletion4.sam", "w")
    SAM_file.write(header)
    snp_num = 1
    SAM_file.write(make_test(snp_num, full_gene_name[5], "CG", "20H1M1I1M1D3M20H", 45 * 3 - 2, 46 * 3)) 
    SAM_file.close()
    pysam.sort('-o', 'test/insertion_deletion4.bam', 'test/insertion_deletion4.sam')

def deletion_test4():
    SAM_file = open("test/deletion4.sam", "w")
    SAM_file.write(header)
    snp_num = 1
    SAM_file.write(make_test(snp_num, full_gene_name[6], "FNP", "20H1M12D6M20H", 433 * 3, 439 * 3)) #test with deletion mutations
    SAM_file.close()
    pysam.sort('-o', 'test/deletion4.bam', 'test/deletion4.sam')

def deletion_test5():
    SAM_file = open("test/deletion5.sam", "w")
    SAM_file.write(header)
    snp_num = 1
    SAM_file.write(make_test(snp_num, full_gene_name[6], "FN", "20H3M12D1M20H", 433 * 3 - 2, 438 * 3 - 2)) #test with deletion mutations
    SAM_file.close()
    pysam.sort('-o', 'test/deletion5.bam', 'test/deletion5.sam')

def frameshift_test():
    SAM_file = open("test/frameshift.sam", "w")
    SAM_file.write(header)
    snp_num = 1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "ITHKRRISRTRPRRSDP*T", "20H24M1I32M20H", 524 * 3 - 2, 542 * 3 - 1))   #Res
    snp_num += 1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "ITHKRRISRTRPGGLTRER", "20H24M1I9M1D23M20H", 524 * 3 - 2, 542 * 3))   #Res*
    snp_num += 1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "ITHKRRISRTGPGGLTRER", "20H24M1I6M1D26M20H", 524 * 3 - 2, 542 * 3))   #Sus
    snp_num += 1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "ITHKRRISRTRPRGLTRER", "20H24M1I15M1D17M20H", 524 * 3 - 2, 542 * 3))  #Res*
    snp_num += 1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "ITHKRRISRTRPRQGLTRER", "20H24M1I15M2I18M20H", 524 * 3 - 2, 542 * 3)) #Res*
    snp_num += 1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "ITHKRRISRTRPRRSDP*R", "20H24M1I29M1D1M20H", 524 * 3 - 2, 542 * 3 - 2))#Sus
    snp_num += 1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "ITHKRRSSALGPGGLTRER", "20H18M1D3M1I34M20H", 524 * 3 - 2, 542 * 3 - 1))#Sus
    snp_num += 1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "ITHKRCSSALGPGGLTRER", "20H18M1D3M1I34M20H", 524 * 3 - 2, 542 * 3 - 1))#Res
    snp_num += 1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "ITHKRCISRTRPRRSDP*T", "20H24M1I32M20H", 524 * 3 - 2, 542 * 3 - 1))   #Res
    snp_num += 1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "ITHKRRISRPRPRRSDP*T", "20H24M1I32M20H", 524 * 3 - 2, 542 * 3 - 1))   #Sus
    snp_num += 1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "ITHKRRISRTRPGGLTRER", "20H3M1D17M1I3M1I9M1D23M20H", 524 * 3 - 2, 542 * 3))   #Res* and 12+ fs
    snp_num += 1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "ITHKRRISRTRPRGLTRER", "20H3M1D17M1I3M1I15M1D17M20H", 524 * 3 - 2, 542 * 3))  #Res* and 12+ fs
    SAM_file.close()
    pysam.sort('-o', 'test/frameshift.bam', 'test/frameshift.sam')

def nonstop_test():
    SAM_file = open("test/nonstop.sam", "w")
    SAM_file.write(header)
    snp_num = 1
    SAM_file.write(make_test(snp_num, full_gene_name[8], "GSAS", "20H12M20H", 479 * 3 - 2, 482 * 3))       #sus
    snp_num += 1
    SAM_file.write(make_test(snp_num, full_gene_name[8], "GSAE", "20H9M1D2M20H", 479 * 3 - 2, 482 * 3))    #res
    snp_num += 1
    SAM_file.write(make_test(snp_num, full_gene_name[8], "GSA*K", "20H10M1I2M20H", 479 * 3 - 2, 482 * 3))  #sus
    snp_num += 1
    SAM_file.write(make_test(snp_num, full_gene_name[8], "GSALK", "20H10M1I2M20H", 479 * 3 - 2, 482 * 3))  #res
    SAM_file.close()
    pysam.sort('-o', 'test/nonstop.bam', 'test/nonstop.sam')
    
def deletion_test6():
    SAM_file = open("test/deletion6.sam", "w")
    SAM_file.write(header)
    snp_num = 1
    SAM_file.write(make_test(snp_num, full_gene_name[8], "AA*T", "20H5M1D7M20H", 24 * 3 - 2, 28 * 3 - 2))    #valid
    snp_num += 1
    SAM_file.write(make_test(snp_num, full_gene_name[8], "AARQ", "20H5M2D7M20H", 24 * 3 - 2, 28 * 3 - 1))    #FStillend
    snp_num += 1
    SAM_file.write(make_test(snp_num, full_gene_name[8], "AA*S", "20H5M1D4M1I2M20H", 24 * 3 - 2, 27 * 3))    #FStillend
    snp_num += 1
    SAM_file.write(make_test(snp_num, full_gene_name[8], "AQ*T", "20H2M1D10M20H", 24 * 3 - 2, 28 * 3 - 2))   #FStillend
    SAM_file.close()
    pysam.sort('-o', 'test/deletion6.bam', 'test/deletion6.sam')

def edge_cases_test():
    SAM_file = open("test/edge_cases.sam", "w")
    SAM_file.write(header)
    snp_num = 1

    # MEG_6094
    #  H   K   R   R   I   S / F   A   L   G   P   G
    # CAC AAA CGT CGT ATC TCC/TTC GCA CTC GGC CCA GGC
    # 526 527 528 529 530   531   532 533 534 535 536

    # MEG_6090
    #  H   K   R   R   L   S / C   A   L   G   P   G
    # CAC AAG CGC CGA CTG TCG/TGC GCG CTG GGG CCC GGC
    # 445 446 447 448 449   450   451 452 453 454 455

    # Case 1: |MMM|MMM|MMM| 
    #             1   1   1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "ATCTTCGCA", "20H9M20H", 530 * 3 - 2, 532 * 3, True))
    snp_num += 1
    # Case 2: |MMM|III|MMM|
    #             1   2   4
    SAM_file.write(make_test(snp_num, full_gene_name[7], "ATCTTCTCC", "20H3M3I3M20H", 530 * 3 - 2, 531 * 3, True))
    snp_num += 1
    # Case 3: |MMM|MII|IMM|MMM|
    #             1       4   1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "ATCTTCTCCGCA", "20H4M3I5M20H", 530 * 3 - 2, 532 * 3, True))
    snp_num += 1
    # Case 4: |MMM|IMM|M II|MMM|
    #             1     4  4   1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "ATCTTCCCCGCA", "20H3M1I3M2I3M20H", 530 * 3 - 2, 532 * 3, True))
    snp_num += 1
    # Case 5: |MMM|IMM|M MI|IMM|MMM|
    #             1     4      4   1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "TCCCGCACGTTCGGC", "20H3M1I4M2I5M20H", 531 * 3 - 2, 534 * 3, True))
    snp_num += 1
    # Case 6: |MMM|III|III|III|MMM|
    #             1   2   4   5   5
    SAM_file.write(make_test(snp_num, full_gene_name[7], "TCCCGCACGTTCGCA", "20H3M9I3M20H", 531 * 3 - 2, 532 * 3, True))
    snp_num += 1
    # Case 7: |MMM|MMI|III|IIM|MMM|
    #             1           6   1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "ATCTCCCGCACCGCA", "20H5M6I4M20H", 530 * 3 - 2, 532 * 3, True))
    snp_num += 1
    # Case 8: |MMM|MMI|III|M MM|IIM|MMM|
    #             1         8      4   1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "ATCTCCCGCCGCCGACTC", "20H5M4I3M2I4M20H", 530 * 3 - 2, 533 * 3, True))
    snp_num += 1
    # Case 9: |MMM|MMI|III|M II|MMM|
    #             1         8  4   1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "ATCTCCCGCCCTGCA", "20H5M4I1M2I3M20H", 530 * 3 - 2, 532 * 3, True))
    snp_num += 1
    # Case 10:  |MMM|DDD DDD DDD MMM|
    #               1   9   9   9   1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "ATCGGC", "20H3M9D3M20H", 530 * 3 - 2, 534 * 3, True))
    snp_num += 1
    # Case 11:  |MMM|MMD M|MM M|MM DDM|MMM|
    #               1   10   11   e   16  1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "ATCTCGCACTCCCCA", "20H5M1D6M2D4M20H", 530 * 3 - 2, 535 * 3, True))
    snp_num += 1
    # Case 12:  |MMM|MMD M|MM M|DD MMM|MMM|
    #               1   10   11   15  1   1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "ATCTCGCACGGCCCA", "20H5M1D4M2D6M20H", 530 * 3 - 2, 535 * 3, True))
    snp_num += 1
    # Case 13:  |MMM|MMD DDM |MMM|
    #               1   1013 14  1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "ATCTTACTC", "20H5M3D4M20H", 530 * 3 - 2, 533 * 3, True))
    snp_num += 1
    # Case 14:  |MMM|DMM DDM|MMM|
    #               1   10  12  1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "ATCCCACTC", "20H3M1D2M2D4M20H", 530 * 3 - 2, 533 * 3, True))
    snp_num += 1
    # Case 15:  |MMM|DMM M|DD MMM|MMM|
    #               1   10   12  1   1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "CGTTCTGCACTC", "20H3M1D3M2D6M20H", 529 * 3 - 2, 533 * 3, True))
    snp_num += 1
    # Case 16:  |MMM|MDM M|MD DMM|MMM|
    #               1   10   11  15  1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "CGTCTATTCGCA", "20H4M1D3M2D5M20H", 528 * 3 - 2, 532 * 3, True))
    snp_num += 1
    # Case 17:  |MMM|MDD DMI|IIM|MMM|
    #               1   1013    4   1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "CGTACTTTCGCA", "20H4M3D1M3I4M20H", 529 * 3 - 2, 532 * 3, True))
    snp_num += 1
    # Case 18:  |MMM|MII|IMD DDM  |MMM|
    #               1       3  13 14  1 
    SAM_file.write(make_test(snp_num, full_gene_name[7], "ATCTTTTCACTC", "20H4M3I1M3D4M20H", 530 * 3 - 2, 533 * 3, True))
    snp_num += 1
    # Case 19: S|SSS|MMM|MMM|MMM|MMM|MMM|SSS|S
    #                   1   1   1   1   1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "GGGGCGTATCTTCGCACTCGGGG", "4S15M4S", 529 * 3 - 2, 533 * 3, True))
    snp_num += 1
    # Case 20: SSS|SMM|IMM|M MM|M MM|M DMM|SSS|S
    #                       4    e    e   16
    SAM_file.write(make_test(snp_num, full_gene_name[7], "GGGGTCTCCCGCACTCGCGGGG", "4S2M1I9M1D2M4S", 530 * 3 - 1, 534 * 3, True))
    snp_num += 1
    # Case 21: SSS|SMM|MDM M|MM M|MM M|MM I|MMS|SSS
    #                     10   11   e    e 4
    SAM_file.write(make_test(snp_num, full_gene_name[7], "GGGGGTACTCCGCACTCCGGGGGG", "4S3M1D10M1I2M4S", 529 * 3 - 1, 534 * 3 - 1, True))
    snp_num += 1
    # Case 22: SSS|SMM|DMM M|MM M|MM I|MMM|SSS|S
    #                     10   11   e 4   1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "GGGGGTTCTCCGCACCTCGGGG", "4S2M1D8M1I3M4S", 529 * 3 - 1, 533 * 3, True))
    snp_num += 1
    # Case 23: S|SSS|MMI|M MM|M MM|M DMS|SSS
    #                     4    e    
    SAM_file.write(make_test(snp_num, full_gene_name[7], "GGGGTCCCGCACTCGGGGG", "4S2M1I7M1D1M4S", 530 * 3 - 2, 534 * 3 - 1, True))
    snp_num += 1
    # Case 24: S|SSS|MMM|III|MMM|MMM|DDD|MMM|SSS|S
    #                   1   2   4   1   9   1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "GGGGATCTTCTCCGCAGGCGGGG", "4S3M3I6M3D3M4S", 530 * 3 - 2, 534 * 3, True))
    snp_num += 1
    # Case 25: S|SSS|MII|III|IMM|MDD DDD DMM|SSS|S
    #                           6  1013 e13 14
    SAM_file.write(make_test(snp_num, full_gene_name[7], "GGGGTTCACCTCCGGCGGGG", "4S1M6I3M6D2M4S", 531 * 3 - 2, 534 * 3, True))
    snp_num += 1
    # Case 26: S|SSS|MII|MM M|III|IMM|DMM DDD DDM|SSS|S
    #                      4         7   10  13  14
    SAM_file.write(make_test(snp_num, full_gene_name[7], "GGGGAAAATCTTTTCCCACGGGG", "4S1M2I3M4I2M1D2M5D1M4S", 530 * 3 - 2, 534 * 3, True))
    snp_num += 1
    # Case 27: S|SSS|MMI|III|M II|MDM M|MD DDM M|MD DMM|SSS|S
    #                         8  4   10  11 e e    n   16
    SAM_file.write(make_test(snp_num, full_gene_name[7], "GGGGCGACACTACCTATCGCTCGGGG", "4S2M4I1M2I1M1D3M3D3M2D2M4S", 528 * 3 - 2, 533 * 3, True))
    snp_num += 1
    # Case 28: S|SSS|MMM|III|III|III|MMD DM|D DDD MMM|SSS|S
    #                   1   2   4   5   3    12  9   1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "GGGGATCGCAGATTTCTCCGGCGGGG", "4S3M9I2M2D1M4D3M4S", 530 * 3 - 2, 534 * 3, True))
    snp_num += 1
    # Case 29: SSS|SMI|M MM|M II|MMM|MII|MM M|IMS|SSS
    #                   e    4  4   1      4
    SAM_file.write(make_test(snp_num, full_gene_name[7], "GGGGGCTATTTCTCCGTGCACATGGGG", "4S1M1I4M2I4M2I3M1I1M4S", 529 * 3 - 1, 533 * 3 - 1, True))
    snp_num += 1
    # Case 30: SSS|SMD M|MM M|MM M|MD DMM|MMM|SSS|S
    #                      10   11   e   16  1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "GGGGGATTTCCGCTCGGCGGGG", "4S1M1D8M2D5M4S", 529 * 3 - 1, 534 * 3, True))
    snp_num += 1
    # Case 31: SSS|SMD M|MD DMM|MMM|IMM|S SS|S
    #                      10  12  1   special
    SAM_file.write(make_test(snp_num, full_gene_name[6], "GGGGACGGTATCTTCGGGG", "4S1M1D2M2D5M1I2M4S", 446 * 3 - 1, 450 * 3 - 1, True))
    snp_num += 1
    # Case 32: SSS|SMD M|MM M|IMM|SSS|S
    #                      10    4
    SAM_file.write(make_test(snp_num, full_gene_name[7], "GGGGGATCTTCCGGGG", "4S1M1D4M1I2M4S", 529 * 3 - 1, 531 * 3, True))
    snp_num += 1
    # Case 33: SSS|SMM|MDD MM|M MM|M IM|SS S|S
    #                     10   11   e
    SAM_file.write(make_test(snp_num, full_gene_name[6], "GGGGAACCGTATCTTGGGG", "4S3M2D6M1I1M4S", 446 * 3 - 1, 450 * 3 - 2, True))
    snp_num += 1
    # Case 34: S|SSS|MMM|IIM|MD DMM|SSS|S
    #                   1      4   16
    SAM_file.write(make_test(snp_num, full_gene_name[7], "GGGGCGTTCATTCGGGG", "4S3M2I2M2D2M4S", 529 * 3 - 2, 531 * 3, True))
    snp_num += 1
    # Case 35: S|SSS|MMI|M MDM|SSS|S
    #                     4   16
    SAM_file.write(make_test(snp_num, full_gene_name[7], "GGGGTTTCGAGGGG", "4S2M1I2M1D1M4S", 531 * 3 - 2, 532 * 3, True))
    snp_num += 1
    # Case 36: S|SSS|MMM|III|MDD MM|M SS|SS
    #                   1   2   3    
    SAM_file.write(make_test(snp_num, full_gene_name[6], "GGGGATCTTCTGCAGGGG", "4S3M3I1M2D3M4S", 449 * 3 - 2, 451 * 3, True))
    snp_num += 1
    # Case 37: S|SSS|MII|MM I|III|MMM|DMM DDD DDM|SSS|S
    #                      4 n   4   1   10  13  14
    SAM_file.write(make_test(snp_num, full_gene_name[7], "GGGGAAAATCTTTTCCCACGGGG", "4S1M2I2M4I3M1D2M5D1M4S", 530 * 3 - 2, 534 * 3, True))
    snp_num += 1
    # MEG_6094 extreme edge case of C insertion suppression
    # Case 38: SSS|SM I|III|M MM|M MM|M MM|M MM|M MM|SSS|S 
    #               C'C TCC'C|GC A|CT C|GG C|CC A|GG
    # Modified:SSS|SM M|II  I|MM M|MM M|MM M|MM M|MMS|SSS
    #                        2    4    1    1    1
    SAM_file.write(make_test(snp_num, full_gene_name[7], "GGGGCCTCCCGCACTCGGCCCAGGGGGG", "4S1M4I15M4S", 531 * 3 - 1, 536 * 3 - 1, True))
    snp_num += 1
    # Case 39: MMI|M MM|III|M MM|M II|IMM|IIM
    #               4        e    e         7
    SAM_file.write(make_test(snp_num, full_gene_name[7], "CGTTATCATCTTCTCTGCTAA", "20H2M1I3M3I4M3I2M2I1M20H", 529 * 3 - 2, 532 * 3, True))
    snp_num += 1
    # Case 40: MMM|DDD DDM MM|M MM|M DMM
    #             1   9   10   11   e   16
    SAM_file.write(make_test(snp_num, full_gene_name[7], "CGTTTCCGCATC", "20H3M5D7M1D2M20H", 528 * 3 - 2, 533 * 3, True))
    snp_num += 1
    # Case 41: MMM|DMM M|MM M|MD DDD  DMM
    #             1   10   11   e  13 13 14
    SAM_file.write(make_test(snp_num, full_gene_name[7], "CGTGTATCTCTC", "20H3M1D7M5D2M20H", 528 * 3 - 2, 533 * 3, True))
    snp_num += 1
    # Case 42: MMM|DMM M|MM DDD DDD M|MM DDM
    #             1   10   11   13 13   e   16
    SAM_file.write(make_test(snp_num, full_gene_name[7], "CACAACGTTCCA", "20H3M1D5M6D3M2D1M20H", 526 * 3 - 2, 532 * 3, True))
    snp_num += 1
    # Case 43: MMM|DMM M|MD DDD DDM M|MM DDM
    #             1   10   11   13 e    e   16
    SAM_file.write(make_test(snp_num, full_gene_name[7], "CACAACGCTCCA", "20H3M1D4M6D4M2D1M20H", 526 * 3 - 2, 532 * 3, True))
    snp_num += 1
    # Case 44: MMM|DMM M|MM MDD|DMM DDM
    #             1   10   11  15  10  12
    SAM_file.write(make_test(snp_num, full_gene_name[7], "CACAACGTCTCC", "20H3M1D6M3D2M2D1M20H", 526 * 3 - 2, 531 * 3, True))
    snp_num += 1

    SAM_file.close()
    pysam.sort('-o', 'test/edge_cases.bam', 'test/edge_cases.sam')

def main():
    # Read config file
    config = configparser.ConfigParser()
    config.read('config.ini')

    # Get variant information
    for gene in SeqIO.parse(config['SOURCE_FILES']['SNP_INFO_FASTA'], 'fasta'):
        # Find index of last pipe ('|') before variant list
        index = -1
        for pipe in range(5):
            index = gene.name[index+1:].find('|') + index + 1
        name = gene.name[:index]
        variants = gene.name[index+1:]

        # Create Gene object
        if 'Must:' in variants:
            if 'Nuc:' in variants:                  gene_variant = Gene.IntrinsicrRNA(name,gene.seq,variants)
            else:                                   gene_variant = Gene.IntrinsicProtein(name,gene.seq,variants)
        elif 'Hyper:' in variants:                  gene_variant = Gene.Hypersusceptible(name,gene.seq,variants)
        elif 'FS-' in variants:
            if 'suppression' in variants:           gene_variant = Gene.Suppressible(name,gene.seq,variants)
            elif 'MEG_6142' in name:                gene_variant = Gene.NormalProtein(name,gene.seq,variants)
            else:                                   gene_variant = Gene.Frameshift(name,gene.seq,variants)
        else:
            if 'Nuc:' in variants:                  gene_variant = Gene.NormalrRNA(name,gene.seq,variants)
            else:                                   gene_variant = Gene.NormalProtein(name,gene.seq,variants)
        gene_dict.update({name + "|RequiresSNPConfirmation":gene_variant})

    for gene in SeqIO.parse(config['SOURCE_FILES']['MEGARES_FASTA'], 'fasta'):
        if gene.name not in gene_dict.keys() and gene.name in header:
            gene.name

    if not os.path.exists('test/'):
        os.makedirs('test/')

    # nonstop_test()
    # frameshift_test()
    # test1()
    # test2()
    # insertion_deletion_test1()
    # insertion_deletion_test2()
    insertion_deletion_test3()
    insertion_deletion_test4()
    # insertion_test1()
    # insertion_test2()
    # insertion_test3()
    # insertion_test4()
    # insertion_test5()
    # insertion_test6()
    # deletion_test1()
    # deletion_test2()
    # deletion_test3()
    # deletion_test4()
    # deletion_test5()
    # deletion_test6()
    # edge_cases_test()

if __name__ == '__main__':
    main()