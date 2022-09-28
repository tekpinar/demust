#!/bin/bash

demust maps -i normPred_evolCombi.txt -d gemme -o pyrin-gemme-aa-414-781-coolwarm-combi --offset 413 -b 1 -e 368
demust maps -i normPred_evolEpi.txt -d gemme -o pyrin-gemme-aa-414-781-coolwarm-epi --offset 413 -b 1 -e 368
demust maps -i normPred_evolInd.txt -d gemme -o pyrin-gemme-aa-414-781-coolwarm-ind --offset 413 -b 1 -e 368

#Highlight the wild type amino acids with small dots. 
demust maps -i normPred_evolCombi.txt -d gemme -o pyrin-gemme-aa-414-781-coolwarm-combi-with-wt --offset 413 -b 1 -e 368 -s rcsb_pdb_4CG4-414-781.fasta


#Parse Riesselman, 2016 dataset.
#demust riesselman -d dataRiesselman/Johannessen2016.csv -e SCH_Average -o MK01_Johannessen2016_SCH_Average
#You need to add the dataset to the package. 

#Some example usages of common demust functionalities
#Remove gaps introduced during the alignment process.
#Only the gap column in the reference sequence are removed.
demust removegaps -i PF00397_rp15-aligned.fasta -o PF00397_rp15-aligned-nogaps.fasta


#Now, let's obtain the necessary experimental data in a single line format from Riesselman dataset.
demust riesselman -d Fields2012_singles.csv -e linear --otype dat -o res_Fields2012_singles_singleline.csv

#Convert a matrix type single point mutation data to a single line simple format.
#Each line contains a mutation and its effect
#Example:
#M14A -0.4567
#M15A 0.003
demust convert -i YAP1_normPred_evolCombi.txt --itype gemme -o YAP1_normPred_evolCombi_dms.txt --otype dms -f YAP1.fasta -b 170
#Here we added -b 170 parameter because the first amino acid corresponding to the experimental fasta is aa 171 (yes 171, not 170).


#After converting your experimental and theoretical data to single line format, you can get Spearman correlation for all single point
#mutations in the experimental file res_Fields2012_singles_singleline.csv
demust compare -i res_Fields2012_singles_singleline.csv --itype singleline -j YAP1_normPred_evolCombi_dms.txt --jtype singleline
