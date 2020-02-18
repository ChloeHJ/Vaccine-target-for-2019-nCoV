# Vaccine-target-for-2019-nCoV
In silico identification of vaccine targets for 2019-nCoV

This repository contains data, code and supplementary tables for "In silico identification of vaccine targets for 2019-nCoV". 

The data/ contain:
- Open reading frame sequences of 2019-nCoV downloaded from NCBI (.fasta)
- Immunogenic peptides characterized positive by in vitro T cell assays deposited in IEDB (.csv)
- Output data tables from NetMHCpan 4.0 on different HLA alleles (.xlsx)
- Positional weight matrices for predicting immunogenicity against 1G4 CD8+ TCR (.RData)

The src/ contain: 
1. Code for analyzing similarity between 2019-nCoV proteome and previously characterized immunogenic peptides in IEDB;
2. Code for de novo search of 9-mer 2019-nCoV peptides that have MHC presentation and TCR recognition potential by NetMHCpan 4.0 and iPred respectively;
3. Code for de novo search of 9-mer 2019-nCoV peptides that have recognition potential against 1G4 CD8+ TCR molecule by Binding/Activation/Killing positional weight matrices.

The Supplementary Tables/ contain: 
- Table 1. nCoV peptides having exact match with immunogenic SARS CoV peptides
- Table 2. nCoV peptides with high sequence similarity with immunogenic IEDB peptides
- Table 3. de novo search on 9-mer nCoV for immunogenic peptides by NetMHCpan and iPred
- Table 4. de novo search on 9-mer nCoV for immunogenic peptides by NetMHCpan and PWM

