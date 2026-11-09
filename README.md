**Exploring Gene Co-expression in Staphylococcus aureus: Biofilm and Genotoxic Stressors Converge on Modules Associated with the SOS Response
**

Edgardo Galan-Vasquez, Roberta dos Reis Ribeiro, Maiana de Oliveira Cerqueira e Costa, Silvia Tenorio-Salgado, Cinthia Valentina Soberanes-Gutiérrez, Marisa Fabiana Nicolás, and Ernesto Perez-Rueda

Code used to build and analyze the WGCNA gene co-expression network for S. aureus, and to test co-expression modules for enrichment in biofilm-related genes.

Contents
WGCNA_S_aureus.R – Builds the co-expression network from RPKM expression data (quality control, soft-threshold selection, module detection, hub gene identification, and export to Cytoscape).
Functional_annotation.py – Annotates genes per module (biofilm-related genes, lncRNAs) and tests each module for biofilm-gene enrichment (hypergeometric test), with summary plots.
Requirements

R: WGCNA, flashClust, doParallel, rrcov

Python: pandas, numpy, matplotlib, seaborn, scipy

Data

Place the required input files (expression matrix, module/annotation tables) in the same folder as the scripts before running them. See the comments at the top of each script for the expected file names.
