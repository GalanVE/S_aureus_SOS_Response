"""
Functional annotation and enrichment analysis of S. aureus co-expression
modules (WGCNA network).

For each module obtained from the WGCNA analysis, this script:
  1. Annotates genes as biofilm-related and/or long non-coding RNA (lncRNA).
  2. Tests each module for a statistically significant over-representation
     of biofilm-related genes (hypergeometric test).
  3. Generates summary plots of module sizes and biofilm enrichment.

Input files expected in the working directory:
  - ModulesStaphy-RPKM-Signed-bicor.csv  (gene -> module color)
  - Microarray_RNAseq.txt                            (biofilm-related genes)
  - NCBIncRNA.txt                                     (lncRNA genes)
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
from collections import Counter
import csv
from scipy.stats import hypergeom
import math

# -----------------------------------------------------------------------------
# 1. Load the gene-to-module assignment table
# -----------------------------------------------------------------------------
modules = pd.read_csv("../Data/ModulesStaphy-RPKM-Signed-bicor.csv", sep='\t')
print(modules.head())
print(len(modules))

# -----------------------------------------------------------------------------
# 2. Flag genes previously reported as biofilm-related
# -----------------------------------------------------------------------------
Biofil = pd.read_csv("../Data/Microarray_RNAseq.txt", sep='\t')
print(Biofil.head())
print(len(Biofil))

modules["BioFilm"] = np.nan
modules['BioFilm'] = modules['BioFilm'].astype(object)
for row in modules.itertuples():
    for row1 in Biofil.itertuples():
        if row.NameGen == row1.NCBI:
            modules.loc[row.Index, 'BioFilm'] = 1

# -----------------------------------------------------------------------------
# 3. Flag genes annotated as long non-coding RNAs (lncRNA)
# -----------------------------------------------------------------------------
#lncRNA = pd.read_csv("NCBIncRNA.txt", sep='\t')
#print(lncRNA.head())
#print(len(lncRNA))

#modules["lncRNA"] = np.nan
#modules['lncRNA'] = modules['lncRNA'].astype(object)
#for row in modules.itertuples():
#    for row1 in lncRNA.itertuples():
#        if row.NameGen == row1.NCBI:
#            modules.loc[row.Index, 'lncRNA'] = 1

#print(modules.head())
#print(modules.count())

# Save the gene table annotated with biofilm/lncRNA information
modules.to_csv('Modulesbiofilm-edgardo.txt', sep='\t', index=False)

# -----------------------------------------------------------------------------
# 4. Summarize module size and annotation counts, per module
# -----------------------------------------------------------------------------
colores = np.unique(modules.Mcolors)
print(len(colores))

d = pd.DataFrame()
for i in colores:
    df_filtered = modules[(modules.Mcolors == i)]
    tamano = df_filtered.index.size
    Biofilmsin = df_filtered['BioFilm'].count()
    lncRNAsin = df_filtered['lncRNA'].count()
    row = {'Modules': i, 'Size': tamano, 'Biofilm': Biofilmsin, 'lncRNA': lncRNAsin}
    d = pd.concat([d, pd.DataFrame([row])])

print(d)

# -----------------------------------------------------------------------------
# 5. Hypergeometric enrichment test: biofilm-related genes per module
# -----------------------------------------------------------------------------
d["Pvalue-Biofilm"] = np.nan
N = 3025  # total number of genes in the genome
m = len(Biofil.NCBI.unique())  # total biofilm-related genes in the genome
for index1, row in d.iterrows():
    x = d.Biofilm[index1]  # biofilm-related genes found in the module
    k = d.Size[index1]     # module size
    p = hypergeom.sf((x - 1), N, k, m)
    d.loc[index1, 'Pvalue-Biofilm'] = p

print(d)

d["Biofilm_logP_value"] = np.nan
d["Biofilm_logP_value"] = d["Pvalue-Biofilm"].apply(lambda x: -math.log10(x))
print(d)

# Keep only modules with a significant enrichment (-log10(p) > 1.3, i.e. p < 0.05)
mask = d['Biofilm_logP_value'] > 1.3
dfil = d[mask][['Modules', 'Biofilm_logP_value']].copy()
dfil = dfil.sort_values(by=['Biofilm_logP_value'], ascending=False)
print(dfil)

# -----------------------------------------------------------------------------
# 6. Biofilm enrichment plot (exploratory version, auto-colored)
# -----------------------------------------------------------------------------
ax = dfil.set_index('Modules')[['Biofilm_logP_value']].plot(
    kind='bar', figsize=(30, 8), color=['#EE3224'], alpha=0.5, grid=0, fontsize=30)
plt.xticks(rotation=60)
plt.ylabel("-log(Pvalue)", fontsize=35)
plt.xlabel("Modules", fontsize=35)
plt.legend(["Biofilm Genes"], fontsize=30)

for p in ax.patches:
    dat = p.get_height()
    if dat > 1.3:
        ax.annotate(str(''), (p.get_x() + .1, p.get_height()), xytext=(0, 30),
                    textcoords='offset points', ha="center", va="center",
                    arrowprops=dict(facecolor='black', shrink=0.1))

plt.savefig('EnrichmentBiofilmStaphy.png', bbox_inches='tight', dpi=300)

# -----------------------------------------------------------------------------
# 7. Biofilm enrichment plot (final figure, module-specific colors/labels)
# -----------------------------------------------------------------------------
color = ['#000000', '#FF00FF', '#90EE90', '#4169E1']
etiq = ['Black', 'Magenta', 'Lightgreen', 'Royalblue']

fig, ax = plt.subplots(figsize=(15, 3))
ax.bar(dfil['Modules'], dfil['Biofilm_logP_value'], color=color)
plt.yticks(fontsize=16)
plt.xticks(dfil['Modules'], etiq, rotation=60, fontsize=16)
plt.ylabel("-log(Pvalue)", fontsize=20)
plt.xlabel("Modules", fontsize=20)
plt.legend('', frameon=False)

plt.savefig('EnrichmentBiofilmStaphy.png', bbox_inches='tight', dpi=300)

# -----------------------------------------------------------------------------
# 8. Module size plot (number of genes per module)
# -----------------------------------------------------------------------------
color = ['black', 'brown', 'cyan', 'darkgreen', 'darkred', 'darkturquoise',
         'green', 'lightgreen', 'magenta', 'royalblue']
fig, ax = plt.subplots(figsize=(12, 6))
ax.bar(d['Modules'], d['Size'], color=color)

plt.yticks(fontsize=16)
plt.xticks(rotation=60, fontsize=16)
plt.ylabel("No. genes", fontsize=20)
plt.xlabel("Modules", fontsize=20)
plt.legend('', frameon=False)

plt.savefig('ModulesStaphy.png', bbox_inches='tight', dpi=400)
