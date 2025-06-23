#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np
from scipy.stats import pearsonr



methylation_df = pd.read_csv('/Users/humairanoor/Documents/AIMP/TCGAGBM/Analysis_files/GBM_methsurv_AIMP_with_ID.csv',index_col=0)
gene_expression_df = pd.read_csv('/Users/humairanoor/Documents/AIMP/TCGAGBM/Analysis_files/GBM_rnaseq_AIMP_with_ID.csv', index_col=0)


overlapping_samples = methylation_df.index.intersection(gene_expression_df.index)

methylation_overlap = methylation_df.loc[overlapping_samples]
gene_expression_overlap = gene_expression_df.loc[overlapping_samples]

correlation_results = pd.DataFrame(index=methylation_overlap.columns, columns=gene_expression_overlap.columns)
pvalue_results = pd.DataFrame(index=methylation_overlap.columns, columns=gene_expression_overlap.columns)


for cpg in methylation_overlap.columns:
    for gene in gene_expression_overlap.columns:
        methylation_values = methylation_overlap[cpg]
        gene_expression_values = gene_expression_overlap[gene]
        
        valid_idx = methylation_values.notna() & gene_expression_values.notna()
        methylation_values = methylation_values[valid_idx]
        gene_expression_values = gene_expression_values[valid_idx]
        
        if len(methylation_values) >= 2:
            corr, p_value = pearsonr(methylation_values, gene_expression_values)
        else:
            corr, p_value = np.nan, np.nan  
        
        
        correlation_results.loc[cpg, gene] = corr
        pvalue_results.loc[cpg, gene] = p_value


print("Correlation Results:")
print(correlation_results)

print("\nP-Value Results:")
print(pvalue_results)

