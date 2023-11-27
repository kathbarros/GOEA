GO Enrichment Analysis Tool

A Python script for Gene Ontology (GO) enrichment analysis to identify significantly overrepresented GO terms in a target gene set compared to a background gene set.

#Installation

Bash
pip install pandas scipy statsmodels (Seabold S, Perktold J. Statsmodels: Econometric and Statistical Modeling with Python. In: Proceedings of the 9th Python in Science Conference 2010 Jun 28 (Vol. 57, No. 61, pp. 10-25080).

#Usage

Input Data

Provide two CSV files:

• background_file.csv: Contains the background gene set, with each row representing a gene
• target_genes_file.csv: Contains the target gene set, with each row representing a gene
Output

A CSV file named enrichment_results.csv is generated, containing the following columns:

• GO term: The GO term identifier
• Enrichment fold: The fold enrichment of the GO term in the target gene set
compared to the background gene set
• Adjusted p-value: The adjusted p-value using Benjamini-Hochberg (FDR)
correction for multiple testing
• Associated genes: A comma-separated list of target genes associated with the
GO term
• GO description: The GO term description

Script Breakdown

1 Import Libraries: Python

This Python script is specifically designed for Gene Ontology (GO) enrichment analysis of non- conventional species. It utilizes a modified pipeline that first merges gene names with their corresponding GO terms before performing the enrichment analysis. This step is crucial for non-
conventional species due to the potential for missing or incomplete GO annotations.

import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multiple_test

2 Read Input Data: Python

background_df = pd.read_csv('background_file.csv')
target_genes_df = pd.read_csv('target_genes_file.csv')

3 Prepare Data for Enrichment Analysis: Python

background_genes = set(background_df['gene_id'].tolist())
target_genes = set(target_genes_df['gene_id'].tolist())
results = []

4 Iterate through GO Terms: Python

for go_term in set(background_df['GO_term'].tolist()):
    
    # Count genes in each category
    
    background_with_term =
len(background_df[background_df['GO_term'] == go_term])
    target_with_term =
len(target_genes_df[target_genes_df['GO_term'] == go_term])
    background_without_term = len(background_genes) -
background_with_term
    target_without_term = len(target_genes) -
target_with_term
    
    # Perform Fisher's exact test and store results
    
    odds_ratio, p_value =
stats.fisher_exact([background_with_term,
background_without_term], [target_with_term,
target_without_term])
    results.append({
        'GO_term': go_term,
        'odds_ratio': odds_ratio,
        'p_value': p_value
})

5 Create Results Dataframe: Python

results_df = pd.DataFrame(results)

6 Apply Multiple Testing Correction: Python

_, adjusted_pvalues, _, _ =
multiple_test(results_df['p_value'], method='fdr_bh')
results_df['adjusted_p_value'] = adjusted_pvalues

7 Identify Associated Genes: Python

for index, row in results_df.iterrows():
    go_term = row['GO_term']
    associated_genes =
target_genes_df[target_genes_df['GO_term'] == go_term]
['gene_id'].tolist()
    results_df.at[index, 'associated_genes'] = ',
'.join(associated_genes)

8 Merge GO Annotations: Python

results_df = results_df.merge(target_genes_df[['GO_term',
'GO_description']], on='GO_term', how='left')

9 Filter and Save Significant Results: Python

significance_threshold = 0.05
significant_results_df =
results_df[results_df['adjusted_p_value'] <=
significance_threshold]
significant_results_df.to_csv('enrichment_results.csv',
index=False)
