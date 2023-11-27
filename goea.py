import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

background_file = 'background.csv'
target_genes_file = 'subset.csv'

background_df = pd.read_csv(background_file)
target_genes_df = pd.read_csv(target_genes_file)

# extract background and target genes
background_genes = background_df['Gene']
target_genes = target_genes_df['Gene']

# perform gene ontology enrichment analysis
results = []

# each unique GO term in the background
for go_id in background_df['go_id'].unique():
    # count the number of genes in the background with the current GO term
    background_count = background_df[background_df['go_id'] == go_id]['Gene'].nunique()

    # count the number of genes in the target genes with the current GO term
    target_count = target_genes_df[target_genes_df['go_id'] == go_id]['Gene'].nunique()

    # count the number of genes in the background without the current GO term
    background_remain = background_genes.nunique() - background_count

    # count the number of genes in the target genes without the current GO term
    target_remain = target_genes.nunique() - target_count

    # Fisher's exact test
    oddsratio, p_value = fisher_exact([[target_count, target_remain], [background_count, background_remain]])

    results.append({
        'go_id': go_id,
        'enrichment_fold': (target_count / target_genes.nunique()) / (background_count / background_genes.nunique()),
        'p_value': p_value
    })

results_df = pd.DataFrame(results)

# multiple testing correction (Benjamini-Hochberg)
adjusted_p_values = multipletests(results_df['p_value'], method='fdr_bh')[1]
results_df['adjusted_p_value'] = adjusted_p_values

# genes related to each GO term
results_df['genes'] = ''
for i, go_id in enumerate(results_df['go_id']):
    genes = target_genes_df[target_genes_df['go_id'] == go_id]['Gene'].tolist()
    results_df.at[i, 'genes'] = ','.join(genes)

results_df = results_df.merge(target_genes_df[['go_id', 'goName', 'ontology']], on='go_id')

# filter and include only significant results
significance_threshold = 0.05
significant_results_df = results_df[results_df['adjusted_p_value'] <= significance_threshold]

output_file = 'down_go_enrichment.csv'
significant_results_df.to_csv(output_file, index=False)

