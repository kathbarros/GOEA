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

