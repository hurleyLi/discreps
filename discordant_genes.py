import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

gene_count_all = pd.read_csv('inputs/var_count_per_gene.tsv', sep = '\t')

# Filter to keep genes with at least one discordant variant call
gene_count = gene_count_all.query('unique_hg38 > 0 or unique_hg19 > 0').copy()
gene_count_gene = list(gene_count['Gene'])

"""
Number of total concordant / discordant variants across filtered genes:
We did not simply add up the values from the gene_count dataframe, because
there are situations where a variant belong to multiple genes. Therefore we 
counted the number of distinct concordant / discordant variants from raw files,
which will be avaiable on dbGaP.
"""

TOTAL_CONCORDANT = 70156
TOTAL_UNIQUE_HG19 = 6660
TOTAL_UNIQUE_HG38 = 7999

# Fisher exact test for each gene
gene_count['fisher_hg19'] = gene_count.apply(lambda x: fisher_exact(
    [[x['unique_hg19'], x['count_concord']],[TOTAL_UNIQUE_HG19, TOTAL_CONCORDANT]], 
    alternative= 'greater')[1], axis = 1)
gene_count['fisher_hg38'] = gene_count.apply(lambda x: fisher_exact(
    [[x['unique_hg38'], x['count_concord']],[TOTAL_UNIQUE_HG38, TOTAL_CONCORDANT]], 
    alternative= 'greater')[1], axis = 1)

# adjust for multiple testing
gene_count['fisher_hg19_qValue'] = multipletests(list(gene_count['fisher_hg19']), method = 'fdr_bh')[1]
gene_count['fisher_hg38_qValue'] = multipletests(list(gene_count['fisher_hg38']), method = 'fdr_bh')[1]

# get genes with significant enrichment
gene_count['sig_hg19'] = gene_count.apply(
    lambda x: True if x['fisher_hg19_qValue']<0.05 and x['unique_hg19'] > 5 else False, axis = 1)

gene_count['sig_hg38'] = gene_count.apply(
    lambda x: True if x['fisher_hg38_qValue']<0.05 and x['unique_hg38'] > 5 else False, axis = 1)

gene_count_sig = gene_count.query('sig_hg19 == True or sig_hg38 == True')
gene_count_sig.to_csv('results/discordant_genes.tsv', sep = '\t', index = False)
