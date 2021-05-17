import pandas as pd
from statsmodels.stats.multitest import multipletests
from scipy.stats import fisher_exact

inputDir = 'inputs'
resultDir = 'results'

def identify_discreps(version, filter_total_count = 10):
    assert version in ['hg19','hg38'], "Version should be either hg19 (as GRCh37) or hg38 (as GRCh38)"

    # read in count of distinct unique (discordant) variants and concordant variants
    unique_binCount = pd.read_csv(f'{inputDir}/total_unique.{version}.bed.binCount', 
        sep = '\t', low_memory=False, header = None, names= ['chr','start', 'end', 'count_unique'])
    total_binCount = pd.read_csv(f'{inputDir}/total.{version}.bed.binCount',
        sep = '\t', low_memory=False, header = None, names= ['count_total'], usecols= [3])

    binCount = pd.concat([unique_binCount, total_binCount], axis = 1)
    binCount.index = binCount.apply(lambda x: ':'.join( [ str(x['chr']),str(x['start']),str(x['end'])]), axis = 1)
    binCount['non_unique'] = binCount['count_total'] - binCount['count_unique']
    binCount['non_unique'] = binCount.apply(lambda x: x['non_unique'] if x['non_unique'] >0 else 0, axis = 1)
    
    total_unique = binCount['count_unique'].sum()
    total_nonunique = binCount['non_unique'].sum()

    # filter to keep bins with greater than <filter_total_count> variants
    binCount_use = binCount[(binCount['count_total'] > filter_total_count)].copy()
    binCount_use['result'] = binCount_use.apply(
        lambda x: fisher_exact([ [x['count_unique'],x['non_unique']], 
                                 [total_unique, total_nonunique]], 
                  alternative = 'greater'), axis = 1)
    binCount_use['oddsratio'] = binCount_use.apply(lambda x: x['result'][0], axis = 1)
    binCount_use['pvalue'] = binCount_use.apply(lambda x: x['result'][1], axis = 1)
    binCount_use['qValue'] = multipletests(list(binCount_use['pvalue']), method = 'fdr_bh')[1]

    binCount_use_sig = binCount_use.query('qValue < 0.01')
    binCount_use_sig.to_csv(f'{resultDir}/sig_{version}_10kb.tsv', sep = '\t')

if __name__ == "__main__":
    identify_discreps('hg19')
    identify_discreps('hg38')
