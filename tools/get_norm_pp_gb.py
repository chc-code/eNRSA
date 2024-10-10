import pandas as pd
import os, sys, re
import numpy as np

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from pydeseq2.preprocessing import deseq2_norm_fit

def run_deseq2(data, metadata, ref_level, col_group=None, min_reads=10, filter_out_low_expression=False, 
               size_factors_in=None, size_factor_only=False):
    """
    A count matrix of shape ‘number of samples’ x ‘number of genes’, containing read counts (non-negative integers),
    Metadata (or “column” data) of shape ‘number of samples’ x ‘number of variables’, containing sample annotations that will be used to split the data in cohorts.  rows = sample, columns = gene
    index should be samplename
    meta data, first col is the sample name (index), following columns usually condidion and group
    min_reads:  if a gene has less than min_reads in all samples, it will be removed
    ref_level,  the level used as reference in the comparison
    """
    
    # data filtering
    
    counttable = data.iloc[:, 2:]
    gene_col = data.iloc[:, 1] # use Transcript as index
    counttable.index = data.iloc[:, 0]
    
    counttable = counttable.T
    metadata.index = counttable.index
    
    # remove samples for which annotations are missing and exclude genes with very low levels of expression
    col_condition = 'condition'
    sample_to_keep = list(~metadata[col_condition].isna())
    counttable = counttable.loc[sample_to_keep, :]
    metadata = metadata.loc[sample_to_keep, :]
    
    if filter_out_low_expression:
        genes_to_keep = counttable.columns[counttable.sum(axis=0) >= min_reads]
        counttable = counttable[genes_to_keep]

    inference = DefaultInference(n_cpus=8)
    dds = DeseqDataSet(
        counts=counttable,
        metadata=metadata,
        design_factors=col_condition, # compare samples based on condition
        ref_level=[col_condition, ref_level],
        refit_cooks=True,
        inference=inference,
    )
    if size_factors_in is None:  # gb change
        # tested, the size factors are exactly the same as the R deseq2 code
        dds.fit_size_factors()
        size_factors = dds.obsm['size_factors']
    else:
        size_factors = size_factors_in
        dds.obsm["size_factors"] = size_factors

    if size_factor_only:
        return None, size_factors
    
    # if each group only have one sample
    n_sam = counttable.shape[0]
    if n_sam == 2:
        log2fc = np.log2(data.iloc[:, 3] / size_factors[1]) / (data.iloc[:, 2] / size_factor[0])
        res_df = data.iloc[:, [0, 1]]
        res_df['log2fc'] = log2fc
        return res_df, size_factors

    if size_factors_in is None:  # gb change
        dds.deseq2()
    else:
        counts = dds.X
        (dds.logmeans, dds.filtered_genes) = deseq2_norm_fit(counts)
        dds.layers["normed_counts"] = counts / size_factors[:, None]

        dds.fit_genewise_dispersions()
        dds.fit_dispersion_trend()
        dds.fit_dispersion_prior()
        dds.fit_MAP_dispersions()
        dds.fit_LFC()
        dds.calculate_cooks()
        dds.refit()

    # avialble attributes = X: count tdata, obs = design factors, 
    # obsm = sample matrix, keys = design_matrix, size_factors
    # varm = gene matrix, keys = dispersions, LFC

    # get stat
    stat_res = DeseqStats(dds, inference=inference)
    stat_res.summary()
    
    # the logFC is the same as R, but the p-value and padj are slightly different
    res_df = stat_res.results_df
    # here must convert gene_col to list, otherwise, the value will be all NaN, because the index does not match
    res_df.insert(0, 'Gene', list(gene_col))
    res_df.reset_index('Transcript', inplace=True)
    res_df = res_df.sort_values('padj')
    
    return res_df, size_factors




def change_pp_gb_with_case(rep1, rep2, data, out_dir, window_size, factor_flag, factor1=None, factor2=None):
    # data = count_pp_gb.txt then collapse transcripts to genes
    # if factor1, factor2 is not None, it should be a list
    
    # columns = "Transcript\tGene"  + ppc_[sam_list], gbc_sam1, gbd_sam1, gbc_sam2, gbd_sam2 ....
    if rep2 == 0:
        factor2 = []
    n_prev_cols = 2 + 4  # 2 = Transcript, Gene, 4 = chr, start, end and strand
    n_sam = rep1 + rep2
    col_idx = {
        'ppc': {
            'cond1': list(range(n_prev_cols, n_prev_cols + rep1)),
            'cond2': list(range(n_prev_cols + rep1, n_prev_cols + n_sam)),
            'combined': list(range(n_prev_cols, n_prev_cols + n_sam)),
        },
        'gbc': {
            'cond1': list(range(n_prev_cols + n_sam,             n_prev_cols + n_sam + rep1 * 2, 2)),
            'cond2': list(range(n_prev_cols + n_sam + rep1 * 2, n_prev_cols + n_sam * 3, 2)),
            'combined': list(range(n_prev_cols + n_sam, n_prev_cols + n_sam * 3, 2)),
        },
        'gbd': {
            'cond1': list(range(n_prev_cols + n_sam + 1,             n_prev_cols + n_sam + rep1 * 2, 2)),
            'cond2': list(range(n_prev_cols + n_sam + rep1 * 2 + 1,  n_prev_cols + n_sam * 3, 2)),
            'combined': list(range(n_prev_cols + n_sam + 1, n_prev_cols + n_sam * 3, 2)),
        },
    }
    
    # cond1 = ctrl, cond2 = case
    idx_ppc_cond1 = col_idx['ppc']['cond1']
    idx_ppc_cond2 = col_idx['ppc']['cond2']
    
    idx_gbc_cond1 = col_idx['gbc']['cond1']
    idx_gbc_cond2 = col_idx['gbc']['cond2']

    idx_gbd_cond1 = col_idx['gbd']['cond1']
    idx_gbd_cond2 = col_idx['gbd']['cond2']
    
    idx_ppc_combined = col_idx['ppc']['combined']
    idx_gbc_combined = col_idx['gbc']['combined']
    idx_gbd_combined = col_idx['gbd']['combined']
    
    print(idx_ppc_combined)
    print(list(enumerate(data.columns)))
    print(col_idx)
    
    sam_list = [re.sub('^ppc_', '', _) for _ in data.columns[idx_ppc_combined]]
    
    data['tmp_ppc_pass'] = data.iloc[:, idx_ppc_combined].apply(lambda x: all(x > 0), axis=1)
    
    for sn, idx_gbc_tmp, idx_gbd_tmp in zip(range(n_sam), idx_gbc_combined, idx_gbd_combined):
        sum_col = data.iloc[:, idx_gbc_tmp].sum()
        data[f'tmp_gbd_pass_{sn+1}'] = data.iloc[:, idx_gbd_tmp] * 10_000_000 / sum_col
    cols_gbd_pass = [f'tmp_gbd_pass_{sn+1}' for sn in range(n_sam)]
    data['tmp_gbd_pass'] = data[cols_gbd_pass].apply(lambda x: x.sum()/n_sam > 0.004, axis=1)
    
    
    # already compared with R result, the same
    data_pass = data[data['tmp_ppc_pass'] & data['tmp_gbd_pass']] # data0 in perl code
    data_drop = data.drop(data_pass.index) # data_0 in perl code

    data_pass = data_pass.drop(columns=[_ for _ in data_pass.columns if _[:4] == 'tmp_'])
    data_drop = data_drop.drop(columns=[_ for _ in data_drop.columns if _[:4] == 'tmp_'])

    data_pass.iloc[:, [0, 1]].to_csv(f'{out_dir}/intermediate/active_gene.txt', sep='\t', index=False, header=False)
    
    data_pass_pp = data_pass.iloc[:, [0, 1, *idx_ppc_combined]] # datapp in R
    data_drop_pp = data_drop.iloc[:, [0, 1, *idx_ppc_combined]] # data_pp in R
    
    data_pass_gb = data_pass.iloc[:, [0, 1, *idx_gbc_combined]] # datagb in R
    data_drop_gb = data_drop.iloc[:, [0, 1, *idx_gbc_combined]] # data_gb in R
    
    if rep2 == 0:
        # case only, the meta data is just the dummy value, won't be used for deseq, just for get the sizefactors
        metadata = pd.DataFrame({
            'condition': ['control'] + ['case'] * (rep1 - 1),
        })
    else:
        metadata = pd.DataFrame({
            'condition': ['control'] * rep1 + ['case'] * rep2,
        })
    ref_level = 'control'

    # with normalization factors
    if factor_flag == 1:
        factor1 = np.array(factor1)
        factor2 = np.array(factor2)
        norm_factors = np.concatenate([factor1, factor2])
        size_factors = 1 / norm_factors

    else:  # without normalization factors
        size_factors = None
    
    if rep2 > 0:
        # gb change
        res_df, size_factors = run_deseq2(data_pass_gb, metadata, ref_level, size_factors_in=size_factors)
        res_df_full = pd.concat([res_df, data_drop_gb.iloc[:, [0, 1]]])
        res_df_full.to_csv(f'{out_dir}/known_gene/gb_change.txt', sep='\t', index=False)
        
        # pp change
        res_df, _ = run_deseq2(data_pass_pp, metadata, ref_level=ref_level, size_factors_in=size_factors)
        res_df_full = pd.concat([res_df, data_drop_pp.iloc[:, [0, 1]]])
        res_df_full.to_csv(f'{out_dir}/known_gene/pp_change.txt', sep='\t', index=False)
    else:
        size_factors = 1 # size factor is 1 for the only sample

    # get the normalized data
    norm_factors = 1 / size_factors
    ppc_norm = data.iloc[:, idx_ppc_combined] * norm_factors
    ppd_norm = ppc_norm / window_size
    ppd_norm.columns = [f'ppd_{_}' for _ in sam_list]
    gbc_norm = data.iloc[:, idx_gbc_combined] * norm_factors
    gbd_norm = data.iloc[:, idx_gbd_combined] * norm_factors
    data_normalized = pd.concat([data.iloc[:, :n_prev_cols], ppc_norm, ppd_norm, gbc_norm, gbd_norm], axis=1)
    data_normalized.to_csv(f'{out_dir}/known_gene/normalized_pp_gb.txt', sep='\t', index=False)
    
    # save normalization factors
    nf = pd.DataFrame({'sample': sam_list, 'nfactor': norm_factors})
    nf.to_csv(f'{out_dir}/intermediate/nf.txt', sep='\t', index=False)

out_dir = '/Users/files/tmp1/nrsa/demo_data'
data = pd.read_csv('/Users/files/tmp1/nrsa/demo_data/intermediate/count_pp_gb.txt', sep='\t')
window_size = 50
n_gene_cols = 2


change_pp_gb_with_case(2, 2, data, out_dir, window_size, 0)
