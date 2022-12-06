#!/usr/bin/env skeleton
# -*- coding: utf-8 -*-
#%%
"""
Created on Mon Sep 19 21:30:44 2022

@author: chenqu
"""

import numpy as np
import pandas as pd
import scanpy as sc
import scipy as sp
from scipy import stats
from collections import Counter
from collections.abc import Iterable
import rpy2
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector


#%%
"""
function to compute average gene expression in bins along pseudotime  

adata: cell adata
bin_no: number of bins to be divided along pseudotime
genes: genes for the computation 
pseudotime_col: column in adata.obs where pseudotime is stored

[return]: gene_summary, a dataframe with genes as rows, and pseudotime bins as columns, and averaged gene expression as the data
    
"""
def bin_expression(adata, bin_no, genes, pseudotime_col):
    # define bins
    bins = np.linspace(0, 1, bin_no+1)
    
    # get gene expression
    x = np.array(adata[:,genes].X.todense())
    # get pseudotime
    y = np.array(adata.obs[pseudotime_col])
    
    # calculate average gene expression in each bin
    gene_summary = pd.DataFrame(columns = bins[:-1], index = genes)
    for i in range(gene_summary.shape[1]):
        time = bins[i]
        select = np.array(bins[i] <= y) & np.array(y < bins[i+1])
        gene_summary.loc[:,time]= np.mean(x[select,:],axis=0)
    
    return gene_summary

#%%
"""
function to compute chatterjee correlation of gene expression with pseudotime  

adata: cell adata
genes:genes selected to compute the correlation
pseudotime_col: column in adata.obs where pseudotime is stored

[return]: cor_res, 
    a dataframe with genes as rows, 
    with cor_res (correlation statistics), 
    pval (p-value), 
    adj_pval (p-value adjusted by BH method) as columns 
    
"""

def chatterjee_corr(adata, genes, pseudotime_col, tie_breaking='theoretical'):
    from scipy.stats import rankdata
    
    # get gene expression
    x = np.array(adata[:,genes].X.todense())
    # get pseudotime
    y = list(adata.obs[pseudotime_col])
    
    # compute chatterjee correlation
    # ref: Sourav Chatterjee (2021) A New Coefficient of Correlation, Journal of the American Statistical Association, 116:536, 2009-2022, DOI: 10.1080/01621459.2020.1758115
    
    if tie_breaking == 'theoretical':
        n_cell, n_gene = x.shape
        stat = np.zeros(n_gene)
        x_sorted = x[np.argsort(y), :]
        r = rankdata(x_sorted, method='max', axis=0)
        l = n_cell + 1 - rankdata(x_sorted, method='min', axis=0)
        stat = 1 - n_cell * np.sum(np.abs(np.diff(r, axis=0)), axis=0) / 2 / np.sum(l * (n_cell - l), axis=0)
        '''
        
        for j in range(n_gene):
            print(j, '/', n_gene, end='\r', flush=True)
            x_j = x_sorted[:, j]
            mx = np.reshape(np.repeat(x_j, n_cell) >= np.tile(x_j, n_cell), (n_cell, n_cell))
            r = np.sum(mx, axis=1)
            l = np.sum(mx, axis=0)
            stat[j] = 1 - n_cell * np.sum(np.abs(np.diff(r))) / 2 / np.sum(l * (n_cell - l))
        '''
    elif tie_breaking == 'random':
        
        # add small perturbation for random tie breaking
        np.random.seed(0) # set random seed for consistent results
        x = x + np.random.randn(x.shape[0], x.shape[1]) * 1e-15  
        stat = 1 - np.sum(np.abs(np.diff(rankdata(x[np.argsort(y), :], axis=0), axis=0)), axis=0) * 3 / (x.shape[0] ** 2 - 1)
        stat = np.array(stat).flatten()
    else:
        raise Exception("tie_breaking should be 'theoretical' or 'random'.")
    
    pval = 1 - sp.stats.norm.cdf(stat, loc=0, scale=np.sqrt(2/5/x.shape[0]))
    
    # put results into dataframe cor_res
    cor_res = pd.DataFrame({'cor_stat': stat, 'pval': pval})
    cor_res.index = genes
    
    # compute adjusted pval using BH method
    stats = importr('stats') 
    cor_res.loc[:,'adj_pval']=stats.p_adjust(FloatVector(cor_res.loc[:,'pval']), method = 'BH')

    # sort genes based on adjusted pval
    cor_res= cor_res.sort_values(by='adj_pval')
    
    return cor_res


#%%
"""
match function finds positions of first match for elements in small_list in
    big_list.
small_list - a list, or numpy.ndarray, or pandas.Series object
big_list - a list, or numpy.ndarray, or pandas.Series object
nomatch - value to include if no match is found, -1 by default
[return] - a list of indices of first matches of small_list in big_list
"""
def match(small_list, big_list, nomatch=-1, sortable=True):
    if sortable:
        order = np.array(np.argsort(big_list))
        big_sorted = np.array(big_list)[order]
        small_list = np.array(small_list)
        l = np.searchsorted(big_sorted, small_list, side='left')
        insert_at_last = l == len(order)
        l[insert_at_last] = 0
        ifnomatch = insert_at_last | (big_sorted[l]!=small_list)
        ret = order[l]
        if np.any(ifnomatch):
            ret[ifnomatch] = nomatch
    else:
        ret = np.array([big_list.index(item) for item in small_list])
    return ret

#%%
"""
lookup function matches a vector of `lookup_value` to the `match_col` column of
    `dataframe` and returns the corresponding values in `result_col` colume of
    `dataframe`. 
lookup_value - a vector of values to lookup, can be a panda.Series, 
    numpy.ndarray or list
dataframe - the lookup table/array, can be a pandas.DataFrame, panda.Series or
    numpy.ndarray or list
match_col - column of `dataframe` to use as values to look up from, can be an 
    int or a string, if equal to -1, rownames of `dataframe` are used
result_col - column of `dataframe` to use as the result vector, can be an int
    or a string, if equal to -1, rownames of `dataframe` are used, if equals
    None, the matched ordering is returned.
[return] - a numpy.ndarray of lookup results, with unmatched positions filled
    with NaN.
"""
def lookup(lookup_value, dataframe, match_col=0, result_col=None):
    isIterable = isinstance(lookup_value, Iterable)
    lookup_value = pd.Series(lookup_value) if isIterable else pd.Series([lookup_value])
    dataframe = pd.DataFrame(dataframe)
    if type(match_col) is int:
        tmp = dataframe.iloc[:, match_col] if match_col >= 0 else dataframe.index
    else:
        tmp = dataframe[match_col]
       
    if result_col is None:
        ret = np.array(match(lookup_value, tmp))
        return(ret if isIterable else ret[0])
    elif type(result_col) is int:
        tmp2 = dataframe.iloc[:, result_col] if result_col >= 0 else dataframe.index
    else:
        tmp2 = dataframe[result_col]
    tmp2 = np.append(tmp2, np.nan)
    m = match(lookup_value, tmp, dataframe.shape[0])
    return(tmp2[m] if isIterable else tmp2[m][0])