import pandas as pd
import numpy as np
from collections import defaultdict
from scipy.stats import hypergeom
from typing import Callable

import databases as db


def outlier_enrichment(df : pd.DataFrame, sign : str, func : Callable[pd.DataFrame, dict], multiplicative : str, tau_mult_sig_val : dict) -> tuple[dict|int, dict|int, dict|int, dict|int]:
    """
    This generic function passes a subselection of outliers to some external function, e.g. 'fraction_three_shared',
    categorizing outliers as significant for both chimeric and multiplicative scale, only chimeric (reported), or only multiplicative.
    The selection of data depends on whether you're interested in negative or positive outliers.
    """
    if sign == "positive":
        overlap = func(df[(df.adjusted_interaction_score_epsilon_or_tau > 0.08) & (df[multiplicative] > tau_mult_sig_val['pos'])])
        only_reported = func(df[(df.adjusted_interaction_score_epsilon_or_tau > 0.08) & (df[multiplicative] < tau_mult_sig_val['pos'])])
        only_mult = func(df[(df.adjusted_interaction_score_epsilon_or_tau < 0.08) & (df[multiplicative] > tau_mult_sig_val['pos'])])
        all_reported = func(df[df.adjusted_interaction_score_epsilon_or_tau > 0.08])
        all_mult = func(df[df[multiplicative] > tau_mult_sig_val['pos']])
    elif sign == "negative":
        overlap = func(df[(df.adjusted_interaction_score_epsilon_or_tau < -0.08) & (df[multiplicative] < tau_mult_sig_val['neg'])])
        only_reported = func(df[(df.adjusted_interaction_score_epsilon_or_tau < -0.08) & (df[multiplicative] >= tau_mult_sig_val['neg'])])
        only_mult = func(df[(df.adjusted_interaction_score_epsilon_or_tau >= -0.08) & (df[multiplicative] < tau_mult_sig_val['neg'])])   
        all_reported = func(df[df.adjusted_interaction_score_epsilon_or_tau < -0.08])
        all_mult = func(df[df[multiplicative] < tau_mult_sig_val['neg']])

    return overlap, only_reported, only_mult, all_reported, all_mult


def outlier_enrichment_both_tails(df, func, multiplicative, tau_mult_sig_val):
    """
    This generic function passes a subselection of outliers to some external function,
    categorizing outliers as significant for both cumulant and multiplicative scale, only cumulant (reported), or only multiplicative.
    The selection of data depends on whether you're interested in negative or positive outliers.
    """

    all_reported = func(df[(df.adjusted_interaction_score_epsilon_or_tau > 0.08) | 
                                (df.adjusted_interaction_score_epsilon_or_tau < -0.08)])
    all_mult = func(df[(df[multiplicative] > tau_mult_sig_val['pos']) | 
                        (df[multiplicative] < tau_mult_sig_val['neg'])])
    
    return all_reported, all_mult


def fraction_physical_twoplus(df : pd.DataFrame) -> dict:
    """
    Returns a dict of the number of gene triplets with 2+ physical interactions along with the total number of triplets
    """
    d = {}
    d['tot'] = len(df.twoplus_physical_interactions)
    d['int'] = np.sum(df.twoplus_physical_interactions)
    d['frac'] = d['int']/d['tot']
    return d

def fraction_physical_three(df : pd.DataFrame) -> dict:
    """
    Returns a dict of the number of gene triplets with 3 physical interactions along with the total number of triplets
    """
    d = {}
    d['tot'] = len(df.three_physical_interactions)
    d['int'] = np.sum(df.three_physical_interactions)
    d['frac'] = d['int']/d['tot']
    return d

def fraction_three_shared(df : pd.DataFrame) -> dict:
    """
    Returns a dict of the number of gene triplets in which all 3 share at elast 1 interactor
    """
    d = {}
    d['tot'] = len(df.three_shared_physical_interactions)
    d['int'] = np.sum(df.three_shared_physical_interactions)
    d['frac'] = d['int']/d['tot']
    return d

def fraction_coex_twoplus(df : pd.DataFrame) -> dict:
    """
    Returns a dict of the number of gene triplets with 3 coexpresison interactions along with the total number of triplets
    """
    d = {}
    d['tot'] = len(df.twoplus_coex_interactions)
    d['int'] = np.sum(df.twoplus_coex_interactions)
    d['frac'] = d['int']/d['tot']
    return d


def fraction_coex_oneplus(df : pd.DataFrame) -> dict:
    """
    Returns a dict of the number of gene triplets with 3 coexpresison interactions along with the total number of triplets
    """
    d = {}
    d['tot'] = len(df.oneplus_coex_interactions)
    d['int'] = np.sum(df.oneplus_coex_interactions)
    d['frac'] = d['int']/d['tot']
    return d

def fraction_coex_three(df : pd.DataFrame) -> dict:
    """
    Returns a dict of the number of gene triplets with 3 coexpresison interactions along with the total number of triplets
    """
    d = {}
    d['tot'] = len(df.three_coex_interactions)
    d['int'] = np.sum(df.three_coex_interactions)
    d['frac'] = d['int']/d['tot']
    return d

def fraction_threeway_shared_go(df : pd.DataFrame) -> dict:
    """
    Returns a dict of the number of gene triplets with 3 shared go terms along with the total number of triplets
    """
    d = {}
    d['tot'] = len(df.threeway_shared_go)
    d['int'] = np.sum(df.threeway_shared_go)
    d['frac'] = d['int']/d['tot']
    return d

def get_hypergeom_params(sample : dict, population : dict) -> list[int]:
    """
    Returns parameters for hypergeom.cdf in a list [k, M, n, N]

    Hypergeometric distribution parameters according to Wikipedia (with scipy.stats hypergeom parameters in parentheses)
    N (M) is the population size,
    K (n) is the number of success states in the population,
    n (N) is the number of draws (i.e. quantity drawn in each trial),
    k (k) is the number of observed successes,
    hypergeom.cdf(k, M, n, N, loc=0)

    """
    assert sample['tot'] <= population['tot'], "you've mis-used the hypergeom function"
    return [sample['int'], population['tot'], population['int'], sample['tot']]

def perform_hypergeom_test(df, sign, func, multiplicative, tau_mult_sig_val):
    """
    This function automates using hypergeom.cdf for both positive and negative outliers, specified by 'sign' parameter.
    The 'func' parameter is used to categories gene triplets by physical interactions or GO categories and must return
    a dict with keys 'int' and 'tot'.
    """
    genome_wide = func(df)
    overlap, only_reported, only_mult, all_reported, all_mult = outlier_enrichment(df, sign, func, multiplicative, tau_mult_sig_val)
    # since we are looking for positive enrichments above population baseline, take 1 - CDF(x)
    overlap_htest = 1-hypergeom.cdf(*get_hypergeom_params(overlap, genome_wide))
    only_reported_htest = 1-hypergeom.cdf(*get_hypergeom_params(only_reported, genome_wide))
    only_mult_htest = 1-hypergeom.cdf(*get_hypergeom_params(only_mult, genome_wide))

    all_reported_htest = 1-hypergeom.cdf(*get_hypergeom_params(all_reported, genome_wide))
    all_mult_htest = 1-hypergeom.cdf(*get_hypergeom_params(all_mult, genome_wide))

    results = pd.DataFrame.from_dict({"type" : ['overlap', 'only_reported', 'only_mult', 'all_reported', 'all_mult'],
                                    "pval" : [overlap_htest, only_reported_htest, only_mult_htest, all_reported_htest, all_mult_htest]})
    return results


def alleles_2_go_enrichment(df):
    """
    This function calculates number of instances in which all three genes belong to the same GO category,
    specified as 'int' key (short for interaction) of a dict
    """
    gene_2_go, goid_2_term = db.get_go_info()
    d = {}
    #go_hit_2plus = 0
    go_hit_3x = 0
    for i,r in df.iterrows():
        alleles = sorted(r['alleles'].split(","))
        alleles = [db.gene_stem_name(i.upper()) for i in alleles]

        go_counts = defaultdict(int)
        for a in alleles:
            if a in gene_2_go:
                # many genes are involved in many GO categories; iterate through these
                for go in gene_2_go[a]:
                    go_counts[go] += 1
    
        counts = np.array([i[1] for i in go_counts.items()])
        #print(np.max(counts))
        if len(counts) > 0:
            assert np.max(counts) <= 3
        #if np.sum(np.any(counts > 1)):
        #    go_hit_2plus += 1
        #if np.sum(np.any(counts == 3)):
        if np.sum(counts == 3) >= 1:
            go_hit_3x += 1

    d['tot'] = len(df)
    d['int'] = go_hit_3x
    d['frac'] = d['int']/d['tot'] 
    return d