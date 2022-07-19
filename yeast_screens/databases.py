
import pandas as pd
import numpy as np
from collections import defaultdict


def gene_stem_name(gene):
    return gene.split('-')[0]

def find_unique_interactions(df, left_gene, right_gene):
    """
    Assumes a dataframe has 2 columns, one for each of the physically interacting genes
    """
    gene_physical_pairwise_interactions = set()
    for genes in zip(df[left_gene], df[right_gene]):
        gene_physical_pairwise_interactions.add( tuple(sorted((gene_stem_name(genes[0].upper()), gene_stem_name(genes[1].upper())))) )

    return gene_physical_pairwise_interactions


def get_physical_interactions_BIOGRID():
    db_dir = "/Users/bjarnold/Princeton_DataX/Epistasis/higher_order_reanalysis/yeast_screens/database/BIOGRID"
    db_interactions = pd.read_csv(f"{db_dir}/BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-4.4.211.tab3.txt", sep="\t")
    db_interactions = db_interactions[['Official Symbol Interactor A', 'Official Symbol Interactor B', 'Experimental System Type']]
    db_interactions = db_interactions.rename(columns={'Official Symbol Interactor A':'official_symbol_interactor_a',
                                    'Official Symbol Interactor B':'official_symbol_interactor_b',
                                    'Experimental System Type':'experimental_system_type'})
    db_interactions = db_interactions[db_interactions.experimental_system_type == "physical"]
    db_interactions = db_interactions.reset_index(drop=True)

    gene_physical_pairwise_interactions = find_unique_interactions(db_interactions, 'official_symbol_interactor_a', 'official_symbol_interactor_b')

    return gene_physical_pairwise_interactions

    

def get_physical_interactions_yeastGenomeDotOrg():
    db_dir = "/Users/bjarnold/Princeton_DataX/Epistasis/higher_order_reanalysis/yeast_screens/database/yeastgenome_dot_org"
    db_interactions = pd.read_csv(f"{db_dir}/interaction_data.tab", sep="\t", header=None)
    db_interactions = db_interactions.rename(columns={0:"feature_name_bait",
                                    1:"standard_gene_name_bait",
                                    2:"feature_name_hit",
                                    3:"standard_gene_name_hit",
                                    4:"experiment_type",
                                    5:"genetic_or_physical_interaction",
                                    6:"source",
                                    7:"man_curated_or_high_thruput",
                                    8:"notes",
                                    9:"phenotype",
                                    10:"refrerence",
                                    11:"citation"})

    db_interactions = db_interactions.loc[(db_interactions.standard_gene_name_bait.notna()) & (db_interactions.standard_gene_name_hit.notna())]
    db_interactions = db_interactions[db_interactions.genetic_or_physical_interaction == "physical interactions"]
    db_interactions = db_interactions.reset_index(drop=True)

    gene_physical_pairwise_interactions = find_unique_interactions(db_interactions, 'standard_gene_name_bait', 'standard_gene_name_hit')
    """gene_physical_pairwise_interactions = set()
    for genes in zip(db_interactions.standard_gene_name_bait, db_interactions.standard_gene_name_hit):
        gene_physical_pairwise_interactions.add( tuple(sorted((gene_stem_name(genes[0].upper()), gene_stem_name(genes[1].upper())))) )"""

    return gene_physical_pairwise_interactions

def get_go_info():
    db_dir = "/Users/bjarnold/Princeton_DataX/Epistasis/higher_order_reanalysis/yeast_screens/database/yeastgenome_dot_org"

    df_gene_2_go = pd.read_csv(f"{db_dir}/go_slim_mapping.tab", sep="\t", header=None)
    df_gene_2_go = df_gene_2_go.rename(columns = {0:"ORF", 
                                                    1:"Gene", 
                                                    2:"SGDID", 
                                                    3:"GO_Aspect", 
                                                    4:"GO_Slim_term", 
                                                    5:"GOID", 
                                                    6:"feature_type"})
    #df_go_def = pd.read_csv(f"{db_dir}/go_terms.tab", sep="\t", header=None)

    gene_2_go = defaultdict(list)
    for i,r in df_gene_2_go.iterrows():
        go_id = r["GOID"].split(":")[1]
        gene_2_go[r['Gene']].append( go_id )

    return gene_2_go
