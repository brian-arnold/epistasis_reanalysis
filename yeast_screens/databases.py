
import pandas as pd
import numpy as np
from collections import defaultdict
import glob, os


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
    # see here for explanation of experimental evidence codes: https://wiki.thebiogrid.org/doku.php/experimental_systems
    db_dir = "/Users/bjarnold/Princeton_DataX/Epistasis/higher_order_reanalysis/yeast_screens/database/BIOGRID"
    db_interactions = pd.read_csv(f"{db_dir}/BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-4.4.211.tab3.txt", sep="\t")
    db_interactions = db_interactions[['Official Symbol Interactor A', 'Official Symbol Interactor B', 'Experimental System Type', 'Experimental System']]
    db_interactions = db_interactions.rename(columns={'Official Symbol Interactor A':'official_symbol_interactor_a',
                                    'Official Symbol Interactor B':'official_symbol_interactor_b',
                                    'Experimental System' : 'experimental_system',
                                    'Experimental System Type':'experimental_system_type'})
    #db_interactions = db_interactions[db_interactions.experimental_system_type == "physical"]
    #db_interactions = db_interactions[db_interactions.experimental_system == "Affinity Capture-MS"]
    db_interactions = db_interactions.reset_index(drop=True)

    #gene_physical_pairwise_interactions = find_unique_interactions(db_interactions, 'official_symbol_interactor_a', 'official_symbol_interactor_b')

    return db_interactions

    

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


def count_physical_interactions(df, gene_physical_pairwise_interactions):
    num_physical_interactions = {}
    twoplus_physical_interactions = {}
    three_physical_interactions = {}

    for i,r in df.iterrows():
        # alleles in gene_physical_pairwise_interactions are sorted
        alleles = sorted(r['alleles'].split(","))
        alleles = [gene_stem_name(i.upper()) for i in alleles]
        allele_pairs = [tuple([alleles[0],alleles[1]]), 
                        tuple([alleles[0],alleles[2]]), 
                        tuple([alleles[1],alleles[2]])]
        num_physical_interactions[r['alleles']] = 0
        twoplus_physical_interactions[r['alleles']] = 0
        three_physical_interactions[r['alleles']] = 0

        for p in allele_pairs:
            if p in gene_physical_pairwise_interactions:
                num_physical_interactions[r['alleles']] += 1

        if num_physical_interactions[r['alleles']] >= 2:
            twoplus_physical_interactions[r['alleles']] = 1
        if num_physical_interactions[r['alleles']] == 3:
            three_physical_interactions[r['alleles']] = 1

    return num_physical_interactions, twoplus_physical_interactions, three_physical_interactions


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


def get_entrezID_2_geneName():

    entrezID_2_geneName_file = "/Users/bjarnold/Princeton_DataX/Epistasis/higher_order_reanalysis/yeast_screens/database/coexpressdb/entrezid_conv/Saccharomyces_cerevisiae.gene_info" 

    df = pd.read_csv(entrezID_2_geneName_file, sep="\t")
    df.Symbol = df.Symbol.str.upper()
    # take only gene name; numbers following hyphens represent alleles
    df.loc[:,'Symbol2'] = [i[0] for i in df.Symbol.str.split("-")]
    entrezID_2_geneName = dict(zip(df.GeneID, 
                                    df.Symbol2))
    
    return entrezID_2_geneName



def get_coexpression_gene_pairs(z_score_threshold):
    
    coexpress_dir = "/Users/bjarnold/Princeton_DataX/Epistasis/higher_order_reanalysis/yeast_screens/database/coexpressdb/union"
    """
    This directory contains one file per gene, each named with Entrez ID. Each file has a list of each other gene along with a normalized z score
    to measure the degree of coexpression.
    """
    entrezID_2_geneName = get_entrezID_2_geneName()

    coexpression_gene_pairs = set()
        #for genes in zip(df[left_gene], df[right_gene]):
        #    gene_physical_pairwise_interactions.add( tuple(sorted((gene_stem_name(genes[0].upper()), gene_stem_name(genes[1].upper())))) )

    for file in glob.glob(f"{coexpress_dir}/*"):
        entrez_id_gene1 = int(os.path.basename(file))
        for i in open(file, 'r').readlines():
            i_parse = i.strip().split('\t')
            entrez_id_gene2 = int(i_parse[0])
            coex_z_score = float(i_parse[1])
            if coex_z_score >= z_score_threshold:
                gene_pair = tuple(sorted((entrezID_2_geneName[entrez_id_gene1], entrezID_2_geneName[entrez_id_gene2])))
                coexpression_gene_pairs.add( gene_pair )

    return coexpression_gene_pairs




gene_2_go = get_go_info()
