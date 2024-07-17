
import pandas as pd
import numpy as np
from collections import defaultdict
import glob, os
import re

# global variable to store the location of database files
BASE_PATH = None

def set_base_path(path : str) -> None :
    global BASE_PATH
    BASE_PATH = path
    return


def gene_stem_name(gene : str) -> str:
    return gene.split('-')[0]

def find_unique_interactions(df : pd.DataFrame, left_gene : str, right_gene : str) -> set:
    """
    Assumes a dataframe has 2 columns, one for each of the physically interacting genes.
    This function returns a set with elements that are tuples of the two genes in the interaction.
    """
    physical_pairwise_interactions_set = set()
    for genes in zip(df[left_gene], df[right_gene]):
        physical_pairwise_interactions_set.add( tuple(sorted((gene_stem_name(genes[0].upper()), gene_stem_name(genes[1].upper())))) )

    return physical_pairwise_interactions_set


def get_physical_interactions_BIOGRID() -> pd.DataFrame:
    # see here for explanation of experimental evidence codes: https://wiki.thebiogrid.org/doku.php/experimental_systems
    db_dir = f"{BASE_PATH}/BIOGRID"
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


def get_physical_interactions_yeastGenomeDotOrg() -> set:
    db_dir = f"{BASE_PATH}/yeastgenome_dot_org"
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


def count_interactions_in_set(df : pd.DataFrame, interaction_set : set) -> tuple[dict, dict, dict, dict]:
    """
    takes as input the main dataframe containing fitness and interaction values for each gene triplet and a set of physical interactions
    for each triplet of genes in the 'alleles' column of the dataframe, we count the number of physical interactions and store in a dict
    """

    num_physical_interactions = {}
    oneplus_physical_interactions = {}
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
            if p in interaction_set:
                num_physical_interactions[r['alleles']] += 1

        if num_physical_interactions[r['alleles']] >= 1:
            oneplus_physical_interactions[r['alleles']] = 1
        if num_physical_interactions[r['alleles']] >= 2:
            twoplus_physical_interactions[r['alleles']] = 1
        if num_physical_interactions[r['alleles']] == 3:
            three_physical_interactions[r['alleles']] = 1

    return num_physical_interactions, oneplus_physical_interactions, twoplus_physical_interactions, three_physical_interactions

def collect_interactions_in_dict(df : pd.DataFrame, left_gene : str, right_gene : str) -> dict:
    """
    Assumes a dataframe has 2 columns, one for each of the physically interacting genes
    """
    physical_pairwise_interactions_dict = defaultdict(set)
    for genes in zip(df[left_gene], df[right_gene]):
        gene1 = gene_stem_name(genes[0].upper())
        gene2 = gene_stem_name(genes[1].upper())
        physical_pairwise_interactions_dict[gene1].add(gene2)
        physical_pairwise_interactions_dict[gene2].add(gene1)

    return physical_pairwise_interactions_dict

def count_shared_interactions_in_dict(df : pd.DataFrame, interaction_dict : dict, num_shared_interactions : int) -> dict:
    """
    calculates whether a set of three genes has at least 1 interactor in common, where the interactor is some gene 
    that isn't included in the set of three genes
    """
    three_shared_physical_interactions = {}

    for i,r in df.iterrows():
        # alleles in gene_physical_pairwise_interactions are sorted
        alleles = sorted(r['alleles'].split(","))
        alleles = [gene_stem_name(i.upper()) for i in alleles]

        three_shared_physical_interactions[r['alleles']] = 0
        # collect all interactors for this gene triplet
        interactors = set()
        for a in alleles:
            for i in interaction_dict[a]:
                if i not in alleles:
                    interactors.add(i)
        # find whether all three genes in set share at least one interactor
        cnt=0
        for i in interactors:
            if i in interaction_dict[alleles[0]] and i in interaction_dict[alleles[1]] and i in interaction_dict[alleles[2]]:
                cnt+=1
                if cnt >= num_shared_interactions:
                    three_shared_physical_interactions[r['alleles']] = 1
                    break

    return three_shared_physical_interactions


def get_go_info():
    db_dir = f"{BASE_PATH}/yeastgenome_dot_org"

    df_gene_2_go = pd.read_csv(f"{db_dir}/go_slim_mapping.tab", sep="\t", header=None)
    df_gene_2_go = df_gene_2_go.rename(columns = {0:"ORF", 
                                                    1:"Gene", 
                                                    2:"SGDID", 
                                                    3:"GO_Aspect", 
                                                    4:"GO_Slim_term", 
                                                    5:"GOID", 
                                                    6:"feature_type"})
    #df_go_def = pd.read_csv(f"{db_dir}/go_terms.tab", sep="\t", header=None)
    # Li et al 2010 ("The cellular robustness by genetic redundancy in budding yeast") only look at biological process
    df_gene_2_go = df_gene_2_go[df_gene_2_go.GO_Aspect == "P"]

    gene_2_go = defaultdict(list)
    goid_2_term = defaultdict()
    for i,r in df_gene_2_go.iterrows():
        go_id = r["GOID"].split(":")[1]
        gene_2_go[r['Gene']].append( go_id )
        goid_2_term[go_id] = r['GO_Slim_term']
    
    return gene_2_go, goid_2_term

def get_go_protein_complexes():
    db_dir = f"{BASE_PATH}/yeastgenome_dot_org"

    with open(f"{db_dir}/go_protein_complex_slim.tab", 'r') as f:
        gene_2_protein_complex = defaultdict(list)
        for line in f:
            line = line.split('\t')
            m = re.search(r'^Component: (.+)/\d+$', line[0])
            #print(line[0])
            #print(m.group(1))
            complex_term = m.group(1)

            line[1] = line[1].strip()
            line[1] = line[1].replace("Verified/","") # data are pipe separated, with some elements containing this string
            line[1] = line[1].replace("||","|") # removing above Verified/ string created double pipes
            line[1] = line[1].strip("|") # removing Verified/ above creates trailing pipe
            #print(line[1])
            for gene in line[1].split('|'):
                #print(gene)
                gene_parsed = gene.split('/')
                gene_common_name = gene_parsed[1]
                #print(gene_common_name)
                
                gene_2_protein_complex[gene_common_name].append(complex_term)

    return gene_2_protein_complex

def count_shared_go(df : pd.DataFrame) -> dict:

    gene_2_go, goid_2_term = get_go_info()

    shared_go = {}
    for i,r in df.iterrows():
        alleles = sorted(r['alleles'].split(","))
        alleles = [gene_stem_name(i.upper()) for i in alleles]

        go_counts = defaultdict(int)
        shared_go[r['alleles']] = 0
        for a in alleles:
            if a in gene_2_go:
                # many genes are involved in many GO categories; iterate through these
                for go in gene_2_go[a]:
                    go_counts[go] += 1
    
        counts = np.array([i[1] for i in go_counts.items()])
        #print(np.max(counts))
        if len(counts) > 0:
            assert np.max(counts) <= 3
        if np.sum(counts == 3) >= 1:
            shared_go[r['alleles']] = 1
    return shared_go


def get_entrezID_2_geneName():

    entrezID_2_geneName_file = f"{BASE_PATH}/coexpressdb/entrezid_conv/Saccharomyces_cerevisiae.gene_info" 

    df = pd.read_csv(entrezID_2_geneName_file, sep="\t")
    df.Symbol = df.Symbol.str.upper()
    # take only gene name; numbers following hyphens represent alleles
    df.loc[:,'Symbol2'] = [i[0] for i in df.Symbol.str.split("-")]
    entrezID_2_geneName = dict(zip(df.GeneID, 
                                    df.Symbol2))
    
    return entrezID_2_geneName


def get_expression_gene_pairs(z_score_threshold : int) -> tuple[set, set]:
    
    coexpress_dir = f"{BASE_PATH}/coexpressdb/union"
    """
    This directory contains one file per gene, each named with Entrez ID. Each file has a list of each other gene along with a normalized z score
    to measure the degree of coexpression.
    """
    entrezID_2_geneName = get_entrezID_2_geneName()

    coexpression_gene_pairs_set = set()
    divexpression_gene_pairs_set = set()
        #for genes in zip(df[left_gene], df[right_gene]):
        #    gene_physical_pairwise_interactions.add( tuple(sorted((gene_stem_name(genes[0].upper()), gene_stem_name(genes[1].upper())))) )

    for file in glob.glob(f"{coexpress_dir}/*"):
        entrez_id_gene1 = int(os.path.basename(file))
        for i in open(file, 'r').readlines():
            i_parse = i.strip().split('\t')
            entrez_id_gene2 = int(i_parse[0])
            coex_z_score = float(i_parse[1])
            gene_pair = tuple(sorted((entrezID_2_geneName[entrez_id_gene1], entrezID_2_geneName[entrez_id_gene2])))
            if coex_z_score >= abs(z_score_threshold):
                coexpression_gene_pairs_set.add( gene_pair )
            elif coex_z_score <= (-1)*abs(z_score_threshold):
                divexpression_gene_pairs_set.add( gene_pair )

    return coexpression_gene_pairs_set, divexpression_gene_pairs_set



# gene_2_protein_complex = get_go_protein_complexes()