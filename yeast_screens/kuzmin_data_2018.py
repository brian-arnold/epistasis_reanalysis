import pandas as pd
import numpy as np
from collections import defaultdict


def load_kuzmin_2018_s1(kuzmin_2018_data_dir):

    kuzmin_2018_s1 = pd.read_csv(f"{kuzmin_2018_data_dir}/Data File S1_Raw genetic interaction dataset.tsv", sep="\t")
    
    kuzmin_2018_s1 = kuzmin_2018_s1.rename(columns = {"Query strain ID" : "query_strain_id",
                    "Query allele name": "query_allele_name",
                    "Array strain ID": "array_strain_id",
                    "Array allele name": "array_allele_name",
                    "Combined mutant type": "combined_mutant_type",
                    "Raw genetic interaction score (epsilon)": "raw_interaction_score_epsilon",
                    "Adjusted genetic interaction score (epsilon or tau)": "adjusted_interaction_score_epsilon_or_tau",
                    "P-value": "pval",
                    "Query single/double mutant fitness": "query_single_double_mutant_fitness",
                    "Array single mutant fitness": "array_single_mutant_fitness",
                    "Combined mutant fitness": "double_triple_mutant_fitness",
                    "Combined mutant fitness standard deviation": "double_triple_mutant_fitness_std"})

    kuzmin_2018_s1 = kuzmin_2018_s1[["query_allele_name", 
                        "array_allele_name", 
                        "combined_mutant_type", 
                        "raw_interaction_score_epsilon", 
                        "adjusted_interaction_score_epsilon_or_tau", 
                        "pval", 
                        "query_single_double_mutant_fitness", 
                        "array_single_mutant_fitness", 
                        "double_triple_mutant_fitness"]]

    kuzmin_2018_s1.array_allele_name = kuzmin_2018_s1.array_allele_name.str.replace("Δ","") 
    
    kuzmin_2018_s1.query_allele_name = kuzmin_2018_s1.query_allele_name.str.replace("Δ","")
    kuzmin_2018_s1.query_allele_name = kuzmin_2018_s1.query_allele_name.str.replace("+",",")

    kuzmin_2018_s1[['query1_allele_name', 'query2_allele_name']] = kuzmin_2018_s1.query_allele_name.str.split(",", expand=True)
    # combine allele, DO NOT SORT, order matters here since the first two alleles are the query alleles, and for trigenic mutants, 
    # the fitness of the query double mutant has to be taken from a specific place
    kuzmin_2018_s1['alleles'] = [','.join(tup) for tup in zip(kuzmin_2018_s1.query1_allele_name, kuzmin_2018_s1.query2_allele_name, kuzmin_2018_s1.array_allele_name)]

    kuzmin_2018_s1_trigenic = kuzmin_2018_s1[kuzmin_2018_s1.combined_mutant_type == 'trigenic'].drop(columns=['combined_mutant_type'])
    kuzmin_2018_s1_trigenic = kuzmin_2018_s1_trigenic.reset_index(drop=True)

    kuzmin_2018_s1_digenic = kuzmin_2018_s1[kuzmin_2018_s1.combined_mutant_type == 'digenic'].drop(columns=['combined_mutant_type'])
    kuzmin_2018_s1_digenic.alleles = kuzmin_2018_s1_digenic.alleles.str.replace("ho,","")
    kuzmin_2018_s1_digenic.alleles = kuzmin_2018_s1_digenic.alleles.str.replace(",ho","")
    kuzmin_2018_s1_digenic = kuzmin_2018_s1_digenic.reset_index(drop=True)

    #dict(zip(kuzmin_2018_s1_digenic.array_allele_name, kuzmin_2018_s1_digenic.array_single_mutant_fitness))

    kuzmin_2018_s1_digenic.query_allele_name = kuzmin_2018_s1_digenic.query_allele_name.str.replace("\+ho","")
    kuzmin_2018_s1_digenic.query_allele_name = kuzmin_2018_s1_digenic.query_allele_name.str.replace("ho\+","")

    return kuzmin_2018_s1, kuzmin_2018_s1_digenic, kuzmin_2018_s1_trigenic 


def load_costanzo_data(costanzo_et_al_data_dir):

    costanzo_nxn = pd.read_csv(f"{costanzo_et_al_data_dir}/SGA_NxN.txt", sep="\t")
    costanzo_exn = pd.read_csv(f"{costanzo_et_al_data_dir}/SGA_ExN.txt", sep="\t")
    costanzo_exe = pd.read_csv(f"{costanzo_et_al_data_dir}/SGA_ExE.txt", sep="\t")

    costanzo_data_column_rename = {"Query Strain ID" : "query_strain_id",
                            "Query allele name": "query_allele_name",
                            "Array Strain ID": "array_strain_id",
                            "Array allele name": "array_allele_name",
                            "Arraytype/Temp" : "Arraytype_Temp",
                            "Genetic interaction score (ε)" : "genetic_interaction_e",
                            "P-value": "pval",
                            "Query single mutant fitness (SMF)": "query_smf",
                            "Array SMF": "array_smf",
                            "Double mutant fitness": "dbl_mutant_fitness",
                            "Double mutant fitness standard deviation": "dbl_mutant_fitness_std"}

    costanzo_nxn = costanzo_nxn.rename(columns = costanzo_data_column_rename)
    costanzo_exn = costanzo_exn.rename(columns = costanzo_data_column_rename)
    costanzo_exe = costanzo_exe.rename(columns = costanzo_data_column_rename)
    costanzo = pd.concat([costanzo_nxn, costanzo_exn, costanzo_exe])

    costanzo = costanzo.astype({'query_strain_id': 'str',
                        'query_allele_name': 'str',
                        'array_strain_id': 'str',
                        'array_allele_name': 'str',
                        'Arraytype_Temp': 'str'})

    costanzo = costanzo[['query_allele_name', 'array_allele_name', 'query_smf', 'array_smf', 'dbl_mutant_fitness', 'genetic_interaction_e']]

    costanzo['alleles'] = [','.join(sorted(tup)) for tup in zip(costanzo.query_allele_name, costanzo.array_allele_name)]

    costanzo = costanzo.reset_index(drop=True)
    
    return costanzo


def consolidate_fitnesses_across_2018_tables(kuzmin_2018_s1_trigenic,
                                            f_k_SMF,
                                            f_i_j_SMF,
                                            f_ij_DMF,
                                            f_ik_jk_DMF,
                                            e_ik_jk_DMF):

    f_ij, f_ik, f_jk = defaultdict(), defaultdict(), defaultdict()
    f_i, f_j, f_k = defaultdict(), defaultdict(), defaultdict()
    e_ik_kuz, e_jk_kuz = defaultdict(), defaultdict()

    for i, row in kuzmin_2018_s1_trigenic.iterrows():
        alleles = row['alleles'].split(",")
        a = row['alleles']

        # ij is query double mutant
        ij = ",".join([alleles[0],alleles[1]])
        ik = ",".join([alleles[0],alleles[2]])
        jk = ",".join([alleles[1],alleles[2]])

        ji = ",".join(([alleles[1],alleles[0]]))
        ki = ",".join(([alleles[2],alleles[0]]))
        kj = ",".join(([alleles[2],alleles[1]]))


        i,j,k = alleles[0], alleles[1], alleles[2]
        
        try:
            f_ij[a] = f_ij_DMF[ij]
        except KeyError:
            try:
                f_ij[a] = f_ij_DMF[ji]
            except KeyError:
                try:
                     f_ij[a] = f_ik_jk_DMF[ij]
                except KeyError:
                    try:
                        f_ij[a] = f_ik_jk_DMF[ji]
                    except KeyError:
                        f_ij[a] = float('nan')

        try:
            f_ik[a] = f_ik_jk_DMF[ik]
        except KeyError:
            try:
                f_ik[a] = f_ik_jk_DMF[ki]
            except KeyError:
                f_ik[a] = float('nan')

        try:
            f_jk[a] = f_ik_jk_DMF[jk]
        except KeyError:
            try:
                f_jk[a] = f_ik_jk_DMF[kj]
            except KeyError:
                f_jk[a] = float('nan')
        
        
        # get single mutant fitnesses
        # Supplement specifically states that f_i and f_j are taken from a previous study
        try:
            f_i[a] = f_i_j_SMF[i]
        except KeyError:
            try:
                f_i[a] = f_k_SMF[i]
            except KeyError:
                f_i[a] = float('nan')
        
        try:
            f_j[a] = f_i_j_SMF[j]
        except KeyError:
            try:
                f_j[a] = f_k_SMF[j]
            except KeyError:
                f_j[a] = float('nan')
        
        # take f_k from current study Kuzmin 2018
        try:
            f_k[a] = f_k_SMF[k]
        except KeyError:
            f_k[a] = float('nan')
                
        # get pairwise epistasis
        try:
            e_ik_kuz[a] = e_ik_jk_DMF[ik]
        except KeyError:
            try:
                e_ik_kuz[a] = e_ik_jk_DMF[ki]
            except KeyError:
                e_ik_kuz[a] = float('nan')

        try:
            e_jk_kuz[a] = e_ik_jk_DMF[jk]
        except KeyError:
            try:
                e_jk_kuz[a] = e_ik_jk_DMF[kj]
            except KeyError:
                e_jk_kuz[a] = float('nan')

    return f_i, f_j, f_k, f_ij, f_ik, f_jk, e_ik_kuz, e_jk_kuz