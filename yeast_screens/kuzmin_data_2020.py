import pandas as pd
import numpy as np
from collections import defaultdict



def load_kuzmin_2020_s1(data_dir : str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    # Load and process Table S1 from Kuzmin et al. 2020
    # More information on the information stored in each column can be found in the supplementary text of Kuzmin et al 2020

    kuzmin_2020_s1 = pd.read_csv(f"{data_dir}/aaz5667-Table-S1.txt", sep="\t")
    
    kuzmin_2020_s1 = kuzmin_2020_s1.rename(columns = {"Query strain ID" : "query_strain_id",
                                                        "Query allele name": "query_allele_name",
                                                        "Array strain ID": "array_strain_id",
                                                        "Array allele name": "array_allele_name",
                                                        "Combined mutant type": "combined_mutant_type",
                                                        "Raw genetic interaction score (epsilon)": "raw_interaction_score_epsilon",
                                                        "Adjusted genetic interaction score (epsilon or tau)": "adjusted_interaction_score_epsilon_or_tau",
                                                        "P-value": "pval",
                                                        "Query single/double mutant fitness": "query_single_double_mutant_fitness",
                                                        "Array single mutant fitness": "array_single_mutant_fitness",
                                                        "Double/triple mutant fitness": "double_triple_mutant_fitness",
                                                        "Double/triple mutant fitness standard deviation": "double_triple_mutant_fitness_std"})

    kuzmin_2020_s1.query_allele_name = kuzmin_2020_s1.query_allele_name.str.replace("_","")
    kuzmin_2020_s1.array_allele_name = kuzmin_2020_s1.array_allele_name.str.replace("_","")

    kuzmin_2020_s1[['query1_allele_name', 'query2_allele_name']] = kuzmin_2020_s1.query_allele_name.str.split("+", expand=True)
    kuzmin_2020_s1['alleles'] = kuzmin_2020_s1[['query1_allele_name', 'query2_allele_name', 'array_allele_name']].agg(','.join, axis=1)

    kuzmin_2020_s1_trigenic = kuzmin_2020_s1[kuzmin_2020_s1.combined_mutant_type == 'trigenic'].drop(columns=['combined_mutant_type'])
    kuzmin_2020_s1_trigenic = kuzmin_2020_s1_trigenic.reset_index(drop=True)
    kuzmin_2020_s1_digenic = kuzmin_2020_s1[kuzmin_2020_s1.combined_mutant_type == 'digenic'].drop(columns=['combined_mutant_type'])
    # ho will either be the 1st or second query allele, never the array allele
    kuzmin_2020_s1_digenic.alleles = kuzmin_2020_s1_digenic.alleles.str.replace("ho,","")
    kuzmin_2020_s1_digenic = kuzmin_2020_s1_digenic.reset_index(drop=True)

    return kuzmin_2020_s1, kuzmin_2020_s1_digenic, kuzmin_2020_s1_trigenic


def load_kuzmin_2020_s2(data_dir : str) -> pd.DataFrame:

    kuzmin_2020_s2 = pd.read_csv(f"{data_dir}/aaz5667-Table-S2.txt", sep="\t")

    kuzmin_2020_s2 = kuzmin_2020_s2.rename(columns = {"Query strain ID" : "query_strain_id",
                        "Query allele name": "query_allele_name",
                        "Array strain ID": "array_strain_id",
                        "Array allele name": "array_allele_name",
                        "Combined mutant type": "combined_mutant_type",
                        "Adjusted genetic interaction score (epsilon or tau)": "adjusted_interaction_score_epsilon_or_tau",
                        "P-value": "pval",
                        "Digenic, Modified trigenic, or Novel trigenic": "digenic_modifiedTrigenic_novelTrigenic"})

    kuzmin_2020_s2[['query1_allele_name', 'query2_allele_name']] = kuzmin_2020_s2.query_allele_name.str.split("+", expand=True)
    kuzmin_2020_s2['three_alleles'] = kuzmin_2020_s2[['query1_allele_name', 'query2_allele_name', 'array_allele_name']].agg(','.join, axis=1)
    kuzmin_2020_s2.three_alleles = kuzmin_2020_s2.three_alleles.str.replace("_","")

    return kuzmin_2020_s2



def load_kuzmin_2020_s3(data_dir : str) -> tuple[pd.DataFrame, pd.DataFrame]:
    kuzmin_2020_s3 = pd.read_csv(f"{data_dir}/aaz5667-Table-S3.txt", sep="\t")
    # note on s3: query allele fitness for digenic crosses are all NaN, only array allele fitness available
    kuzmin_2020_s3 = kuzmin_2020_s3.rename(columns = {"Query strain ID" : "query_strain_id",
                            "Query allele name": "query_allele_name",
                            "Array strain ID": "array_strain_id",
                            "Array allele name": "array_allele_name",
                            "Combined mutant type": "combined_mutant_type",
                            "Raw genetic interaction score (epsilon)": "raw_interaction_score_epsilon",
                            "Adjusted genetic interaction score (epsilon or tau)": "adjusted_interaction_score_epsilon_or_tau",
                            "P-value": "pval",
                            "Query single/double mutant fitness": "query_single_double_mutant_fitness",
                            "Array single mutant fitness": "array_single_mutant_fitness",
                            "Double/triple mutant fitness": "double_triple_mutant_fitness",
                            "Double/triple mutant fitness standard deviation": "double_triple_mutant_fitness_std"})

    kuzmin_2020_s3.query_allele_name = kuzmin_2020_s3.query_allele_name.str.replace("_","")
    kuzmin_2020_s3.array_allele_name = kuzmin_2020_s3.array_allele_name.str.replace("_","")

    kuzmin_2020_s3_digenic = kuzmin_2020_s3[kuzmin_2020_s3.combined_mutant_type == "digenic"]

    return kuzmin_2020_s3, kuzmin_2020_s3_digenic


def load_kuzmin_2020_s5(data_dir : str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    kuzmin_2020_s5 = pd.read_csv(f"{data_dir}/aaz5667-Table-S5.txt", sep="\t")
    # S5 table with single and double mutant fitnesses
    kuzmin_2020_s5 = kuzmin_2020_s5.rename(columns = {"Allele1" : "allele1",
                                                "Allele2" : "allele2",
                                                "Mutant type" : "mutant_type",
                                                "Fitness" : "fitness",
                                                "St.dev." : "std"})

    kuzmin_2020_s5.allele1 = kuzmin_2020_s5.allele1.str.replace("_","")
    kuzmin_2020_s5.allele2 = kuzmin_2020_s5.allele2.str.replace("_","")
    kuzmin_2020_s5['alleles'] = kuzmin_2020_s5[['allele1', 'allele2']].agg(','.join, axis=1)
    kuzmin_2020_s5.alleles = kuzmin_2020_s5.alleles.str.replace("ho,","")
    kuzmin_2020_s5.alleles = kuzmin_2020_s5.alleles.str.replace("ho","")
    kuzmin_2020_s5.alleles = kuzmin_2020_s5.alleles.str.strip(",")

    kuzmin_2020_s5_dblMut = kuzmin_2020_s5[kuzmin_2020_s5.mutant_type == "Double mutant"]
    kuzmin_2020_s5_singMut = kuzmin_2020_s5[kuzmin_2020_s5.mutant_type == "Single mutant"]

    return kuzmin_2020_s5, kuzmin_2020_s5_singMut, kuzmin_2020_s5_dblMut


def load_kuzmin_2020_s10(data_dir : str) -> pd.DataFrame:
    kuzmin_2020_s10 = pd.read_csv(f"{data_dir}/aaz5667-Table-S10.txt", sep="\t")
    kuzmin_2020_s10 = kuzmin_2020_s10.rename(columns={"Seq div rate": "seq_div_rate"})
    return kuzmin_2020_s10



def consolidate_fitnesses_across_tables(kuzmin_2020_s1_trigenic, kuzmin_2020_s1_dblMutFit, kuzmin_2020_s1_singMutFit, kuzmin_2020_s1_epsilon, kuzmin_2020_s1_epsilon_pvals, kuzmin_2020_s3_singMutFit, kuzmin_2020_s5_singMutFit, kuzmin_2020_s5_dblMutFit):

    f_ij, f_ik, f_jk = defaultdict(), defaultdict(), defaultdict()
    f_i, f_j, f_k = defaultdict(), defaultdict(), defaultdict()
    e_ik_kuz, e_ik_kuz_pval = defaultdict(), defaultdict()
    e_jk_kuz, e_jk_kuz_pval = defaultdict(), defaultdict()

    for i, row in kuzmin_2020_s1_trigenic.iterrows():
    #if i <= 100:
        alleles = row['alleles'].split(",")
        ij = alleles[0] + "," + alleles[1]
        ji = alleles[1] + "," + alleles[0]

        ik = alleles[0] + "," + alleles[2]
        ki = alleles[2] + "," + alleles[0]

        jk = alleles[1] + "," + alleles[2]
        kj = alleles[2] + "," + alleles[1]

        i,j,k = alleles[0], alleles[1], alleles[2]
       
        try:
            f_ij[row['alleles']] = kuzmin_2020_s5_dblMutFit[ij]
        except KeyError:
            try:
                f_ij[row['alleles']] = kuzmin_2020_s5_dblMutFit[ji]
            except:
                f_ij[row['alleles']] = float('nan')
        try:
            f_ik[row['alleles']] = kuzmin_2020_s1_dblMutFit[ik]
        except KeyError:
            try:
                f_ik[row['alleles']] = kuzmin_2020_s1_dblMutFit[ki]
            except:
                f_ik[row['alleles']] = float('nan')
        try:
            f_jk[row['alleles']] = kuzmin_2020_s1_dblMutFit[jk]
        except KeyError:
            try:
                f_jk[row['alleles']] = kuzmin_2020_s1_dblMutFit[kj]
            except:
                f_jk[row['alleles']] = float('nan')
        
        
        # get single mutant fitnesses
        try:
            f_i[row['alleles']] = kuzmin_2020_s3_singMutFit[i]
        except KeyError:
            f_i[row['alleles']] = float('nan')
        try:
            f_j[row['alleles']] = kuzmin_2020_s5_singMutFit[j]
        except KeyError:
            f_j[row['alleles']] = float('nan')  
        try:
            f_k[row['alleles']] = kuzmin_2020_s1_singMutFit[k]
        except KeyError:
            f_k[row['alleles']] = float('nan')      

       
       # get pairwise epistasis
        try:
            e_ik_kuz[row['alleles']] = kuzmin_2020_s1_epsilon[ik]
        except KeyError:
            try:
                e_ik_kuz[row['alleles']] = kuzmin_2020_s1_epsilon[ki]
            except:
                e_ik_kuz[row['alleles']] = float('nan')

        try:
            e_jk_kuz[row['alleles']] = kuzmin_2020_s1_epsilon[jk]
        except KeyError:
            try:
                e_jk_kuz[row['alleles']] = kuzmin_2020_s1_epsilon[kj]
            except:
                e_jk_kuz[row['alleles']] = float('nan')

       # get pairwise epistasis pvalues
        try:
            e_ik_kuz_pval[row['alleles']] = kuzmin_2020_s1_epsilon_pvals[ik]
        except KeyError:
            try:
                e_ik_kuz_pval[row['alleles']] = kuzmin_2020_s1_epsilon_pvals[ki]
            except:
                e_ik_kuz_pval[row['alleles']] = float('nan')

        try:
            e_jk_kuz_pval[row['alleles']] = kuzmin_2020_s1_epsilon_pvals[jk]
        except KeyError:
            try:
                e_jk_kuz_pval[row['alleles']] = kuzmin_2020_s1_epsilon_pvals[kj]
            except:
                e_jk_kuz_pval[row['alleles']] = float('nan')                

    return f_i, f_j, f_k, f_ij, f_ik, f_jk, e_ik_kuz, e_jk_kuz, e_ik_kuz_pval, e_jk_kuz_pval