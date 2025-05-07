#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 23 16:10:13 2025

@author: rafael
"""

import os
import pandas as pd

# Configuração dos caminhos
files = {
    "Escherichia coli": {
        "genes": "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Genes/E_coli/group_counts_total.tsv",
        "phages": "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Phages/E.coli/E_coli_count_total_phages_sizes.tsv"
    },
    "Salmonella Enterica": {
        "genes": "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Genes/Salmonella/group_counts_total.tsv",
        "phages": "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Phages/Salmonella_enterica/Salmonella_count_total_phages_sizes.tsv"
    },
    "Klebsiella pneumoniae": {
        "genes": "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Genes/Klebsiella/group_counts_total.tsv",
        "phages": "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Phages/Klebsiella_p/Klebsiella_count_total_phages_sizes.tsv"
    }
}

output_dir = "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Diversity_Results_Rarefied/"
os.makedirs(output_dir, exist_ok=True)

def process_bacteria(name, paths):
    # Carrega e processa dados de genes
    genes_df = pd.read_csv(paths["genes"], sep="\t")
    genes_agg = genes_df.groupby(["SAMPLE", "Group"], as_index=False)["Gene_Count"].sum()
    
    # Carrega dados de fagos e obtém lista de amostras
    phages_df = pd.read_csv(paths["phages"], sep="\t")
    phage_samples = phages_df["Sample"].unique()
    
    # Filtra genes para amostras com fagos
    filtered_genes = genes_agg[genes_agg["SAMPLE"].isin(phage_samples)]
    
    # Salva o resultado
    output_path = os.path.join(output_dir, f"{name.replace(' ', '_')}_filtered_genes.tsv")
    filtered_genes.to_csv(output_path, sep="\t", index=False)
    
    return filtered_genes.head()

# Executa para cada bactéria
sample_outputs = {}
for name, paths in files.items():
    sample_outputs[name] = process_bacteria(name, paths)