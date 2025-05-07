#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt

# ===========================================================
# FUNÇÃO PARA CALCULAR ÍNDICES DE DIVERSIDADE (SHANNON & SIMPSON)
# ===========================================================
def calculate_diversity_indices(proportions):
    proportions = proportions / proportions.sum()
    proportions = proportions[proportions > 0]

    shannon = -np.sum(proportions * np.log(proportions))
    simpson = 1 - np.sum(proportions ** 2)

    return {"Shannon": shannon, "Simpson": simpson}

# ===========================================================
# FUNÇÃO PARA RAREFAÇÃO DA DIVERSIDADE
# ===========================================================
def rarefaction_diversity(df, sample_sizes=[101], n_iter=33):
    results = []

    for sample_size in sample_sizes:
        shannon_values = []
        simpson_values = []

        for _ in range(n_iter):
            sampled_samples = np.random.choice(df.index.unique(), size=min(sample_size, len(df.index)), replace=False)
            sampled_data = df.loc[sampled_samples]
            proportions = sampled_data.sum() / sampled_data.sum().sum()
            proportions = proportions[proportions > 0]

            diversity = calculate_diversity_indices(proportions)
            shannon_values.append(diversity["Shannon"])
            simpson_values.append(diversity["Simpson"])

        results.append({
            "Sample_Size": sample_size,
            "Shannon": np.mean(shannon_values),
            "Simpson": np.mean(simpson_values)
        })

    return pd.DataFrame(results)

# ===========================================================
# CONFIGURAÇÃO DOS ARQUIVOS BRUTOS
# ===========================================================
files = {
    "Escherichia coli": {
        "genes": "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Diversity_Results_Rarefied/Escherichia_coli_filtered_genes.tsv",
        "phages": "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Phages/E.coli/E_coli_count_total_phages_sizes.tsv"
    },
    "Salmonella Enterica": {
        "genes": "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Diversity_Results_Rarefied/Salmonella_Enterica_filtered_genes.tsv",
        "phages": "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Phages/Salmonella_enterica/Salmonella_count_total_phages_sizes.tsv"
    },
    "Klebsiella pneumoniae": {
        "genes": "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Diversity_Results_Rarefied/Klebsiella_pneumoniae_filtered_genes.tsv",
        "phages": "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Phages/Klebsiella_p/Klebsiella_count_total_phages_sizes.tsv"
    }
}

output_dir = "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Diversity_Results_Rarefied/"
os.makedirs(output_dir, exist_ok=True)

all_genes_results = []
all_phages_results = []

for bacterium, paths in files.items():
    print(f"📊 Processando {bacterium}...")

    genes_df = pd.read_csv(paths["genes"], sep="\t")
    genes_df = genes_df.groupby(["SAMPLE", "Group"], as_index=False).agg({"Gene_Count": "sum"})
    genes_pivot = genes_df.pivot(index="SAMPLE", columns="Group", values="Gene_Count").fillna(0)

    genes_rarefied = rarefaction_diversity(genes_pivot, sample_sizes=[101], n_iter=333333)
    genes_rarefied["Bacterium"] = bacterium
    all_genes_results.append(genes_rarefied)

    phages_df = pd.read_csv(paths["phages"], sep="\t")

    if "Taxon_Level" in phages_df.columns and "Taxon_Name" in phages_df.columns:
        phages_df = phages_df[phages_df["Taxon_Level"] == "Genus"]
        phages_df = phages_df.rename(columns={"Taxon_Name": "Genus"})
    elif "Genus" not in phages_df.columns:
        raise ValueError(f"O arquivo {paths['phages']} não possui coluna 'Genus'.")

    phages_df = phages_df.groupby(["Sample", "Genus"], as_index=False).agg({"Count": "sum"})
    phages_pivot = phages_df.pivot(index="Sample", columns="Genus", values="Count").fillna(0)

    phages_rarefied = rarefaction_diversity(phages_pivot, sample_sizes=[101], n_iter=333333)
    phages_rarefied["Bacterium"] = bacterium
    all_phages_results.append(phages_rarefied)

bacterium_prefix = {
    "Escherichia coli": "e_coli",
    "Salmonella Enterica": "salmonella",
    "Klebsiella pneumoniae": "klebsiella"
}

for bacterium in files.keys():
    prefix = bacterium_prefix[bacterium]
    genes_rarefied = [x for x in all_genes_results if x["Bacterium"].iloc[0] == bacterium][0]
    phages_rarefied = [x for x in all_phages_results if x["Bacterium"].iloc[0] == bacterium][0]

    shannon_genes = genes_rarefied['Shannon'].iloc[0]
    shannon_phages = phages_rarefied['Shannon'].iloc[0]
    simpson_genes = genes_rarefied['Simpson'].iloc[0]
    simpson_phages = phages_rarefied['Simpson'].iloc[0]

    # === Gráfico Shannon ===
    plt.figure(figsize=(10, 6))
    sns.barplot(x=['Genes', 'Phages'], y=[shannon_genes, shannon_phages], hue=['Genes', 'Phages'], palette='mako', legend=False)
    plt.title(f'Index Shannon - {bacterium}', fontsize=14, style='italic')
    plt.ylabel('Value', fontsize=12)
    plt.xticks(fontsize=12)
    plt.savefig(os.path.join(output_dir, f'{prefix}_shannon.png'), dpi=700, bbox_inches='tight')
    plt.close()

    # === Gráfico Simpson ===
    plt.figure(figsize=(10, 6))
    sns.barplot(x=['Genes', 'Phages'], y=[simpson_genes, simpson_phages], hue=['Genes', 'Phages'], palette='mako', legend=False)
    plt.title(f'Index Simpson - {bacterium}', fontsize=14, style='italic')
    plt.ylabel('Value', fontsize=12)
    plt.xticks(fontsize=12)
    plt.savefig(os.path.join(output_dir, f'{prefix}_simpson.png'), dpi=700, bbox_inches='tight')
    plt.close()

# === Exportar resultados ===
genes_all_df = pd.concat(all_genes_results, ignore_index=True)
genes_all_df["Source"] = "Genes"

phages_all_df = pd.concat(all_phages_results, ignore_index=True)
phages_all_df["Source"] = "Phages"

final_df = pd.concat([genes_all_df, phages_all_df], ignore_index=True)
final_df = final_df[["Bacterium", "Source", "Sample_Size", "Shannon", "Simpson"]]

tsv_output_path = os.path.join(output_dir, "diversity_indices_summary.tsv")
final_df.to_csv(tsv_output_path, sep="\t", index=False)

print(f"📁 Resultados salvos em: {tsv_output_path}")
