#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 18 09:20:21 2025

@author: rafael

Analisa diversidade beta (Bray-Curtis) em genes e fagos com:

✅ Normalização por CSS (Cumulative Sum Scaling)
✅ Agregação por grupo bacteriano para comparação justa
✅ Teste de Kruskal-Wallis para identificar diferenças significativas por gene/fago entre os grupos
✅ Teste post-hoc de Dunn para comparações entre pares de bactérias (com filtro de p < 0.05)
✅ Geração de TSV com valores totais normalizados por bactéria por feature
✅ Geração de TSV com comparações par a par (post-hoc)
✅ Cálculo de distância Bray-Curtis entre medianas de grupos (genes e fagos)
✅ Geração de heatmaps representando a beta diversidade entre os grupos
✅ Geração de gráficos PCoA (2D e 3D) com exibição de stress
"""

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import kruskal
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import MDS
import scikit_posthocs as sp
from mpl_toolkits.mplot3d import Axes3D
from skbio.stats.distance import permanova, DistanceMatrix
# ─────────────────────────────── CONFIGURAÇÃO ───────────────────────────────

bacteria_names = {
    "E_coli": r"$\it{Escherichia\ coli}$",
    "Salmonella": r"$\it{Salmonella\ enterica}$",
    "Klebsiella": r"$\it{Klebsiella\ pneumoniae}$"
}

# ─────────────────────────────── CONFIGURAÇÃO ───────────────────────────────

input_files = {
    bacteria_names["E_coli"]: {
        "genes": "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Diversity_Results_Rarefied/Escherichia_coli_filtered_genes.tsv",
        "phages": "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Phages/E.coli/E_coli_count_total_phages_sizes.tsv"
    },
    bacteria_names["Salmonella"]: {
        "genes": "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Diversity_Results_Rarefied/Salmonella_Enterica_filtered_genes.tsv",
        "phages": "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Phages/Salmonella_enterica/Salmonella_count_total_phages_sizes.tsv"
    },
    bacteria_names["Klebsiella"]: {
        "genes": "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Diversity_Results_Rarefied/Klebsiella_pneumoniae_filtered_genes.tsv",
        "phages": "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Phages/Klebsiella_p/Klebsiella_count_total_phages_sizes.tsv"
    }
}

output_dir = "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Phages_Beta/results"
os.makedirs(output_dir, exist_ok=True)

output_dir = "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Phages_Beta/results"
os.makedirs(output_dir, exist_ok=True)


def normalize_css(df):
    df = df.loc[:, (df > 0).any(axis=0)]
    df = df.div(df.sum(axis=1), axis=0)
    return df.fillna(0)

def carregar_dados_somados():
    dados_genes = []
    dados_phages = []

    for bact, paths in input_files.items():
        genes_df = pd.read_csv(paths["genes"], sep="\t")
        g_pivot = genes_df.pivot_table(index="Sample", columns="Group", values="Gene_Count", aggfunc="sum", fill_value=0)
        g_norm = normalize_css(g_pivot)
        g_norm["__bacteria__"] = bact
        dados_genes.append(g_norm)

        phage_df = pd.read_csv(paths["phages"], sep="\t")
        if "Taxon_Level" in phage_df.columns:
            phage_df = phage_df[phage_df["Taxon_Level"] == "Genus"].rename(columns={"Taxon_Name": "Genus"})
        p_pivot = phage_df.pivot_table(index="Sample", columns="Genus", values="Count", aggfunc="sum", fill_value=0)
        p_norm = normalize_css(p_pivot)
        p_norm["__bacteria__"] = bact
        dados_phages.append(p_norm)

    df_genes = pd.concat(dados_genes).fillna(0)
    df_phages = pd.concat(dados_phages).fillna(0)
    return df_genes, df_phages


def kruskal_com_posthoc(df, tipo):
    grupos = df["__bacteria__"]
    df = df.drop(columns="__bacteria__")
    kruskal_resultados = []
    todas_comparacoes = []

    for feat in df.columns:
        try:
            valores = [df[grupos == g][feat].values for g in grupos.unique()]
            stat, pval = kruskal(*valores)
            if pval < 0.05:
                kruskal_resultados.append((feat, stat, pval))
                dunn = sp.posthoc_dunn(df[[feat]].assign(grupo=grupos), val_col=feat, group_col="grupo", p_adjust='bonferroni')
                dunn_filtered = dunn.where(dunn < 0.05)
                for g1 in dunn_filtered.index:
                    for g2 in dunn_filtered.columns:
                        if g1 != g2 and pd.notnull(dunn_filtered.loc[g1, g2]):
                            todas_comparacoes.append((feat, g1, g2, dunn_filtered.loc[g1, g2]))
        except Exception as e:
            print(f"Erro em {feat}: {e}")

    pd.DataFrame(kruskal_resultados, columns=["Feature", "Statistic", "P_value"]).to_csv(
        os.path.join(output_dir, f"{tipo}_kruskal_significant.tsv"), sep="\t", index=False)

    pd.DataFrame(todas_comparacoes, columns=["Feature", "Group1", "Group2", "P_value"]).to_csv(
        os.path.join(output_dir, f"{tipo}_dunn_posthoc.tsv"), sep="\t", index=False)


def matriz_braycurtis(df, tipo, dpi=600):
    grupos = df.pop("__bacteria__")
    medianas = df.groupby(grupos).median()
    matriz = pd.DataFrame(
        squareform(pdist(medianas, metric="braycurtis")),
        index=medianas.index,
        columns=medianas.index
    )
    matriz.to_csv(os.path.join(output_dir, f"{tipo}_braycurtis_matriz.tsv"), sep="\t")

    plt.figure(figsize=(6,5))
    sns.heatmap(matriz, annot=True, cmap="magma", square=True)
    plt.title(f"Beta diversidade ({tipo}) - Bray-Curtis")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{tipo}_braycurtis_heatmap.png"), dpi=dpi)
    plt.close()

    mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
    coords = mds.fit_transform(matriz.values)
    stress = mds.stress_
    pd.DataFrame([[tipo, stress]], columns=["Tipo", "Stress"]).to_csv(os.path.join(output_dir, f"{tipo}_mds_stress.tsv"), sep="\t", index=False)

    mds_df = pd.DataFrame(coords, columns=["Dim1", "Dim2"], index=matriz.index)
    plt.figure(figsize=(6,5))
    sns.scatterplot(x="Dim1", y="Dim2", hue=mds_df.index, data=mds_df, s=100)
    for i in mds_df.index:
        plt.text(mds_df.loc[i, "Dim1"] + 0.01, mds_df.loc[i, "Dim2"], i)
    plt.title(f"PCoA - {tipo} (Bray-Curtis)\nStress = {stress:.4f}")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{tipo}_pcoa.png"), dpi=500)
    plt.close()

    mds3d = MDS(n_components=3, dissimilarity='precomputed', random_state=42)
    coords3d = mds3d.fit_transform(matriz.values)
    fig = plt.figure(figsize=(7,6))
    ax = fig.add_subplot(111, projection='3d')
    for i, label in enumerate(matriz.index):
        ax.scatter(coords3d[i, 0], coords3d[i, 1], coords3d[i, 2], label=label, s=100)
        ax.text(coords3d[i, 0], coords3d[i, 1], coords3d[i, 2], label)
    ax.set_title(f"PCoA 3D - {tipo}")
    ax.set_xlabel("Dim1")
    ax.set_ylabel("Dim2")
    ax.set_zlabel("Dim3")
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{tipo}_pcoa_3d.png"), dpi=600)
    plt.close()

def permanova_e_dispersion(df, tipo):
    grupos = df["__bacteria__"]
    
    # Criar índice único por amostra + grupo
    df["__sample_id__"] = [f"{g}_{i}" for g, i in zip(grupos, df.index)]
    df = df.set_index("__sample_id__")
    df_metadata = pd.DataFrame({"grupo": grupos.values}, index=df.index)
    df = df.drop(columns="__bacteria__")

    # Matriz de distâncias
    dist_matrix = pd.DataFrame(
        squareform(pdist(df.values, metric="braycurtis")),
        index=df.index,
        columns=df.index
    )
    dm = DistanceMatrix(dist_matrix.values, ids=dist_matrix.index)

    # PERMANOVA
    result_perm = permanova(dm, df_metadata, column='grupo', permutations=999)

    result_perm_dict = {
        "Source": ["grupo"],
        "Test_statistic": [result_perm['test statistic']],
        "p-value": [result_perm['p-value']]
    }

    # Adiciona 'Permutations' se existir
    if 'permutations' in result_perm:
        result_perm_dict["Permutations"] = [result_perm['permutations']]

    # Cria DataFrame e salva
    result_perm_df = pd.DataFrame(result_perm_dict)
    result_perm_df.to_csv(os.path.join(output_dir, f"{tipo}_permanova.tsv"), sep="\t", index=False)

    # Beta-dispersion
    grupo_labels = df_metadata["grupo"].unique()
    dispersoes = []
    for grupo in grupo_labels:
        grupo_idx = df_metadata[df_metadata["grupo"] == grupo].index
        sub_df = df.loc[grupo_idx]
        mediana = sub_df.median()
        dists = sub_df.apply(lambda row: pdist([row.values, mediana.values], metric="braycurtis")[0], axis=1)
        for amostra, dist in dists.items():
            dispersoes.append((grupo, amostra, dist))

    dispersoes_df = pd.DataFrame(dispersoes, columns=["Grupo", "Amostra", "Distância"])
    dispersoes_df.to_csv(os.path.join(output_dir, f"{tipo}_beta_dispersion.tsv"), sep="\t", index=False)

    plt.figure(figsize=(7,5))
    sns.boxplot(data=dispersoes_df, x="Grupo", y="Distância")
    plt.title(f"Dispersão Beta - {tipo}")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{tipo}_beta_dispersion_boxplot.png"), dpi=300)
    plt.close()



# ─────────────────────────────── MAIN ───────────────────────────────

df_genes, df_phages = carregar_dados_somados()
print("✔️ Dados normalizados por amostra e agrupados por bactéria carregados!")

kruskal_com_posthoc(df_genes.copy(), "genes")
kruskal_com_posthoc(df_phages.copy(), "phages")

matriz_braycurtis(df_genes.copy(), "genes", dpi=300)
matriz_braycurtis(df_phages.copy(), "phages", dpi=300)

permanova_e_dispersion(df_genes.copy(), "genes")
permanova_e_dispersion(df_phages.copy(), "phages")

print("✅ Análise completa: Kruskal, Bray-Curtis, PERMANOVA e dispersão beta finalizados!")