#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 24 08:23:44 2025
@author: rafael

Analisa dispersão beta com:
✅ Levene (homogeneidade de variância)
✅ ANOVA (diferença de dispersão)
✅ PERMDISP-like (teste de permutação)
✅ Boxplots salvos
"""

import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols
from scipy.stats import levene
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Caminhos para os arquivos
df_genes = pd.read_csv("/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Phages_Beta/results/genes_beta_dispersion.tsv", sep="\t")
df_phages = pd.read_csv("/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Phages_Beta/results/phages_beta_dispersion.tsv", sep="\t")

def permutacao_permdisp(df, n_perm=999999, seed=42):
    np.random.seed(seed)
    grupo_original = df["Grupo"].values
    distancias = df["Distância"].values
    grupos_unicos = np.unique(grupo_original)

    def ss_between(grupos, distancias):
        medias = {g: distancias[grupos == g].mean() for g in grupos}
        grand_mean = distancias.mean()
        ssb = sum(
            [(distancias[grupos == g].shape[0]) * ((media - grand_mean) ** 2)
             for g, media in medias.items()]
        )
        return ssb

    ssb_observado = ss_between(grupo_original, distancias)

    permutacoes = []
    for _ in range(n_perm):
        grupo_perm = np.random.permutation(grupo_original)
        ssb_perm = ss_between(grupo_perm, distancias)
        permutacoes.append(ssb_perm)

    permutacoes = np.array(permutacoes)
    p_valor = (np.sum(permutacoes >= ssb_observado) + 1) / (n_perm + 1)
    return p_valor

def testar_dispersao(df, tipo):
    print(f"\n🔬 Teste de dispersão beta - {tipo.upper()}")

    # Levene
    grupos = df["Grupo"].unique()
    grupos_vals = [df[df["Grupo"] == g]["Distância"].values for g in grupos]
    stat_levene, p_levene = levene(*grupos_vals)
    print(f"Levene: stat={stat_levene:.4f}, p={p_levene:.4f}")

    # ANOVA
    modelo = ols("Distância ~ Grupo", data=df).fit()
    anova_tabela = sm.stats.anova_lm(modelo, typ=2)
    print(anova_tabela)

    # PERMDISP-like
    p_perm = permutacao_permdisp(df)
    print(f"🌀 PERMDISP-like (999 permutações) → p = {p_perm:.4f}")

    # Boxplot
    plt.figure(figsize=(7, 5))
    sns.boxplot(data=df, x="Grupo", y="Distância")
    plt.title(f"Dispersão Beta - {tipo}")
    plt.tight_layout()
    plt.savefig(f"boxplot_dispersao_{tipo}.png", dpi=300)
    plt.close()

# Rodar para ambos
testar_dispersao(df_genes, "genes")
testar_dispersao(df_phages, "phages")
