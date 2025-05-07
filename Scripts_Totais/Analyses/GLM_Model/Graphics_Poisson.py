import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# === Função genérica para gerar heatmaps ===
def gerar_heatmap(matrix, title, output_file, italic_x=False, italic_y=False):
    plt.figure(figsize=(18, 10))
    sns.set(style="white")

    ax = sns.heatmap(
        matrix,
        cmap="YlOrRd",
        linewidths=0.5,
        linecolor="lightgray",
        cbar_kws={'label': 'Fold-change (e^coef)'},
        square=False,
        annot=False
    )

    plt.title(title, fontsize=14, pad=20)

    if italic_x:
        ax.set_xticklabels(
            [label.get_text() for label in ax.get_xticklabels()],
            rotation=45, ha='right', fontsize=10, style='italic')
    else:
        ax.set_xticklabels(
            [label.get_text() for label in ax.get_xticklabels()],
            rotation=45, ha='right', fontsize=10)

    if italic_y:
        ax.set_yticklabels(
            [label.get_text() for label in ax.get_yticklabels()],
            rotation=0, fontsize=10, style='italic')

    plt.tight_layout()
    plt.savefig(output_file, dpi=600)
    print(f"✅ Heatmap salvo em: {output_file}")
    plt.close()

# === Diretórios de entrada ===
datasets = {
    "E.coli": "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Results/Poisson/Poisson_correct/E._coli",
    "Salmonella enterica": "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Results/Poisson/Poisson_correct/Salmonella_enterica",
    "Klebsiella pneumoniae": "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Results/Poisson/Poisson_correct/Klebsiella"
}

# === Definindo os thresholds específicos para cada bactéria ===
thresholds_per_bacteria = {
    "E.coli": {
        "phage-gene": 0.75,
        "gene-gene": 0.75,
        "phage-combo": 0.95,
        "phage-phage": 0.85
    },
    "Salmonella enterica": {
        "phage-gene": 0.70,
        "gene-gene": 0.75,
        "phage-combo": 0.80,
        "phage-phage": 0.50
    },
    "Klebsiella pneumoniae": {
        "phage-gene": 0.60,
        "gene-gene": 0.75,
        "phage-combo": 0.70,
        "phage-phage": 0.30
    }
}

for bacteria, base_dir in datasets.items():
    input_files = {
        "phage-gene": os.path.join(base_dir, "phage-gene.tsv"),
        "gene-gene": os.path.join(base_dir, "gene-gene.tsv"),
        "phage-combo": os.path.join(base_dir, "phage-combo.tsv"),
        "phage-phage": os.path.join(base_dir, "phage-phage.tsv")
    }

    for tipo, path in input_files.items():
        df = pd.read_csv(path, sep="\t", skip_blank_lines=True)
        df.columns = df.columns.str.strip()

        if tipo != "phage-phage":
            threshold = thresholds_per_bacteria[bacteria].get(tipo, 0.75)
            if "coef" in df.columns:
                df = df[df['coef'] >= df['coef'].quantile(threshold)]
                df['coef'] = np.exp(df['coef'])  # converte para escala de multiplicação

        if tipo == "phage-gene":
            df = df.dropna(subset=["gene", "phage", "coef"])
            df_grouped = df.groupby(["gene", "phage"], as_index=False).agg({"coef": "mean"})
            matrix = df_grouped.pivot(index="gene", columns="phage", values="coef")
            out_path = os.path.join(base_dir, "heatmap_phage_gene.png")
            title = f"Positive correlations between phage genera and resistance genes in $\\it{{{bacteria}}}$"
            if not matrix.empty:
                gerar_heatmap(matrix, title, out_path, italic_x=True, italic_y=False)
            else:
                print(f"⚠️  Matriz vazia para {tipo} ({bacteria}) — nada será plotado.")

        elif tipo == "gene-gene":
            df = df.dropna(subset=["gene1", "gene2", "coef"])
            df_grouped = df.groupby(["gene1", "gene2"], as_index=False).agg({"coef": "mean"})
            matrix = df_grouped.pivot(index="gene1", columns="gene2", values="coef")
            out_path = os.path.join(base_dir, "heatmap_gene_gene.png")
            title = f"Gene-gene correlations in $\\it{{{bacteria}}}$"
            if not matrix.empty:
                gerar_heatmap(matrix, title, out_path, italic_x=True, italic_y=True)
            else:
                print(f"⚠️  Matriz vazia para {tipo} ({bacteria}) — nada será plotado.")

        elif tipo == "phage-combo":
            df = df.dropna(subset=["gene", "predictors", "coef"])
            df["coef"] = np.exp(df["coef"])
            df[["phage1", "phage2"]] = df['predictors'].str.split('+', expand=True)
            df_grouped = df.groupby(["gene", "predictors"], as_index=False).agg({"coef": "mean"})
            matrix = df_grouped.pivot(index="gene", columns="predictors", values="coef")
            out_path = os.path.join(base_dir, "heatmap_phage_combo.png")
            title = f"Combo correlations with phage genera in $\\it{{{bacteria}}}$"
            if not matrix.empty:
                gerar_heatmap(matrix, title, out_path, italic_x=True, italic_y=False)
            else:
                print(f"⚠️  Matriz vazia para {tipo} ({bacteria}) — nada será plotado.")

        elif tipo == "phage-phage":
            threshold = thresholds_per_bacteria[bacteria].get(tipo, 0.75)
            df = df.dropna(subset=["phage1", "phage2", "odds_ratio"])
            df = df[df['odds_ratio'] >= df['odds_ratio'].quantile(threshold)]
            df_grouped = df.groupby(["phage1", "phage2"], as_index=False).agg({"odds_ratio": "mean"})
            matrix = df_grouped.pivot(index="phage1", columns="phage2", values="odds_ratio")
            out_path = os.path.join(base_dir, "heatmap_phage_phage.png")
            title = f"Phage-phage correlations in $\\it{{{bacteria}}}$"
            if not matrix.empty:
                gerar_heatmap(matrix, title, out_path, italic_x=True, italic_y=True)
            else:
                print(f"⚠️  Matriz vazia para {tipo} ({bacteria}) — nada será plotado.")
