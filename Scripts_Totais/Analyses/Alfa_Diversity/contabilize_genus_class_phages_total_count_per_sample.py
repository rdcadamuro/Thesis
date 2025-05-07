#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 15:25:36 2025

@author: rafael

Conta, para cada amostra, quantos fagos pertencem a cada GÊNERO (rank G)
e a cada ESPÉCIE (rank S/S1).

Entrada : /home/rafael/final_results_lysogenic_complemented.tsv
Kraken2 : …/Nodes_all/<sample>/<NODE‑prefix>*_kraken2_report.tsv

Saída   : imprime na tela e grava em
          /home/rafael/phage_genus_species_count.tsv

Colunas : Sample | Genus | Species | Count
"""

import os, glob, csv, pandas as pd
from collections import defaultdict

# ───── caminhos ─────
TSV_IN   = "/home/rafael/final_results_corrected_copy_.tsv"
KRAKEN   = "/home/rafael/Genomes/LVA/E_coli/Analises_Virais/Taxonomia/Blast/WGS_Bacteria/Enterobacter/E_coli/output_nodes/Kraken_Nodes/Nodes_all"
TSV_OUT  = "/home/rafael/E_coli_phage_genus_species_count.tsv"

CACHE = {}  # cache de relatórios lidos

# ─────────────────── utilidades ─────────────────── #
def node_prefix(node: str) -> str:
    """ devolve prefixo até '_cov_' inclusive """
    return node.split("_cov_")[0] + "_cov_" if "_cov_" in node else node

def find_report(sample: str, prefix: str) -> str | None:
    """procura relatório Kraken2 do prefixo numa sample"""
    sdir = os.path.join(KRAKEN, sample)
    subdirs = []
    if os.path.isdir(os.path.join(sdir, "Lysogenic")) or os.path.isdir(os.path.join(sdir, "Lytic")):
        if os.path.isdir(os.path.join(sdir, "Lysogenic")):
            subdirs.append(os.path.join(sdir, "Lysogenic"))
        if os.path.isdir(os.path.join(sdir, "Lytic")):
            subdirs.append(os.path.join(sdir, "Lytic"))
    else:
        subdirs.append(sdir)

    patterns = [
        f"{prefix}*_{sample}.fna_kraken2_report.tsv",
        f"{prefix}*.fna_kraken2_report.tsv",
        f"{prefix}*_{sample}.fna_kraken2_report.txt",
        f"{prefix}*.fna_kraken2_report.txt"
    ]

    for sub in subdirs:
        for pat in patterns:
            hit = glob.glob(os.path.join(sub, pat))
            if hit:
                return hit[0]
    return None

def parse_taxa(report: str) -> tuple[str|None, str|None]:
    """
    Lê um report Kraken2 e devolve (genus_name, species_name ou None).
    Se não houver rank G, retorna None.
    """
    if not report:
        return None, None
    if report in CACHE:
        return CACHE[report]

    g_name = s_name = None
    with open(report) as fh:
        for ln in fh:
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            rank = parts[3].strip()
            tax  = " ".join(parts[5:]).strip()
            if rank == "G" and tax.endswith("virus"):
                g_name = tax
            elif rank in ["S", "S1"] and g_name and tax.startswith(g_name):
                s_name = tax

    CACHE[report] = (g_name, s_name)
    return g_name, s_name

# ────────────────────── main ────────────────────── #
def main():
    df = pd.read_csv(TSV_IN, sep="\t")
    counts = defaultdict(int)  # (sample, genus, species) -> n

    for _, row in df.iterrows():
        sample = str(row["sample"]).strip()
        node   = row["scaffolds"].split(";")[0].strip()
        pref   = node_prefix(node)

        rpt = find_report(sample, pref)
        genus, species = parse_taxa(rpt)

        if genus:
            counts[(sample, genus, species or "")] += 1

    with open(TSV_OUT, "w", newline="") as w:
        writer = csv.writer(w, delimiter="\t")
        writer.writerow(["Sample", "Genus", "Species", "Count"])
        for (s, g, sp), n in sorted(counts.items(), key=lambda x: (x[0][0], x[0][1], -x[1])):
            writer.writerow([s, g, sp, n])

    print(f"✅ Resultado salvo em {TSV_OUT}")

if __name__ == "__main__":
    main()
