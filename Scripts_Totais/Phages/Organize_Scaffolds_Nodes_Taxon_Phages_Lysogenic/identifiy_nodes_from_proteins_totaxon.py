#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Agrupa fragments integrados pelo VIBRANT (Salmonella) por sobreposição + taxonomia
e grava em /home/rafael/salmonella_final_results_corrected.tsv.

⚠️  NOVO: se o diretório de um sample dentro de `kraken_dir` NÃO contiver
          sub‑pastas “Lysogenic” e “Lytic”, o script procura os relatórios
          Kraken2 diretamente dentro desse próprio diretório.
"""

import os, glob, pandas as pd
from collections import defaultdict

# ──────────────── Configurações ────────────────
vibrant_dir = "/home/rafael/Genomes/LVA/E_coli/Analises_Virais/Vibrant/Enterobacter/Klebsiella_p"
kraken_dir  = "/home/rafael/Genomes/LVA/E_coli/Analises_Virais/Vibrant/Enterobacter/Klebsiella_p/New_nodes/Kraken/Kraken_results"
output_file = "/home/rafael/Klebsiella_final_results_corrected.tsv"
gap_threshold = 500

rank_order = ['R','D','D1','K','P','C','C1','O','F','G','S','S1','U']
rank_index = {r:i for i,r in enumerate(rank_order)}

# ──────────── utilidades Kraken ────────────
def get_deepest_taxonomy(report_path: str) -> str|None:
    if not report_path or not os.path.exists(report_path):
        return None
    best, best_depth = None, -1
    with open(report_path) as fh:
        for ln in fh:
            cols = ln.rstrip('\n').split('\t')
            if len(cols) < 6:
                continue
            rank = cols[3].strip()
            depth = rank_index.get(rank, -1)
            if depth >= best_depth:
                best_depth = depth
                best = f"{cols[4].strip()} {cols[5].strip()}"
    return best

def find_kraken_file(sample: str, fragment: str) -> str|None:
    """
    Procura o relatório Kraken2 para *fragment*:
      1. em <sample>/Lysogenic/
      2. em <sample>/Lytic/
      3. diretamente em <sample>/
    Padrões testados:
      {fragment}_{sample}.fna_kraken2_report.tsv
      {fragment}.fna_kraken2_report.tsv
      {fragment}_kraken2_report.tsv
    """
    sample_dir = os.path.join(kraken_dir, sample)
    search_dirs = [
        os.path.join(sample_dir, "Lysogenic"),
        os.path.join(sample_dir, "Lytic"),
        sample_dir                         # ← procura aqui por último
    ]

    patterns = [
    f"{fragment}_{sample}.fna_kraken2_report.tsv",
    f"{fragment}.fna_kraken2_report.tsv",
    f"{fragment}_kraken2_report.tsv",
    f"{fragment}_{sample}.fna_kraken2_report.txt",
    f"{fragment}.fna_kraken2_report.txt",
    f"{fragment}_kraken2_report.txt"
    ]

    for sdir in search_dirs:
        if not os.path.isdir(sdir):
            continue
        for pat in patterns:
            hits = glob.glob(os.path.join(sdir, pat))
            if hits:
                return hits[0]
    return None

# ──────────────── processamento VIBRANT ────────────────
def process_sample(sample_path: str, sample: str) -> list[dict]:
    coord_tsv = os.path.join(sample_path,
                             f"VIBRANT_results_{sample}",
                             f"VIBRANT_integrated_prophage_coordinates_{sample}.tsv")
    if not os.path.exists(coord_tsv):
        print(f"[!] coordenadas não encontradas: {coord_tsv}")
        return []

    df = pd.read_csv(coord_tsv, sep='\t')
    if df.empty:
        return []

    # 1. grupos por overlap
    df = df.sort_values(['scaffold','nucleotide start']).reset_index(drop=True)
    groups, cur, chr_, end_ = [], [], None, -1e18
    for i,r in df.iterrows():
        s,e,c = r['nucleotide start'], r['nucleotide stop'], r['scaffold']
        if c == chr_ and s <= end_ + gap_threshold:
            cur.append(i); end_ = max(end_, e)
        else:
            if cur: groups.append(cur)
            cur, chr_, end_ = [i], c, e
    if cur: groups.append(cur)

    # 2. sub‑agrupamento por taxonomia
    out = []
    for idxs in groups:
        buckets = defaultdict(list)
        for i in idxs:
            frag = df.at[i,'fragment']
            rpt  = find_kraken_file(sample, frag)
            tax  = get_deepest_taxonomy(rpt) or "N/A"
            buckets[tax].append(i)

        for tax, ilist in buckets.items():
            starts = df.loc[ilist,'nucleotide start']
            ends   = df.loc[ilist,'nucleotide stop']
            scafs  = df.loc[ilist,'scaffold'].unique()
            out.append({
                "sample"     : sample,
                "scaffolds"  : "; ".join(scafs),
                "fragments"  : "; ".join(df.loc[ilist,'fragment']),
                "start"      : int(starts.min()),
                "end"        : int(ends.max()),
                "length"     : int(ends.max() - starts.min()),
                "taxonomy"   : tax,
                "overlapping": len(ilist) > 1
            })
    return out

# ──────────────────────────────── main ────────────────────────────────
def main():
    results = []
    for d in os.listdir(vibrant_dir):
        if not d.startswith("VIBRANT_"):
            continue
        sample = d.split("_",1)[1]
        print(f"▶ processando sample {sample}")
        results.extend(process_sample(os.path.join(vibrant_dir, d), sample))

    pd.DataFrame(results).to_csv(output_file, sep='\t', index=False)
    print(f"\n✅  Resultado salvo em {output_file} ({len(results)} linhas)")

if __name__ == "__main__":
    main()
