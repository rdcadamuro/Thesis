#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Kraken2 Batch Runner com múltiplos diretórios de entrada e saída emparelhados,
recriando estrutura de subpastas com base no nome do sample.

@author: rafael
"""

import os
import subprocess

def run_kraken2_on_files(input_dirs, output_dirs, kraken_db):
    assert len(input_dirs) == len(output_dirs), "❌ As listas de entrada e saída devem ter o mesmo comprimento."

    for input_dir, output_dir in zip(input_dirs, output_dirs):
        if not os.path.exists(input_dir):
            print(f"❌ Diretório de entrada não encontrado: {input_dir}")
            continue

        for root, _, files in os.walk(input_dir):
            if not files:
                continue

            # Nome relativo da subpasta (ex: MS0002)
            rel_path = os.path.relpath(root, input_dir)
            output_subdir = os.path.join(output_dir, rel_path)
            os.makedirs(output_subdir, exist_ok=True)

            empty_entries = []

            for file in files:
                if file.endswith((".fna", ".fasta")):
                    input_file = os.path.join(root, file)
                    output_file = os.path.join(output_subdir, f"{file}_kraken2_output.txt")
                    report_file = os.path.join(output_subdir, f"{file}_kraken2_report.txt")

                    if os.path.exists(output_file) and os.path.exists(report_file):
                        print(f"⏩ Pulando {file}: já existem outputs.")
                    else:
                        kraken_command = [
                            "kraken2",
                            "--db", kraken_db,
                            "--output", output_file,
                            "--report", report_file,
                            "--threads", "64",
                            input_file
                        ]

                        try:
                            print(f"🚀 Kraken2 em {file} → {os.path.basename(output_subdir)}")
                            subprocess.run(kraken_command, check=True)
                            print(f"✅ Finalizado: {file}")
                        except subprocess.CalledProcessError as e:
                            print(f"❌ Erro Kraken2 em {file}: {e}")
                            continue

                    # Verifica se está vazio ou sem classificações
                    if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
                        empty_entries.append((file, "Saída vazia"))
                    else:
                        with open(output_file, "r") as f:
                            classified = [line for line in f if line.startswith("C")]
                            if not classified:
                                empty_entries.append((file, "Sem classificações (sem 'C')"))

            # Relatório de vazios/unclassified
            if empty_entries:
                tsv_path = os.path.join(output_subdir, "kraken2_empty_or_unclassified.tsv")
                with open(tsv_path, "w") as tsv:
                    tsv.write("Arquivo\tMotivo\n")
                    for fname, motivo in empty_entries:
                        tsv.write(f"{fname}\t{motivo}\n")
                print(f"📄 Relatório salvo: {tsv_path}")

# === 🧬 CONFIGURAÇÃO ===
input_dirs = [
    "/home/rafael/Genomes/LVA/E_coli/Analises_Virais/Vibrant/Enterobacter/Klebsiella_p/New_nodes/New_nodes_correct/MS6364"
]

output_dirs = [
    "/home/rafael/Genomes/LVA/E_coli/Analises_Virais/Vibrant/Enterobacter/Klebsiella_p/New_nodes/Kraken/Kraken_results/MS6364"
]

kraken_db = "/home/rafael/kraken2_db"

run_kraken2_on_files(input_dirs, output_dirs, kraken_db)
