#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script para unificar (merge) os arquivos FNA (fagos) usando apenas a coluna "scaffolds"
do TSV base (/home/rafael/final_results_corrected.tsv) para determinar quais nodes
representam o mesmo fago. O merge será realizado apenas se a coluna "scaffolds" contiver
mais de um node. Os arquivos FNA são buscados nas pastas "Lysogenic" e "Lytic", e o arquivo
unificado é salvo na pasta "Unified" dentro do diretório do sample, com o nome formado pela 
concatenação de todos os nodes unificados e o sample.
"""

import os
import glob
import pandas as pd

# Caminhos
tsv_path = "/home/rafael/salmonella_final_results_lysogenic_complemented.tsv"
base_fna_dir = "/home/rafael/Genomes/LVA/E_coli/Analises_Virais/Vibrant/Enterobacter/Salmonella/New_nodes/New_Nodes_Correct"

def find_fna_file(sample, node):
    """
    Procura o arquivo FNA correspondente ao node nas pastas Lysogenic e Lytic.
    Usa o padrão: {node}_*_{sample}.fna
    Retorna o primeiro arquivo encontrado ou None.
    """
    folders = ["Lysogenic", "Lytic"]
    for folder in folders:
        search_path = os.path.join(base_fna_dir, sample, folder, f"{node}_*_{sample}.fna")
        print(f"🔍 Procurando por: {search_path}")
        files = glob.glob(search_path)
        if files:
            print(f"📁 Arquivo encontrado: {files[0]}")
            return files[0]
    print(f"[!] Arquivo FNA não encontrado para node '{node}' no sample '{sample}'")
    return None

def merge_fna_files(file_list, output_path):
    """
    Concatena os conteúdos dos arquivos em file_list e salva em output_path.
    """
    with open(output_path, "w") as outfile:
        for f in file_list:
            with open(f, "r") as infile:
                outfile.write(infile.read())
    print(f"✅ Merged {len(file_list)} files into {output_path}")

def process_tsv_row(row):
    """
    Processa uma linha do TSV base usando a coluna "scaffolds".
    Se houver mais de um node (separados por ";"), realiza o merge dos arquivos FNA
    correspondentes e os salva na pasta "Unified" do sample.
    
    Retorna o caminho do arquivo unificado, ou None se houver apenas um node ou se 
    nenhum arquivo for encontrado.
    """
    sample = str(row['sample']).strip()
    scaffolds_str = row['scaffolds']  # Exemplo: "NODE_17_length_82844_cov_15.438676; NODE_16_length_86950_cov_14.563146; NODE_55_length_31773_cov_13.759559"
    scaffold_list = [x.strip() for x in scaffolds_str.split(";") if x.strip()]
    
    # Realiza merge somente se houver mais de um node
    if len(scaffold_list) < 2:
        print(f"Sample {sample} possui apenas 1 node ({scaffold_list[0]}). Ignorando merge.")
        return None
    
    files_to_merge = []
    for node in scaffold_list:
        fna = find_fna_file(sample, node)
        if fna:
            files_to_merge.append(fna)
    
    if not files_to_merge:
        print(f"[!] Nenhum arquivo encontrado para sample {sample} com nodes: {scaffold_list}")
        return None
    
    # Nome de saída: concatenação de todos os nodes (separados por '_') + sample
    out_base = "_".join(scaffold_list)
    output_filename = f"{out_base}_{sample}.fna"
    
    # Cria a pasta "Unified" dentro do sample, se não existir
    unified_dir = os.path.join(base_fna_dir, sample, "Unified")
    os.makedirs(unified_dir, exist_ok=True)
    output_path = os.path.join(unified_dir, output_filename)
    
    merge_fna_files(files_to_merge, output_path)
    return output_path

def main():
    # Carrega o TSV base
    df = pd.read_csv(tsv_path, sep="\t")
    print("Colunas do TSV:", df.columns.tolist())
    
    for index, row in df.iterrows():
        print(f"\nProcessando linha {index} para sample {row['sample']}")
        merged_file = process_tsv_row(row)
        if merged_file:
            print(f"Arquivo unificado criado: {merged_file}")
        else:
            print(f"Merge não realizado para a linha {index}.")
    
    print("\n✅ Merge de todos os grupos concluído.")

if __name__ == "__main__":
    main()
