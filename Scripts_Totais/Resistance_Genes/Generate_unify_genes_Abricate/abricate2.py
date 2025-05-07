import os
import pandas as pd

def find_tsv_files_in_subdirs(root_dir):
    """Percorre todas as subpastas de Enterobacteriaceae para encontrar arquivos .tsv."""
    enterobacteriaceae_path = os.path.join(root_dir, "Enterobacteriaceae")
    if not os.path.exists(enterobacteriaceae_path):
        print(f"Diretório Enterobacteriaceae não encontrado em: {root_dir}")
        return

    for species_dir in ["E_coli", "Klebsiella", "Salmonella"]:
        species_path = os.path.join(enterobacteriaceae_path, species_dir)
        if not os.path.exists(species_path):
            print(f"Diretório {species_dir} não encontrado em: {enterobacteriaceae_path}")
            continue

        # Processa todos os arquivos .tsv nas subpastas da espécie
        for dirpath, _, filenames in os.walk(species_path):
            tsv_files = [f for f in filenames if f.endswith((".card.tsv", ".megares.tsv", ".resfinder.tsv", ".vfdb.tsv","_plasmidfinder.tsv"))]
            if tsv_files:
                print(f"Processando {len(tsv_files)} arquivos em: {dirpath}")
                process_tsv_files(dirpath, tsv_files)

def load_tsv_files(directory, files):
    """Carrega arquivos TSV com tratamento robusto de erros"""
    all_data = []
    
    for file in files:
        file_path = os.path.join(directory, file)
        try:
            # Verificação preliminar
            if not os.path.exists(file_path):
                print(f"Arquivo não encontrado: {file_path}")
                continue
                
            if os.path.getsize(file_path) == 0:
                print(f"Arquivo vazio: {file_path}")
                continue

            # Tentativa de leitura com múltiplos encodings
            encodings = ['utf-8', 'utf-8-sig', 'latin1']
            for encoding in encodings:
                try:
                    df = pd.read_csv(
                        file_path,
                        sep='\t',
                        encoding=encoding,
                        on_bad_lines='warn'  # Pula linhas problemáticas
                    )
                    break
                except UnicodeDecodeError:
                    continue
            else:
                print(f"Não foi possível ler {file_path} com os encodings testados")
                continue

            # Verifica se leu dados válidos
            if df.empty:
                print(f"Arquivo sem dados válidos: {file_path}")
                continue

            df["Source_File"] = file
            all_data.append(df)

        except Exception as e:
            print(f"Erro ao processar {file_path}: {str(e)}")
            # Debug avançado
            with open(file_path, 'rb') as f:
                print(f"Primeiros bytes: {f.read(100)}")

    return pd.concat(all_data, ignore_index=True) if all_data else pd.DataFrame()
    """Carrega todos os arquivos .tsv e adiciona uma coluna com o nome do arquivo."""
    all_data = []
    for file in files:
        file_path = os.path.join(directory, file)
        try:
            df = pd.read_csv(file_path, sep="\t")
            df["Source_File"] = file  # Adiciona nome do arquivo como coluna
            all_data.append(df)
        except Exception as e:
            print(f"Erro ao ler {file_path}: {e}")
    return pd.concat(all_data, ignore_index=True) if all_data else pd.DataFrame()

def filter_and_save(data, keyword, output_file):
    """Filtra os dados por banco de dados e salva em um arquivo unificado."""
    filtered_data = data[data['Source_File'].str.contains(keyword, case=False, regex=False)]
    if not filtered_data.empty:
        filtered_data.to_csv(output_file, sep="\t", index=False)
        print(f"Arquivo unificado salvo em: {output_file}")
    else:
        print(f"Nenhum dado encontrado para {keyword.upper()}.")

def process_tsv_files(directory, files):
    """Processa arquivos .tsv e cria versões unificadas por banco de dados."""
    output_path = os.path.join(directory, "Unified")
    os.makedirs(output_path, exist_ok=True)

    combined_data = load_tsv_files(directory, files)
    if combined_data.empty:
        print(f"Nenhum dado válido encontrado em: {directory}")
        return

    # Verifica colunas mínimas necessárias
    required_columns = ['Source_File', 'SEQUENCE', 'START', 'END', 'GENE', 'RESISTANCE']
    available_columns = [col for col in required_columns if col in combined_data.columns]
    if not available_columns:
        print(f"Colunas necessárias não encontradas em: {directory}")
        return

    combined_data_filtered = combined_data[available_columns]

    # Define caminhos de saída para cada banco de dados
    databases = {
        "card": "CARD_unified_results.tsv",
        "resfinder": "RESFINDER_unified_results.tsv",
        "vfdb": "VFDB_unified_results.tsv",
        "megares": "MEGARES_unified_results.tsv",
        "plasmidfinder": "PLASMIDFINDER_unified_results.tsv"
    }

    for db, output_file in databases.items():
        output_path_db = os.path.join(output_path, output_file)
        filter_and_save(combined_data_filtered, db, output_path_db)

# Caminho principal
root_directory = "/home/rafael/Genomes/LVA/E_coli/Results_Abricate_Ecoli"
find_tsv_files_in_subdirs(root_directory)