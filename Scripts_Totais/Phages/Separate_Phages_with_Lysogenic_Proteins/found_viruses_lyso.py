import os
import csv
import re
import shutil  # Novo import necessário para mover arquivos

def extract_virus_genre(report_file):
    """
    Primeiro, percorre o arquivo em busca de um termo (nas linhas com rank 'G' ou 'S')
    que termine com "virus". Se encontrado, retorna imediatamente.
    
    Se nenhum termo for encontrado, realiza uma segunda passagem no arquivo
    procurando por linhas cujo termo comece com "unclassified". Se houver exatamente
    uma ocorrência, retorna esse termo; caso contrário, retorna None.
    """
    try:
        # Primeira passagem: busca por termo que termina com "virus" nas linhas de rank 'G' ou 'S'
        with open(report_file, 'r') as f:
            for line in f:
                columns = line.strip().split('\t')
                if len(columns) < 6:
                    continue
                rank = columns[3].strip()
                term = columns[5].strip()
                if rank in ['G', 'S']:
                    match = re.search(r'\b(\w+virus)\b', term, re.IGNORECASE)
                    if match:
                        return match.group(1)
        
        # Segunda passagem: busca por termos "unclassified" em qualquer linha
        unclassified_terms = []
        with open(report_file, 'r') as f:
            for line in f:
                columns = line.strip().split('\t')
                if len(columns) < 6:
                    continue
                term = columns[5].strip()
                if term.lower().startswith("unclassified"):
                    unclassified_terms.append(term)
        if len(unclassified_terms) == 1:
            return unclassified_terms[0]
    except Exception as e:
        print(f"Erro ao processar {report_file}: {e}")
    return None

def find_report_file(sample_dir, node, samplename):
    """
    Procura o arquivo correspondente ao NODE, verifica tamanho e move para 'Missing' se necessário.
    Retorna uma tupla (caminho_do_arquivo, mensagem_de_erro).
    """
    try:
        for file in os.listdir(sample_dir):
            if file.startswith(node) and file.endswith(f"_{samplename}.fna_kraken2_report.tsv"):
                file_path = os.path.join(sample_dir, file)
                
                # Verifica tamanho do arquivo
                if os.path.getsize(file_path) < 280:
                    # Cria diretório Missing se não existir
                    missing_dir = os.path.join(sample_dir, "Missing")
                    os.makedirs(missing_dir, exist_ok=True)
                    print(f"✅ Diretório 'Missing' criado/verificado em: {missing_dir}")
                    
                    # Move o arquivo
                    dest_path = os.path.join(missing_dir, file)
                    shutil.move(file_path, dest_path)
                    print(f"⚠️ Arquivo pequeno movido: {file} -> {dest_path}")
                    return (None, "Arquivo muito pequeno (movido para Missing)")
                
                return (file_path, None)  # Arquivo válido
        return (None, "Arquivo não encontrado")
    
    except Exception as e:
        print(f"Erro ao procurar arquivos em {sample_dir}: {e}")
        return (None, f"Erro: {str(e)}")

def process_lysogenic_nodes(base_dir, lysogenic_tsv, output_tsv):
    """
    Adiciona verificação de tamanho de arquivo e movimentação para pasta Missing.
    """
    results = []
    not_found_nodes = []
    
    with open(lysogenic_tsv, 'r') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        
        for row in reader:
            samplename = row['Samplename']
            node_fragment = row['NODE_Fragment']
            proteins = row['Proteins']
            node = node_fragment.split("_fragment_")[0]
            sample_dir = os.path.join(base_dir, samplename)
            
            # Busca arquivo com verificação de tamanho
            report_file, issue = find_report_file(sample_dir, node, samplename)
            
            if not report_file:
                print(f"❌ Problema no arquivo: {issue}")
                not_found_nodes.append((samplename, node_fragment, proteins, issue))
                virus_genre = 'Not Found'
            else:
                print(f"🔍 Processando arquivo: {report_file}")
                virus_genre = extract_virus_genre(report_file) or 'Not Found'
                if virus_genre == 'Not Found':
                    not_found_nodes.append((samplename, node_fragment, proteins, "Gênero não identificado"))
            
            results.append({
                'Samplename': samplename,
                'NODE_Fragment': node_fragment,
                'Proteins': proteins,
                'Virus_Genre': virus_genre
            })
    
    # Escreve os resultados no arquivo de saída
    with open(output_tsv, 'w', newline='') as tsvfile:
        writer = csv.DictWriter(tsvfile, fieldnames=['Samplename', 'NODE_Fragment', 'Proteins', 'Virus_Genre'], delimiter='\t')
        writer.writeheader()
        writer.writerows(results)
    
    # Salva os casos não encontrados em outro arquivo
    not_found_file = output_tsv.replace(".tsv", "_not_found.tsv")
    with open(not_found_file, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        writer.writerow(['Samplename', 'NODE_Fragment', 'Proteins', 'Issue'])
        writer.writerows(not_found_nodes)
    
    print(f"✅ Resultados salvos em {output_tsv}")
    print(f"🚨 Detalhes dos Not Found salvos em {not_found_file}")

# Caminhos principais
base_directory = "/home/rafael/Genomes/LVA/E_coli/Analises_Virais/Taxonomia/Blast/WGS_Bacteria/Enterobacter/E_coli/output_nodes/Kraken_Nodes/Nodes_all"
lysogenic_tsv = "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Phages/E.coli/Lysogenic_only/e_coli_lysogenic.phages.tsv"
output_tsv = "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Phages/E.coli/Lysogenic_only/E_coli_count_all_lysogenic.tsv"

# Executa o processamento
process_lysogenic_nodes(base_directory, lysogenic_tsv, output_tsv)
