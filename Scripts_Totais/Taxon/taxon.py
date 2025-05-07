import argparse
import os
import subprocess
from Bio import Entrez
import pandas as pd

# Configurações globais
Entrez.email = "cadamuro.rafael@gmail.com"
Entrez.api_key = "dad248454a51d9d946350764b56750f61a07"  # Sua chave API do NCBI

# Função para executar o BLAST
def executar_blast(arquivo_entrada, output_dir, threads=60):
    nome_amostra = os.path.basename(arquivo_entrada).split('.')[0]
    saida_blast = os.path.join(output_dir, f"resultados_blast_{nome_amostra}.txt")
    saida_formatada = os.path.join(output_dir, f"resultados_formatados_{nome_amostra}.tsv")

    # Checkpoints para verificação de existência dos arquivos
    if os.path.exists(saida_formatada):
        print(f"Resultados já processados para {nome_amostra}. Pulando reanálise...")
        return saida_formatada  # Retornar o caminho do .tsv gerado
    
    if os.path.exists(saida_blast):
        print(f"BLAST já foi executado para {nome_amostra}. Usando resultados existentes...")
    else:
        # Verificar se o arquivo de entrada está vazio ou corrompido
        if os.path.getsize(arquivo_entrada) == 0:
            print(f"Arquivo vazio: {arquivo_entrada}")
            return "corrupted"
        with open(arquivo_entrada, 'r') as f:
            if not f.readline().startswith('>'):
                print(f"Arquivo corrompido: {arquivo_entrada}")
                return "corrupted"
        
        # Executar o BLAST
        comando_blast = f"blastn -query {arquivo_entrada} -db ref_viruses_rep_genomes -out {saida_blast} " \
                        f"-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' " \
                        f"-evalue 1e-5 -num_threads {threads} -max_target_seqs 5"
        try:
            subprocess.run(comando_blast, shell=True, check=True)
        except subprocess.CalledProcessError:
            print(f"Erro ao executar BLAST para {nome_amostra}")
            return "error"

    # Confirmar que o arquivo de saída foi gerado
    if not os.path.exists(saida_blast):
        print(f"Erro: arquivo {saida_blast} não foi gerado.")
        return "error"

    return saida_blast

# Função para processar e formatar os resultados do BLAST
def processar_resultados_blast(saida_blast, output_dir):
    if saida_blast in ["corrupted", "error"]:
        return None
    
    nome_amostra = os.path.basename(saida_blast).split('_')[2].split('.')[0]
    saida_formatada = os.path.join(output_dir, f"resultados_formatados_{nome_amostra}.tsv")

    # Carregar e filtrar os dados do BLAST
    def carregar_e_filtrar_dados(saida_blast):
        df = pd.read_csv(saida_blast, sep='\t', header=None,
                         names=['qseqid', 'sseqid', 'similaridade', 'comprimento', 'mismatch', 'gapopen',
                                'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'stitle'])
        return df[(df['similaridade'] > 80) & (df['evalue'] < 1e-10)]

    # Salvar o arquivo formatado
    def salvar_resultados(df, output_file):
        df.to_csv(output_file, sep='\t', index=False)
        print(f"Resultados formatados salvos em {output_file}")

    # Carregar e filtrar os dados do BLAST
    df_filtrado = carregar_e_filtrar_dados(saida_blast)
    
    if df_filtrado.empty:
        print(f"Nenhum resultado válido encontrado em {saida_blast}")
        return "empty"

    # Salvar os resultados filtrados em formato .tsv
    salvar_resultados(df_filtrado, saida_formatada)

    # Confirmar que o arquivo .tsv foi gerado com sucesso
    if os.path.exists(saida_formatada):
        print(f"Arquivo {saida_formatada} gerado com sucesso.")
        return saida_formatada
    else:
        print(f"Erro ao gerar o arquivo {saida_formatada}")
        return "error"

# Função principal para processar todos os subdiretórios e arquivos .fna
def main():
    parser = argparse.ArgumentParser(description="Executar BLAST e processar resultados para sequências virais.")
    parser.add_argument("-i", "--input_dir", required=True, help="Diretório de entrada com subdiretórios contendo arquivos FASTA/FNA")
    parser.add_argument("-o", "--output_dir", required=True, help="Diretório de saída para resultados")
    parser.add_argument("-t", "--threads", type=int, default=64, help="Número de threads para o BLAST (padrão: 62)")
    
    args = parser.parse_args()
    
    # Criar diretório de saída se não existir
    os.makedirs(args.output_dir, exist_ok=True)

    # Arquivo de registro de amostras corrompidas ou com erro no BLAST
    corrupted_file = os.path.join(args.output_dir, "amostras_corrompidas_ou_com_erro.tsv")
    
    with open(corrupted_file, 'w') as f:
        f.write("amostra\tstatus\n")

    # Percorrer todos os subdiretórios dentro do diretório de entrada
    for root, dirs, files in os.walk(args.input_dir):
        for file in files:
            if file.endswith(('.fna', '.fasta')):
                arquivo_entrada = os.path.join(root, file)
                nome_amostra = root.split("/")[-2]  # Identifica a pasta da amostra (e.g., MS07253)
                tipo_seq = root.split("/")[-1]     # Identifica se é Lytic ou Lysogenic

                # Criar um diretório específico para cada combinação amostra/tipo no diretório de saída
                output_subdir = os.path.join(args.output_dir, nome_amostra, tipo_seq)
                os.makedirs(output_subdir, exist_ok=True)

                print(f"Processando {arquivo_entrada}...")

                # Executar o BLAST
                saida_blast = executar_blast(arquivo_entrada, output_subdir, args.threads)
                
                # Registro de arquivos corrompidos ou com erro no BLAST
                if saida_blast in ["corrupted", "error"]:
                    with open(corrupted_file, 'a') as f:
                        f.write(f"{nome_amostra}\t{tipo_seq}\t{saida_blast}\n")
                    continue

                # Processar os resultados do BLAST e verificar o .tsv gerado
                saida_formatada = processar_resultados_blast(saida_blast, output_subdir)

if __name__ == "__main__":
    main()
