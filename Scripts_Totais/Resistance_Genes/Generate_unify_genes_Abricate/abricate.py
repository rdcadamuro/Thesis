import os
import subprocess
import multiprocessing

# Diretórios de entrada e saída
input_dirs = [
    "/home/rafael/Genomes/LVA/E_coli/Fastas_E_Coli",
]

output_base_dirs = [
    "/home/rafael/Genomes/LVA/E_coli/Results_Abricate_Ecoli",
    ]

# Parâmetros do Abricate
min_cov = 90
min_id = 90

# Lista de bancos de dados a serem processados (incluindo 'megares')
DATABASES = ['card', 'resfinder', 'vfdb', 'megares', 'plasmidfinder']

def create_directory(dir_path):
    """Cria um diretório se ele não existir."""
    os.makedirs(dir_path, exist_ok=True)

def should_process_sample(output_subdir, file_base):
    """Verifica se TODOS os arquivos de saída já existem."""
    return not all(
        os.path.exists(os.path.join(output_subdir, f"{file_base}.{db}.tsv"))
        for db in DATABASES
    )

def process_db(args):
    """Processa um único banco de dados para uma amostra com Abricate."""
    fasta_path, output_subdir, db = args
    file_base = os.path.splitext(os.path.basename(fasta_path))[0]
    output_path = os.path.join(output_subdir, f"{file_base}.{db}.tsv")

    try:
        cmd = f"abricate --db {db} --minid {min_id} --mincov {min_cov} --threads 20 {fasta_path} > {output_path}"
        subprocess.run(cmd, shell=True, check=True)
        print(f"[OK] {db.upper()} salvo em: {output_path}")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Falha ao processar {db.upper()} para {file_base}: {e}")

def process_fasta_files(input_dir, output_base_dir):
    """Processa todos os arquivos .fasta em paralelo."""
    args_list = []
    
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith(".fasta"):
                fasta_path = os.path.join(root, file)
                relative_path = os.path.relpath(root, input_dir)
                output_subdir = os.path.join(output_base_dir, relative_path)
                create_directory(output_subdir)
                file_base = os.path.splitext(os.path.basename(fasta_path))[0]

                # Se algum arquivo estiver faltando, adiciona TODOS os bancos para reprocessamento
                if should_process_sample(output_subdir, file_base):
                    for db in DATABASES:
                        args_list.append((fasta_path, output_subdir, db))

    # Processamento paralelo usando todos os cores disponíveis
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        pool.map(process_db, args_list)

if __name__ == "__main__":
    for input_dir, output_dir in zip(input_dirs, output_base_dirs):
        print(f"\nIniciando processamento para: {input_dir}")
        process_fasta_files(input_dir, output_dir)
        print(f"Concluído: {output_dir}\n")