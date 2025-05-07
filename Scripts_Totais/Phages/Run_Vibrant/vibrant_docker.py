import os
import subprocess

# Configurações do usuário
input_folder = "/home/rafael/Documentos/PhD/WGS_Bacterias_final/WGS_Bacterias/WGS_Listeriaceae/WGS_Listeria/WGS_Listeria_monocytogenes"  # Pasta com os arquivos de entrada
output_folder = "/home/rafael/Genomes/LVA/E_coli/Analises_Virais/Vibrant/Listeria"  # Correto: pasta onde os resultados serão salvos
threads = 64  # Número de threads a serem usadas

# Função para executar o VIBRANT no Docker
def run_vibrant(input_file):
    command = [
        "docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}",
        "-v", f"{input_folder}:/input",
        "-v", f"{output_folder}:/output",
        "staphb/vibrant",
        "VIBRANT_run.py",
        "-i", f"/input/{input_file}",
        "-folder", "/output",
        "-t", str(threads)
    ]
    print(f"Executando: {' '.join(command)}")
    subprocess.run(command, check=True)

# Processar todos os arquivos .fasta na pasta de entrada
for file_name in os.listdir(input_folder):
    if file_name.endswith(".fasta"):
        print(f"Processando arquivo: {file_name}")
        run_vibrant(file_name)

print("✅ Processamento concluído!")
