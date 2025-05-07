import pandas as pd
import os

# Configurações para cada bactéria
bacteria_config = [
    {
        'name': 'E_coli',
        'base_dir': '/home/rafael/Genomes/LVA/E_coli/Analises_Virais/Abricate/Enterobacteriaceae/E_coli/Unified',
        'plasmidial_file': 'PLASMIDFINDER_unified_results.tsv',
        'file_names': [
            'CARD_unified_results.tsv',
            'RESFINDER_unified_results.tsv',
            'VFDB_unified_results.tsv',
            'MEGARES_unified_results.tsv',
            'Bacmet_unified_results.tsv',
            'PLASMIDFINDER_unified_results.tsv'
        ]
    },
    {
        'name': 'Salmonella',
        'base_dir': '/home/rafael/Genomes/LVA/E_coli/Analises_Virais/Abricate/Enterobacteriaceae/Salmonella/Unified',
        'plasmidial_file': 'PLASMIDFINDER_unified_results.tsv',
        'file_names': [
            'CARD_unified_results.tsv',
            'RESFINDER_unified_results.tsv',
            'VFDB_unified_results.tsv',
            'MEGARES_unified_results.tsv',
            'Bacmet_unified_results.tsv',
            'PLASMIDFINDER_unified_results.tsv'
        ]
    },
    {
        'name': 'Klebsiella_pneumoniae',
        'base_dir': '/home/rafael/Genomes/LVA/E_coli/Analises_Virais/Abricate/Enterobacteriaceae/Klebsiella/Klebsiella_pneumoniae/Unified',
        'plasmidial_file': 'PLASMIDFINDER_unified_results.tsv',
        'file_names': [
            'CARD_unified_results.tsv',
            'RESFINDER_unified_results.tsv',
            'VFDB_unified_results.tsv',
            'MEGARES_unified_results.tsv',
            'Bacmet_unified_results.tsv',
            'PLASMIDFINDER_unified_results.tsv'
        ]
    }
]

# 2. Criar um dicionário inicial com categorias
gene_to_group = {
   "mdfA": "Efflux pumps",
    "marA": "Efflux pumps",
    "acrB": "Efflux pumps",
    "acrS": "Efflux pumps",
    "acrE": "Efflux pumps",
    "acrF": "Efflux pumps",
    "mdtP": "Efflux pumps",
    "mdtO": "Efflux pumps",
    "mdtN": "Efflux pumps",
    "oqxA": "Efflux pumps",
    "mdtE": "Efflux pumps",
    "mdtF": "Efflux pumps",
    "mdtK": "Efflux pumps",
    "mdtM": "Efflux pumps",
    "mdtA": "Efflux pumps",
    "mdtB": "Efflux pumps",
    "mdtC": "Efflux pumps",
    "yojI": "Efflux pumps",
    "emrY": "Efflux pumps",
    "emrK": "Efflux pumps",
    "emrA": "Efflux pumps",
    "emrB": "Efflux pumps",
    "acrD": "Efflux pumps",
    "Escherichia_coli_emrE": "Efflux pumps",
    "Escherichia_coli_acrA": "Efflux pumps",
    "tolC": "Efflux pumps",
    "msbA": "Efflux pumps",
    "mdtG": "Efflux pumps",
    "mdtH": "Efflux pumps",
    "TEM-1": "Beta-lactam resistance",
    "CTX-M-15": "Beta-lactam resistance",
    "CTX-M-14": "Beta-lactam resistance",
    "CTX-M-55": "Beta-lactam resistance",
    "CTX-M-32": "Beta-lactam resistance",
     "CTX-M-65": "Beta-lactam resistance",
    "CTX-M-1": "Beta-lactam resistance",
    "TEM-135": "Beta-lactam resistance",
    "FosB": "Fosfomycin resistance",
    "BLA1": "Beta-lactam resistance",
    "dfrA21": "Trimethoprim resistance",
    "aadA4": "Aminoglycoside resistance",
    "TEM-116": "Beta-lactam resistance",
    "TEM-143": "Beta-lactam resistance",
    "AAC(2')-Ia": "Aminoglycoside resistance",
    "catIII": "Phenicol resistance",
    "AAC(6')-Ib7": "Aminoglycoside resistance",
    "Salmonella_enterica_cmlA": "Phenicol resistance",
    "aadA24": "Aminoglycoside resistance",
    "OXA-2": "Beta-lactam resistance",
    "aadA17": "Aminoglycoside resistance",
    "CTX-M-27": "Beta-lactam resistance",
    "TEM-33": "Beta-lactam resistance",
    "TEM-52": "Beta-lactam resistance",
    "blaTEM-104_1": "Beta-lactam resistance",
    "blaCTX-M-1_1": "Beta-lactam resistance",
    "tet(B)_1": "Tetracycline resistance",
    "lnu(G)_1": "Lincosamide resistance",
    "mph(G)_1": "Macrolide resistance",
    "blaNDM-4_1": "Carbapenem resistance",
    "erm(B)_18": "Macrolide resistance",
    "sul2_14": "Sulfonamide resistance",
    "mcr-1.11_1": "Colistin resistance",
    "tet(M)_5": "Tetracycline resistance",
    "qnrB7_1": "Quinolone resistance",
    "blaOXA-2_1": "Beta-lactam resistance",
    "aph(6)-Id_4": "Aminoglycoside resistance",
    "mcr-2.1_1": "Colistin resistance",
    "aph(3'')-Ib_4": "Aminoglycoside resistance",
    "blaCTX-M-55_1": "Beta-lactam resistance",
    "espY2": "Type III Secretion System (T3SS)",
    "shuA": "Shigella toxin",
    "hlyC": "Hemolysin",
    "gspE": "Type III Secretion System (T3SS)",
    "gspF": "Type III Secretion System (T3SS)",
    "gspG": "Type III Secretion System (T3SS)",
    "gspH": "Type III Secretion System (T3SS)",
    "gspI": "Type III Secretion System (T3SS)",
    "gspJ": "Type III Secretion System (T3SS)",
    "espP": "Type III Secretion System (T3SS)",
    "cesT": "Cecropin resistance",
    "eae": "Adhesion",
    "escD": "Type III Secretion System (T3SS)",
    "sepL": "Type III Secretion System (T3SS)",
    "espA": "Type III Secretion System (T3SS)",
    "espD": "Type III Secretion System (T3SS)",
    "espB": "Type III Secretion System (T3SS)",
    "cesD2": "Cecropin resistance",
    "escF": "Type III Secretion System (T3SS)",
    "escG": "Type III Secretion System (T3SS)",
    "cesL": "Cecropin resistance",
    "escV": "Type III Secretion System (T3SS)",
    "escN": "Type III Secretion System (T3SS)",
    "escO": "Type III Secretion System (T3SS)",
    "escP": "Type III Secretion System (T3SS)",
    "sepQ/escQ": "Type III Secretion System (T3SS)",
    "espH": "Type III Secretion System (T3SS)",
    "nleE": "Type III Secretion System (T3SS)",
    "nleB1": "Type III Secretion System (T3SS)",
    "cesD": "Cecropin resistance",
    "escC": "Type III Secretion System (T3SS)",
    "sepD": "Type III Secretion System (T3SS)",
    "escJ": "Type III Secretion System (T3SS)",
    "escI": "Type III Secretion System (T3SS)",
    "sepZ/espZ": "Type III Secretion System (T3SS)",
    "stx2B": "Shigella toxin",
    "stx2A": "Shigella toxin",
    "east1": "Type III Secretion System (T3SS)",
    "nleF": "Type III Secretion System (T3SS)",
    "nleH2": "Type III Secretion System (T3SS)",
    "map": "Motility",
    "tir": "Type III Secretion System (T3SS)",
    "espX2": "Type III Secretion System (T3SS)",
    "nleD": "Type III Secretion System (T3SS)",
    "nleH1": "Type III Secretion System (T3SS)",
    "nleC": "Type III Secretion System (T3SS)",
    "nleB2": "Type III Secretion System (T3SS)",
    "espG": "Type III Secretion System (T3SS)",
    "escT": "Type III Secretion System (T3SS)",
    "escS": "Type III Secretion System (T3SS)",
    "espJ": "Type III Secretion System (T3SS)",
    "nleA/espI": "Type III Secretion System (T3SS)",
    "espY2": "Type III Secretion System (T3SS)",
    "cesAB": "Cecropin resistance",
    "escL": "Type III Secretion System (T3SS)",
    "paa": "Pathogenicity",
    "escE": "Type III Secretion System (T3SS)",
    "espR3": "Type III Secretion System (T3SS)",
    "espR4": "Type III Secretion System (T3SS)",
    "espY3": "Type III Secretion System (T3SS)",
    "espK": "Type III Secretion System (T3SS)",
    "espX6": "Type III Secretion System (T3SS)",
    "espW": "Type III Secretion System (T3SS)",
    "espM2": "Type III Secretion System (T3SS)",
    "shuY": "Shigella toxin",
    "shuX": "Shigella toxin",
    "shuT": "Shigella toxin",
    "shuA": "Shigella toxin",
    "shuS": "Shigella toxin",
    "stcE": "Pathogenicity",
    "espM1": "Type III Secretion System (T3SS)",
    "stxA": "Shigella toxin",
    "stx1B": "Shigella toxin",
    "espY4": "Type III Secretion System (T3SS)",
       "TEM-176": "Beta-lactam resistance",
    "AAC(2')-IIa": "Aminoglycoside resistance",
    "TEM-168": "Beta-lactam resistance",
    "aadA14": "Aminoglycoside resistance",
    "dfrA5_1": "Trimethoprim resistance",
    "aadA17_1": "Aminoglycoside resistance",
    "fosB1_1": "Fosfomycin resistance",
    "tet(C)_2": "Tetracycline resistance",
    "blaTEM-106_1": "Beta-lactam resistance",
    "sul1_9": "Sulfonamide resistance",
    "blaCTX-M-32_2": "Beta-lactam resistance",
    "blaTEM-150_1": "Beta-lactam resistance",
    "blaCTX-M-32_2": "Beta-lactam resistance",
    "blaTEM-116_1": "Beta-lactam resistance",
    "dfrA8_1": "Trimethoprim resistance",
    "tet(H)_2": "Tetracycline resistance",
    "hugA_1": "Shigella toxin",
    "tet(M)_7": "Tetracycline resistance",
    "erm(B)_1": "Macrolide resistance",
    "blaTEM-52C_2": "Beta-lactam resistance",
    "aadA22_1": "Aminoglycoside resistance",
    "mcr-4.1_1": "Colistin resistance",
    "aadA12_1": "Aminoglycoside resistance",
    "dfrA21_1": "Trimethoprim resistance",
    "mph(A)_1": "Macrolide resistance",
    "blaOXA-162_1": "Beta-lactam resistance",
    "blaCTX-M-65_1": "Beta-lactam resistance",
    "blaCTX-M-32_2": "Beta-lactam resistance",
    "aadA14_1": "Aminoglycoside resistance",
    "aph(3')-Ia_10": "Aminoglycoside resistance",
    "sul2_14": "Sulfonamide resistance",
    "faeF": "Adhesion and biofilm formation",
    "faeH": "Adhesion and biofilm formation",
    "faeI": "Adhesion and biofilm formation",
    "faeJ": "Adhesion and biofilm formation",
    "cesF": "Cecropin resistance",
    "nleG7'": "Type III Secretion System (T3SS)",
    "eltA": "Type III Secretion System (T3SS)",
    "eltB": "Type III Secretion System (T3SS)",
    "f17d-C": "Adhesion",
    "f17d-D": "Adhesion",
    "f17d-A": "Adhesion",
    "afaF-VII": "Adhesion",
    "afaA-VIII": "Adhesion",
    "afaB-VIII": "Adhesion",
    "cnf1": "Cytotoxic necrotizing factor",
    "cdtC": "Cytotoxin",
    "cdtA": "Cytotoxin",
    "espF": "Type III Secretion System (T3SS)",
    "etgA": "Pathogenicity",
    "escU": "Type III Secretion System (T3SS)",
    "escR": "Type III Secretion System (T3SS)",
    "cif": "Cytotoxic necrotizing factor",
    "senB": "Adhesion",
    "kpsT": "Capsule biosynthesis",
    "gtrA": "Glycosyl transferase",
    "gtrB": "Glycosyl transferase",

    "nleG7": "Type III Secretion System (T3SS)",
    "Escherichia_coli_ampC": "Beta-lactam resistance",
    "OXA-1": "Beta-lactam resistance",
    "OXA-48": "Beta-lactam resistance",
    "OXA-162": "Beta-lactam resistance",
    "NDM-4": "Beta-lactam resistance",
    "Escherichia_coli_ampC1_beta-lactamase": "Beta-lactam resistance",
    "Escherichia_coli_ampH": "Beta-lactam resistance",
    "ErmB": "Macrolide resistance",
     "blaTEM-168_1": "Beta-lactam resistance",
    "blaCTX-M-27_1": "Beta-lactam resistance",
    "sul2_6": "Sulfonamide resistance",
    "aadA1_2": "Aminoglycoside resistance",
    "blaTEM-143_1": "Beta-lactam resistance",
    "blaTEM-33_1": "Beta-lactam resistance",
    "blaTEM-135_1": "Beta-lactam resistance",
    "blaTEM-122_1": "Beta-lactam resistance",
    "blaTEM-176_1": "Beta-lactam resistance",
    "mcr-3.2_1": "Colistin resistance",
    "aadA1_5": "Aminoglycoside resistance",
    "aac(2')-Ia_1": "Aminoglycoside resistance",
    "catA3_1": "Phenicol resistance",
    "sul1_2": "Sulfonamide resistance",
    "aadA1_3": "Aminoglycoside resistance",
    "lnu(F)_3": "Lincosamide resistance",
    "aadA2_2": "Aminoglycoside resistance",
    "aac(2')-IIa_1": "Aminoglycoside resistance",
    "aadA4_1": "Aminoglycoside resistance",
    "tet(M)_2": "Tetracycline resistance",
    "f17d-G": "Adhesion",
    "espL2": "Type III Secretion System (T3SS)",
    "ibeA": "Invasion",
    "espN": "Type III Secretion System (T3SS)",
    "espX7/nleL": "Type III Secretion System (T3SS)",
    "fliS": "Flagella synthesis",
    "nleA": "Type III Secretion System (T3SS)",
    "estIa": "Type III Secretion System (T3SS)",
    "sfaH": "Adhesion",
    "sfaS": "Adhesion",
    "sfaG": "Adhesion",
    "sfaF": "Adhesion",
    "sfaE": "Adhesion",
    "afaC-VIII": "Adhesion",
    "afaD-VIII": "Adhesion",
    "afaE-VIII": "Adhesion",
    "capE": "Capsule biosynthesis",
    "dep/capD": "Capsule biosynthesis",
    "capA": "Capsule biosynthesis",
    "capC": "Capsule biosynthesis",
    "capB": "Capsule biosynthesis",
    "cya": "Cytotoxin",
    "pagA": "Invasion",
    "lef": "Leucocidin resistance",
    "nheA": "Type III Secretion System (T3SS)",
    "espFu/tccP": "Type III Secretion System (T3SS)",
    "afaA-VIII": "Adhesion",
    "afaB-VIII": "Adhesion",
    "afaF-VII": "Adhesion",
    "ibaA": "Invasion",
    "mdtJ/ebrB/ydgF": "Efflux pumps",
    "fetB/ybbM": "Iron acquisition",
    "zitB/ybgR": "Zinc transporter",
    "ydeI": "Transporter",
    "mdtA/yegM": "Efflux pumps",
    "rcnA/yohM": "Cobalt and nickel transporter",
    "nikA": "Nickel transporter",
    "nikE": "Nickel transporter",
    "gadW/yhiW": "Acid resistance",
    "ibpB": "Stress response",
    "hdeB/yhiC": "Acid resistance",
    "gadE/yhiE": "Acid resistance",
    "fecE": "Iron uptake",
    "asr": "Acid stress response",
    "fetA/ybbL": "Iron uptake",
    "terB": "Tellurite resistance",
    "terD": "Tellurite resistance",
    "merR1": "Mercury resistance",
    "recG": "DNA repair",
    "copR": "Copper resistance",
    "merP": "Mercury resistance",
    "fliG": "Flagella synthesis",
    "gtrA": "Glycosyl transferase",
    "gtrB": "Glycosyl transferase",
    "gtrII": "Glycosyl transferase",
    "sfaA": "Adhesion",
    "afaC-I": "Adhesion",
    "afaB-I": "Adhesion",
    "afaD": "Adhesion",
    "aadA9_1": "Aminoglycoside resistance",
    "efa1": "Stress response",
    "afaC-VII": "Adhesion",
    "ospG": "Toxin production",
    "toxB": "Toxin production",
    "tcpC": "Toxin production",
    "inhA": "Inhibition of cell wall synthesis",
    "etpB": "Type III Secretion System (T3SS)",
    "afaA": "Adhesion",
    "afaD-VII": "Adhesion",
    "afaG-VII": "Adhesion",
    "afaE-VII": "Adhesion",
    "faeG": "Adhesion",
    "afaE-V": "Adhesion",
    "mdtI/ydgE": "Efflux pumps",
    "nikD": "Nickel transporter",
  "nikB": "Nickel transporter",
  "yjaA": "Ion transporter",
  "nikC": "Nickel transporter",
  "BAC0741|tcrA|tr|Q3MNJ6|Q3MNJ6_ENTFC": "Disinfectant resistance",  # Corrigido
  "hdeA/yhiB": "Acid resistance",
  "rcnB/yohN": "Cobalt and nickel transporter",
  "qacC/qacD/smr": "Disinfectant resistance",  # Correto
  "merR2": "Mercury resistance",
  "silR": "Silver resistance",
  "silF": "Silver resistance",
  "tcrB": "Copper resistance",
  "bac0708|chtS|ZP_06677962|sensor": "Metal ion resistance",  # Corrigido
    "draP": "Acid resistance",
    "daaF": "Adhesion",
     "mphG": "Macrolide resistance",
    "msrE": "Macrolide resistance",
    "bepD": "Efflux pumps",
    "eefX": "Efflux pumps",
    "adeT2": "Efflux pumps",
    "cepA": "Beta-lactam resistance",
    "abeM": "Efflux pumps",
      "adeN": "Efflux pumps",
    "adeB": "Efflux pumps",
    "adeG": "Efflux pumps",
    "BAC0708|chtS|ZP_06677962|sensor": "Environmental stress sensor",
    "aadA12": "Aminoglycoside resistance",
    "aadA5": "Aminoglycoside resistance",
    "amvA": "Efflux pumps",
    "aadA22": "Aminoglycoside resistance",
    "APH(3'')-Ib": "Aminoglycoside resistance",
    "APH(6)-Id": "Aminoglycoside resistance",
    "ANT(3'')-IIa": "Aminoglycoside resistance",
    "tet(A)": "Tetracycline resistance",
    "tet(B)": "Tetracycline resistance",
    "tet(C)": "Tetracycline resistance",
    "tet(D)": "Tetracycline resistance",
    "tet(H)": "Tetracycline resistance",
  
    "eptA": "Polymyxin resistance",
    "ugd": "Polymyxin resistance",
    "pmrF": "Polymyxin resistance",
    "dfrA1": "Folate synthesis inhibitors",
    "dfrA5": "Folate synthesis inhibitors",
    "dfrA8": "Folate synthesis inhibitors",
    "dfrA17": "Folate synthesis inhibitors",
    "cmlA1": "Phenicol resistance",
    "catI": "Phenicol resistance",
    "catB2_1": "Phenicol resistance",
    "catB3_2": "Phenicol resistance",
   "floR_2": "Phenicol resistance",
    "catA1_1": "Phenicol resistance",
    "catA2_1": "Phenicol resistance",
    "cmlA1_1": "Phenicol resistance",
    "linG": "Lincosamide resistance",
    "H-NS": "Global regulators",
    "CRP": "Global regulators",
    "baeS": "Two-component systems",
    "baeR": "Two-component systems",
    "cpxA": "Two-component systems",
    "kdpE": "Two-component systems",
    "mdfA": "Efflux pumps",
     "ErmB": "Macrolide resistance",
        "OXA-1": "Beta-lactam resistance",
     "TEM-1": "Beta-lactam resistance",
    "CTX-M-15": "Beta-lactam resistance",
    "tet(A)": "Tetracycline resistance",
    "tet(B)": "Tetracycline resistance",
    "blaCTX-M-14": "Beta-lactam resistance",
    # Efflux pumps
    "mdfA": "Efflux pumps",
   "acrS": "Efflux pumps",
    "acrE": "Efflux pumps",
    "acrF": "Efflux pumps",
    "mdtP": "Efflux pumps",
    "mdtO": "Efflux pumps",
    "mdtN": "Efflux pumps",
    "oqxA": "Efflux pumps",
    "mdtE": "Efflux pumps",
    "mdtF": "Efflux pumps",
    "mdtK": "Efflux pumps",
    # Beta-lactam resistance
    "TEM-1": "Beta-lactam resistance",
    "CTX-M-15": "Beta-lactam resistance",
    "CTX-M-14": "Beta-lactam resistance",
    "CTX-M-55": "Beta-lactam resistance",
    "CTX-M-32": "Beta-lactam resistance",
    "Escherichia_coli_ampC": "Beta-lactam resistance",
    "OXA-1": "Beta-lactam resistance",
    "OXA-48": "Beta-lactam resistance",
    "OXA-162": "Beta-lactam resistance",
    "NDM-4": "Beta-lactam resistance",
    # Macrolide resistance
    "ErmB": "Macrolide resistance",
      "mphG": "Macrolide resistance",
    "msrE": "Macrolide resistance",
    # Aminoglycoside resistance
    "aadA12": "Aminoglycoside resistance",
    "aadA5": "Aminoglycoside resistance",
    "aac(3)-IV": "Aminoglycoside resistance",
    "ANT(3'')-IIa": "Aminoglycoside resistance",
    # Tetracycline resistance
    "tet(A)": "Tetracycline resistance",
    "tet(B)": "Tetracycline resistance",
    "tet(C)": "Tetracycline resistance",
    "tet(D)": "Tetracycline resistance",
    "tet(H)": "Tetracycline resistance",
    # Polymyxins
    "MCR-1": "Colistin resistance",
    "MCR-2": "Colistin resistance",
    "MCR-3.2": "Colistin resistance",
    "MCR-4": "Colistin resistance",
    "MCR-4.2": "Colistin resistance",
    # Fosfomycin
    "FosA7": "Fosfomycin resistance",
    "FosA6": "Fosfomycin resistance",
    # Bacitracin
    "bacA": "Bacitracin resistance",
       # Reguladores e outros
    "sdiA": "Quorum sensing",
    "CRP": "Global regulators",
    "gadW": "Global regulators",
    "gadX": "Global regulators",
    "ramA": "Global regulators",
    # Outros grupos
    "qacH": "Disinfectant resistance",
     "cmlA1": "Phenicol resistance",
    "lnuF": "Lincosamide resistance",
    "lnuG": "Lincosamide resistance",
    # Quinolones
    "QnrS1": "Quinolone resistance",
    "QnrB5": "Quinolone resistance",
    "QnrB58": "Quinolone resistance",
    "QnrD1": "Quinolone resistance",
    # Bleomycin
    "determinant_of_bleomycin_resistance": "Bleomycin resistance",
    # Efflux pumps
    "mdfA": "Efflux pumps",
       "acrS": "Efflux pumps",
    "acrE": "Efflux pumps",
    "acrF": "Efflux pumps",
    "mdtP": "Efflux pumps",
    "mdtO": "Efflux pumps",
    "mdtN": "Efflux pumps",
    "oqxA": "Efflux pumps",
    "mdtE": "Efflux pumps",
    "mdtF": "Efflux pumps",
    "mdtK": "Efflux pumps",
    "mdtM": "Efflux pumps",
    "mdtA": "Efflux pumps",
    "mdtB": "Efflux pumps",
    "mdtC": "Efflux pumps",
    "yojI": "Efflux pumps",
    "emrY": "Efflux pumps",
    "emrK": "Efflux pumps",
    "emrB": "Efflux pumps",
    "acrD": "Efflux pumps",
    "Escherichia_coli_emrE": "Efflux pumps",
    "Escherichia_coli_acrA": "Efflux pumps",
     "msbA": "Efflux pumps",
    "mdtG": "Efflux pumps",
    "mdtH": "Efflux pumps",
    # Beta-lactam resistance
    "TEM-1": "Beta-lactam resistance",
    "CTX-M-15": "Beta-lactam resistance",
    "CTX-M-14": "Beta-lactam resistance",
    "CTX-M-55": "Beta-lactam resistance",
    "CTX-M-32": "Beta-lactam resistance",
    "Escherichia_coli_ampC": "Beta-lactam resistance",
    "OXA-1": "Beta-lactam resistance",
    "OXA-48": "Beta-lactam resistance",
    "OXA-162": "Beta-lactam resistance",
    "NDM-4": "Beta-lactam resistance",
    "Escherichia_coli_ampC1_beta-lactamase": "Beta-lactam resistance",
    "Escherichia_coli_ampH": "Beta-lactam resistance",
    # Macrolide resistance
    "ErmB": "Macrolide resistance",
     "mphG": "Macrolide resistance",
    "msrE": "Macrolide resistance",
    # Aminoglycoside resistance
   
    "aadA9": "Aminoglycoside resistance",
    
    # Tetracycline resistance
    "tet(A)": "Tetracycline resistance",
    "tet(B)": "Tetracycline resistance",
    "tet(C)": "Tetracycline resistance",
    "tet(D)": "Tetracycline resistance",
      "tet(H)": "Tetracycline resistance",
    # Polymyxins
 
    "eptA": "Polymyxin resistance",
    "ugd": "Polymyxin resistance",
    "pmrF": "Polymyxin resistance",
    # Folate synthesis inhibitors
    "dfrA1": "Folate synthesis inhibitors",
    "dfrA5": "Folate synthesis inhibitors",
    "dfrA8": "Folate synthesis inhibitors",
    "dfrA17": "Folate synthesis inhibitors",
    # Phenicol resistance
     "cmlA1": "Phenicol resistance",
    "catI": "Phenicol resistance",
    # Lincosamide resistance
    "linG": "Lincosamide resistance",
    # Reguladores e outros
    "H-NS": "Global regulators",
    "CRP": "Global regulators",
    "baeS": "Two-component systems",
    "baeR": "Two-component systems",
    "cpxA": "Two-component systems",
    "kdpE": "Two-component systems",
    # ... (outros genes já categorizados)
    "evgA": "Two-component systems",
    "evgS": "Two-component systems",
    "emrR": "Efflux pumps",
    "Escherichia_coli_mdfA": "Efflux pumps",
    # Efflux pumps
    "mdfA": "Efflux pumps",
    "acrS": "Efflux pumps",
    "acrE": "Efflux pumps",
    "acrF": "Efflux pumps",
    "mdtP": "Efflux pumps",
    "mdtO": "Efflux pumps",
    "mdtN": "Efflux pumps",
    "oqxA": "Efflux pumps",
    "mdtE": "Efflux pumps",
    "mdtF": "Efflux pumps",
    "mdtK": "Efflux pumps",
    "mdtM": "Efflux pumps",
    "mdtA": "Efflux pumps",
    "mdtB": "Efflux pumps",
    "mdtC": "Efflux pumps",
    "yojI": "Efflux pumps",
    "emrY": "Efflux pumps",
    "emrK": "Efflux pumps",
     "emrB": "Efflux pumps",
    "acrD": "Efflux pumps",
    "Escherichia_coli_emrE": "Efflux pumps",
    "Escherichia_coli_acrA": "Efflux pumps",
     "msbA": "Efflux pumps",
    "mdtG": "Efflux pumps",
    "mdtH": "Efflux pumps",
    "Staphylococcys_aureus_LmrS": "Efflux pumps",
    "Staphylococcus_aureus_norA": "Efflux pumps",
    "mepA": "Efflux pumps",
    # Macrolide resistance
    "ErmB": "Macrolide resistance",
    "mphC": "Macrolide resistance",
     # Aminoglycoside resistance
    "aadA12": "Aminoglycoside resistance",
    "aadA22": "Aminoglycoside resistance",
    "APH(3'')-Ib": "Aminoglycoside resistance",
    "APH(6)-Id": "Aminoglycoside resistance",
    "APH(3')-IIIa": "Aminoglycoside resistance",
    "ANT(3'')-IIa": "Aminoglycoside resistance",
    "ANT(4')-Ib": "Aminoglycoside resistance",
    "AAC(6')-Ie-APH(2'')-Ia": "Aminoglycosides",
    "SAT-4": "Aminoglycoside resistance",
    # Tetracycline resistance
    "tet(A)": "Tetracycline resistance",
    "tet(B)": "Tetracycline resistance",
    "tet(C)": "Tetracycline resistance",
    "tet(D)": "Tetracycline resistance",
     "tet(H)": "Tetracycline resistance",
    "tet(K)": "Tetracycline resistance",
    "tet(38)": "Tetracycline resistance",
    # Beta-lactam resistance
    "PC1_Beta-lactam resistance_(blaZ)": "Beta-lactam resistance",
       # Fosfomycin
    "FosB1": "Fosfomycin resistance",
    "Staphylococcus_aureus_FosB": "Fosfomycin resistance",
    # Lincosamide resistance
    "lnuA": "Lincosamide resistance",
    "LNUA": "Lincosamide resistance",
    "vgaALC": "Lincosamide resistance",
    # Two-component systems
    "arlS": "Two-component systems",
    "arlR": "Two-component systems",
    # Global regulators
    "mgrA": "Global regulators",
    # Others
    "fusC": "Fusidic acid",
      "mupA": "Mupirocin resistance",
      "qacF": "Disinfectant resistance",
    "smvA/emrB": "Disinfectant resistance",
   
    # Efflux pumps
    "acrD": "Efflux pumps",
    "acrF/envD": "Efflux pumps",
    "emrB": "Efflux pumps",
    "emrD": "Efflux pumps",
    "mdtA": "Efflux pumps",
    "mdtB": "Efflux pumps",
    "mdtC": "Efflux pumps",
    "mdtF/yhiV": "Efflux pumps",
    "mdtG/yceE": "Efflux pumps",
    "mdtK/ydhE": "Efflux pumps",
    "lmrS": "Efflux pumps",
    "yddg/emrE": "Efflux pumps",
    "ydeP": "Efflux pumps",
    "yodD": "Efflux pumps",

    # Two-component systems
    "baeR": "Two-component systems",
    "baeS": "Two-component systems",
    "cpxA": "Two-component systems",
    "cpxR": "Two-component systems",
    "pmrG": "Two-component systems",
    "phoB": "Two-component systems",
    "zraS/hydG": "Two-component systems",
    "evgS": "Two-component systems",

    # Metal resistance
    "cueP": "Metal resistance",
    "cueR/ybbI": "Metal resistance",
    "golT": "Metal resistance",
  
       "zraP": "Metal resistance",

    # Virulence factors
    "bhsA/ycfR/comC": "Virulence factors",
    "norA": "Virulence factors",
    "modA": "Virulence factors",
    "modB": "Virulence factors",
    "nikR": "Virulence factors",
       
    # Oxidative stress resistance
    "oxyRkp": "Oxidative stress resistance",
    "soxR": "Oxidative stress resistance",
    "soxS": "Oxidative stress resistance",

    # Global regulators
   "iclR": "Global regulators",
   "rpoS": "Global regulators",
   "EVGS": "Virulence factors", 
    "rcnR/yohL": "Global regulators",
    # Virulence factors
    "bhsA/ycfR/comC": "Virulence factors",
    "norA": "Virulence factors",
    "modA": "Virulence factors", 
    "modB": "Virulence factors",
    "nikR": "Virulence factors",
    "CMY": "Beta-lactam resistance",
    "BLAEC": "Beta-lactam resistance",  # Beta-lactamase de E. coli
   "OXA": "Beta-lactam resistance",    # Beta-lactamase classe OXA
   "PBP4B": "Beta-lactam resistance",  # Proteína de ligação à penicilina
   
   # Aminoglycoside resistance
   "AAC6-PRIME": "Aminoglycoside resistance",  # Acetiltransferase
   
   # Efflux pumps
   "QACG": "Efflux pumps",             # Resistência a compostos quaternários
   "MDTM": "Efflux pumps",             # Sistema de efluxo multidroga
   
   # Streptothricin resistance
   "SAT": "Streptothricin resistance",
    # Beta-lactam resistance
    "BLE": "Beta-lactam resistance",       # Resistência a bleomicina (bleomycin resistance)
    
    # Aminoglycoside resistance
    "GADW": "Aminoglycoside resistance", # Gene associado a resistência a aminoglicosídeos
    
    # Efflux pumps
    "MDTP": "Efflux pumps",               # Sistema de efluxo multidroga
    "QACEDELTA1": "Efflux pumps",         # Bombas de efluxo de compostos quaternários
    
    # Folate pathway resistance
    "DFRA": "Folate pathway resistance",  # Resistência a trimetoprima (dihydropteroate reductase)
    
    # Mupirocin resistance
    "MVRC": "Mupirocin resistance" ,
   
   # Aminoglycoside resistance
   "ANT2-DPRIME": "Aminoglycoside resistance",
   "ANT3-DPRIME": "Aminoglycoside resistance",
    
    # Oxidative stress resistance
    "oxyRkp": "Oxidative stress resistance",
    "soxR": "Oxidative stress resistance",
    "soxS": "Oxidative stress resistance",
    "RLMH": "Oxidative stress resistance",  # Corrigido: agora no grupo certo!
 "CATA": "Oxidative stress resistance",

    # Quorum sensing
    "perR": "Quorum sensing",
   
    # Copper resistance
    "BAC0725|copA|sp|Q59385|COPA_ECOLI": "Copper resistance",
    "cuiD": "Copper resistance",
    "cutA": "Copper resistance",
    "cutE/lnt": "Copper resistance",

    # Virulence factorss
    "dsbA": "Virulence factors",
    "dsbB": "Virulence factors",
    "dsbC": "Virulence factors",

    # Iron transport and resistance
    "fieF/yiip": "Iron and zinc efflux",
    "mntH/yfeP": "Manganese and iron transport",
    "mntP/yebN": "Manganese efflux",
    "mntR": "Manganese transport regulator",
    "yieF": "Iron transport",
    "yqjH": "Iron transport",

    # Phosphate transport
    "pstA": "Phosphate transport",
    "pstB": "Phosphate transport",
    "pstC": "Phosphate transport",
    "pstS": "Phosphate transport",
    "pitA": "Phosphate transport",

    # Resistance to other substances
      "terW": "Tellurite resistance",
       "silA": "Silver resistance",

    # Other transport systems
    "glpF": "Glycerol transport",
    "actP/yjcG": "Formate transport",
    "zupT/ygiE": "Zinc transport",
    "zur/yjbK": "Zinc transport",

    # Antibiotic resistance
    "gesA": "Beta-lactam resistance",
    "gesB": "Beta-lactam resistance",
    "gesC": "Beta-lactam resistance",
    "GESB": "Beta-lactam resistance",
    "GESC": "Beta-lactam resistance",

    # Miscellaneous
    "fabI": "Fatty acid biosynthesis",
    "opmD/nmpC": "Outer membrane porin",
    "ostA/lptD": "Outer membrane biogenesis",
    "hns": "Efflux pumps",
    "acrB": "Efflux pumps",
    "acrD": "Efflux pumps",
    "acrD/yffA": "Efflux pumps",
    "acrF/envD": "Efflux pumps",
    "emrB": "Efflux pumps",
    "emrD": "Efflux pumps",
    "emrE/mvrC": "Efflux pumps",
    "mdfA/cmr": "Efflux pumps",
    "mdtB/yegN": "Efflux pumps",
    "mdtC/yegO": "Efflux pumps",
    "mdtF/yhiV": "Efflux pumps",
    "mdtG/yceE": "Efflux pumps",
    "mdtK/ydhE": "Efflux pumps",
    "mdtM/yjiO": "Efflux pumps",
    "mdtN/yjcR": "Efflux pumps",
    "oqxB": "Efflux pumps",
    "tolC": "Efflux pumps",

    # Metal resistance (Resistência a metais pesados)
    "arsA": "Metal resistance (Arsenic)",
    "arsB": "Metal resistance (Arsenic)",
    "arsC": "Metal resistance (Arsenic)",
    "arsD": "Metal resistance (Arsenic)",
    "arsR": "Metal resistance (Arsenic)",
    "cueO": "Metal resistance (Copper)",
    "zinT/yodA": "Metal resistance (Zinc)",

    # Disinfectant resistance (Resistência a desinfetantes)
    "qacE": "Disinfectant resistance",
    "qacEdelta1": "Disinfectant resistance",
      

    # Global regulators and virulence factors (Reguladores globais e virulência)
    "rpoS": "Global regulators",
    "iclR": "Global regulators",
    "comR/ycfQ": "Global regulators",
    "ncrA": "Virulence factors",
    "ncrY": "Virulence factors",
    "ybtP": "Virulence factors",
    "ybtQ": "Virulence factors",
    "ymgB/ariR": "Biofilm regulators",
     # AcrAB-TolC Efflux pumps System
    "acrA": "Efflux pumps",
    "acrB": "Efflux pumps",
    "acrD": "Efflux pumps",
    "acrE/envC": "Efflux pumps",
    "acrF/envD": "Efflux pumps",
    "acrR/ybaH": "Efflux pumps",

   
    # Multidrug Transporters
    "mdfA/cmr": "Efflux pumps",
    "mdtA": "Efflux pumps",
    "mdtB/yegN": "Efflux pumps",
    "mdtC/yegO": "Efflux pumps",
    "mdtE/yhiU": "Efflux pumps",
    "mdtF/yhiV": "Efflux pumps",
    "mdtG/yceE": "Efflux pumps",
    "mdtK/ydhE": "Efflux pumps",
    "mdtM/yjiO": "Efflux pumps",
    "mdtN/yjcR": "Efflux pumps",

    # Heavy Metal Resistance
    "merA": "Mercury resistance",
    "merR": "Mercury resistance",
    "merT": "Mercury resistance",
    "pcoA": "Copper resistance",
    "pcoB": "Copper resistance",
    "pcoC": "Copper resistance",
    "pcoD": "Copper resistance",
    "pcoE": "Copper resistance",
    "pcoR": "Copper resistance",
    "pcoS": "Copper resistance",
    "cusA/ybdE": "Copper resistance",
    "cusB": "Copper resistance",
    "cusC/ylcB": "Copper resistance",
    "cusF/cusX": "Copper resistance",
    "cusR/ylcA": "Copper resistance",
    "cusS": "Copper resistance",

       # Arsenic and Heavy Metal Regulators
      "phoR": "Phosphate regulation",

    # Regulators and Stress Response
    "robA": "Stress response regulators",
   
    # Additional Transporters
    "sitA": "Iron transport",
    "sitB": "Iron transport",
    "sitC": "Iron transport",
    "sitD": "Iron transport",
    "znuA/yebL": "Zinc transport",
    "znuB/yebI": "Zinc transport",
    "znuC/yebM": "Zinc transport",
    "zntA/yhhO": "Zinc transport",
    "zntR/yhdM": "Zinc transport",

    # Multidrug Regulators
    "emrR": "Efflux pumps",
    "marR": "Efflux pumps",

    # Others
    "silA": "Silver resistance",
    "silB": "Silver resistance",
    "silC": "Silver resistance",
    "silE": "Silver resistance",
    "silP": "Silver resistance",
    "silS": "Silver resistance",
    "tehA": "Tellurite resistance",
    "tehB": "Tellurite resistance",
    "terA": "Tellurite resistance",
    "terC": "Tellurite resistance",
    "terE": "Tellurite resistance",
    "terZ": "Tellurite resistance",
    "tolC": "Efflux pumps",
    # Phosphate transport (Transporte de fosfato)
    "pstC": "Phosphate transport",
    "sugE": "Phosphate transport",

    # Oxidative stress resistance (Resistência ao estresse oxidativo)
    "oxyRkp": "Oxidative stress resistance",
    "soxR": "Oxidative stress resistance",
    "sodB": "Oxidative stress resistance",

       # Zinc transport and regulation (Transporte e regulação de zinco)
    "zupT/ygiE": "Zinc transport",
    "zraR/hydH": "Zinc transport",

    # Metal transporters and virulence factors (Transportadores de metais e virulência)
    "fecD": "Iron transport",
    "cutC": "Copper transport",
    "cutF/nlpE": "Copper transport",
    "modC": "Molybdate transport",
    "modE": "Molybdate transport regulators",

    # Virulence factorss (Isomerases de dissulfeto proteico)
    "dsbB": "Virulence factors",
    "dsbC": "Virulence factors",

    # Acid resistance (Resistência a ácidos)
    "gadA": "Acid resistance",
    "gadB": "Acid resistance",
    "gadC/xasA": "Acid resistance",

    # Biofilm regulators (Reguladores de biofilme)
    "ymgB/ariR": "Biofilm regulators",

    # Functional but less documented genes
    "kdeA": "Efflux pumps",
    "kexD": "Efflux pumps",
    "kmrA": "Efflux pumps",
    "kpnE": "Outer membrane protein",
    "kpnO": "Outer membrane porin",

    # Fatty acid biosynthesis (Biossíntese de ácidos graxos)
    "fabI": "Fatty acid biosynthesis",

    # Heat shock proteins (Proteínas de choque térmico)
    "ibpA": "Heat shock proteins",
       # Beta-lactâmicos (Beta-lactam resistance)
    "blaZ": "Beta-lactam resistance",
    "mecA": "Methicillin resistance",
    "mecC": "Methicillin resistance",
    "blaOXA": "Beta-lactam resistance",
    "blaTEM": "Beta-lactam resistance",
    "blaSHV": "Beta-lactam resistance",
    "blaCTX-M": "Extended-spectrum Beta-lactamase (ESBL)",
    "blaNDM": "Beta-lactam resistance",
    "blaKPC": "Carbapenem resistance",
    "blaVIM": "Beta-lactam resistance",
    "blaIMP": "Beta-lactam resistance",
    "blaGES": "Carbapenem resistance",
    # Aminoglycoside resistance (Resistência a aminoglicosídeos)
   "AAC(6')-Iaa": "Aminoglycoside resistance",
   "AAC(6')-Iy": "Aminoglycoside resistance",
   "AAC(3)-Id": "Aminoglycoside resistance",
   "AAC(3)-IId": "Aminoglycoside resistance",
   "AAC(3)-IIe": "Aminoglycoside resistance",
   "ANT(2'')-Ia": "Aminoglycoside resistance",
   "APH(3')-Ia": "Aminoglycoside resistance",
   "APH(4)-Ia": "Aminoglycoside resistance",

   # Beta-lactam resistance (Resistência a beta-lactâmicos)
   "TEM-150": "Beta-lactam resistance",
   "CTX-M-9": "Beta-lactam resistance",
   "SHV-134": "Beta-lactam resistance",
   "CARB-3": "Beta-lactam resistance",
   "CMY-59": "Beta-lactam resistance",
   "catII_from_Escherichia_coli_K-12": "Beta-lactam resistance",

   # Tetracycline resistance (Resistência a tetraciclinas)
   "tet(G)": "Tetracycline resistance",

   # Macrolide resistance (Resistência a macrolídeos)
   "mphE": "Macrolide resistance",
   "mefB": "Macrolide resistance",

   # Sulfonamide resistance (Resistência a sulfonamidas)
   "dfrA12": "Sulfonamide resistance",
   "dfrA14": "Sulfonamide resistance",
   "dfrA15": "Sulfonamide resistance",
   "dfrA16": "Sulfonamide resistance",
   # Efflux Pumps and Transport Systems (Bombas de efluxo e sistemas de transporte)
    "bcr": "Efflux pumps",
    "cmeB": "Efflux pumps",
    "emmdR": "Efflux pumps",

    # Metal Resistance (Resistência a metais)
    "merB": "Metal resistance",

    # Magnesium Transport (Transporte de magnésio)
    "mgtA": "Magnesium transport",
    # Nitrofurantoin Resistance (Resistência a nitrofurantoína)
    "nfsA": "Nitrofurantoin resistance",  # Enzima relacionada à redução de nitrofurantoína.

    # Fimbriae and Adhesion (Fímbria e adesão)
    "faeC": "Fimbriae formation and adhesion",  # Parte do operon de fímbria Fae (pilus de adesão em E. coli).
    "faeD": "Fimbriae formation and adhesion",  # Associado à montagem de fímbrias em patógenos.
    "faeE": "Fimbriae formation and adhesion",  # Envolvido no transporte e montagem de fímbrias.

    # Virulence factors (Fator de virulência)
    "aslA": "Virulence factors",  # Gene associado à patogenicidade em algumas cepas bacterianas.
    "mgtD": "Magnesium transport",  # Caso esteja relacionado a "corD"

    # Copper Transport (Transporte de cobre)
    "corA": "Copper transport",
    "corB": "Copper transport",
    "corC": "Copper transport",
    "corD": "Copper transport",

    # Oxidative Stress Response (Resposta ao estresse oxidativo)
    "sodA": "Oxidative stress resistance",

    # Folate Pathway Resistance (Resistência no caminho do folato)
    "dfrA1_8": "Folate pathway resistance",

    # Virulence factorss (Fatores de virulência)
    "cdtB": "Virulence factors",
    "sak": "Virulence factors",

    # Adhesion and Biofilm Formation (Aderência e formação de biofilme)
    "ychH": "Adhesion and biofilm formation",

   # Quinolo resistance resistance (Resistência a quinolonas)
   "QnrA1": "Quinolone resistance", 
 # Aminoglycoside resistance (Resistência a aminoglicosídeos)
    "AAC(3)-IV": "Aminoglycoside resistance",
    "AAC(3)-VIa": "Aminoglycoside resistance",
    "aadA2_1": "Aminoglycoside resistance",
    "aadA7_1": "Aminoglycoside resistance",
    "aadA5_1": "Aminoglycoside resistance",
    "str_1": "Aminoglycoside resistance",
    "aph(3')-Ia_1": "Aminoglycoside resistance",
    "aph(3')-Ia_9": "Aminoglycoside resistance",

    # Beta-lactam resistance (Resistência a beta-lactâmicos)
    "blaCMY-2_1": "Beta-lactam resistance",
    "TEM-12": "Beta-lactam resistance",

    # Iron acquisition systems (Sistemas de aquisição de ferro)
    "fyuA": "Iron acquisition system",
    "ybtE": "Iron acquisition system",
    "ybtT": "Iron acquisition system",
    "ybtU": "Iron acquisition system",
    "irp1": "Iron acquisition system",
    "irp2": "Iron acquisition system",
    "ybtA": "Iron acquisition system",
    "ybtX": "Iron acquisition system",
    "ybtS": "Iron acquisition system",

    # Adhesion and biofilm formation (Aderência e formação de biofilme)
    "fimB": "Adhesion and biofilm formation",
    "fimE": "Adhesion and biofilm formation",
    "fimA": "Adhesion and biofilm formation",
    "fimG": "Adhesion and biofilm formation",
    "fdeC": "Adhesion and biofilm formation",
    "ykgK/ecpR": "Adhesion and biofilm formation",
    "yagZ/ecpA": "Adhesion and biofilm formation",
    "yagY/ecpB": "Adhesion and biofilm formation",
    "yagX/ecpC": "Adhesion and biofilm formation",
    "yagW/ecpD": "Adhesion and biofilm formation",
    "yagV/ecpE": "Adhesion and biofilm formation",

    # Type VI secretion systems (Sistemas de secreção tipo VI)
    "gspC": "Type VI secretion system",
    "gspD": "Type VI secretion system",
    "gspK": "Type VI secretion system",
    "gspL": "Type VI secretion system",
    "gspM": "Type VI secretion system",

    # Iron siderophore systems (Sistemas sideróforos de ferro)
    "iroB": "Iron siderophore system",
    "iroC": "Iron siderophore system",
    "iroD": "Iron siderophore system",
    "iroE": "Iron siderophore system",
    "iroN": "Iron siderophore system",

    # Effector proteins (Proteínas efetoras)
    "sopA": "Effector protein",
    "sseL": "Effector protein",
    "sspH1": "Effector protein",
    "espL4": "Effector protein",
    "espX4": "Effector protein",
    "espX5": "Effector protein",
    "espL1": "Effector protein",
    "espX1": "Effector protein",
    "espY1": "Effector protein",
    "espR1": "Effector protein",
    "esaB": "Effector protein",
    "esaC": "Effector protein",
    "esxA": "Effector protein",
# Iron Transport System (Sistema de transporte de ferro)
    "fepA": "Iron transport system",  # Receptor de sideróforo (ferric enterobactin receptor) que facilita a captação de ferro em bactérias.
     

   # Multidrug efflux (Efluxo multidroga)
   
   # Colistin resistance (Resistência a colistina)
   "MCR-9": "Colistin resistance",

   # Outras categorias para genes únicos
   "aadA2": "Aminoglycoside resistance",
   "aadA3": "Aminoglycoside resistance",
   "aadA7": "Aminoglycoside resistance",
   "SAT-1": "Aminoglycoside resistance",
# Type III Secretion System (T3SS) Effectors
    "steA": "Type III Secretion System (T3SS)",
    "sifB": "Type III Secretion System (T3SS)",
    "steB": "Type III Secretion System (T3SS)",
    "sseJ": "Type III Secretion System (T3SS)",
    "sseK2": "Type III Secretion System (T3SS)",
    "sopE2": "Type III Secretion System (T3SS)",
    "sopD2": "Type III Secretion System (T3SS)",
    "sseK1": "Type III Secretion System (T3SS)",
    "steC": "Type III Secretion System (T3SS)",
    "pipB": "Type III Secretion System (T3SS)",
    "sopB/sigD": "Type III Secretion System (T3SS)",
    "sopD": "Type III Secretion System (T3SS)",
    "gogB": "Type III Secretion System (T3SS)",
    "slrP": "Type III Secretion System (T3SS)",
    "sifA": "Type III Secretion System (T3SS)",
    "spiC/ssaB": "Type III Secretion System (T3SS)",

  

    # Fimbrial Adhesins
    "fimF": "Fimbrial Adhesin",
    "fimH": "Fimbrial Adhesin",
    "fimD": "Fimbrial Assembly",
    "fimC": "Fimbrial Assembly",
    "fimI": "Fimbrial Assembly",

    # Stress Resistance / Virulence
       "sodCI": "Oxidative stress resistance",
    # Beta-lactam resistance
    "blaCTX-M-9_1": "Extended-spectrum Beta-lactamase (ESBL)",
   "blaCTX-M-15_1": "Extended-spectrum Beta-lactamase (ESBL)",
   "blaSHV-12_1": "Extended-spectrum Beta-lactamase (ESBL)",
   
   "blaCARB-2_1": "Carbapenem resistance",
 

    # Type III Secretion System (T3SS) Structural and Regulators
    "invH": "Type III Secretion System (T3SS)",
    "invF": "Type III Secretion System (T3SS)",
    "invG": "Type III Secretion System (T3SS)",
    "invE": "Type III Secretion System (T3SS)",
    "invA": "Type III Secretion System (T3SS)",
    "invB": "Type III Secretion System (T3SS)",
    "invC": "Type III Secretion System (T3SS)",
    "invI": "Type III Secretion System (T3SS)",
    "invJ": "Type III Secretion System (T3SS)",
    "spaO": "Type III Secretion System (T3SS)",
    "spaP": "Type III Secretion System (T3SS)",
    "spaQ": "Type III Secretion System (T3SS)",
    "spaR": "Type III Secretion System (T3SS)",
    "spaS": "Type III Secretion System (T3SS)",
    "sicA": "Type III Secretion System (T3SS)",
    "sipB/sspB": "Type III Secretion System (T3SS)",
    "sipC/sspC": "Type III Secretion System (T3SS)",
    "sipD": "Type III Secretion System (T3SS)",
    "sipA/sspA": "Type III Secretion System (T3SS)",
    "sicP": "Type III Secretion System (T3SS)",
    "sptP": "Type III Secretion System (T3SS)",
    "prgH": "Type III Secretion System (T3SS)",
    "prgI": "Type III Secretion System (T3SS)",
    "prgJ": "Type III Secretion System (T3SS)",
    "prgK": "Type III Secretion System (T3SS)",
    "orgA": "Type III Secretion System (T3SS)",
    "orgB": "Type III Secretion System (T3SS)",
    "orgC": "Type III Secretion System (T3SS)",

    # Fimbrial Adhesins
    "fimF": "Fimbrial Adhesin",
    "fimH": "Fimbrial Adhesin",
    "fimD": "Fimbrial Assembly",
    "fimC": "Fimbrial Assembly",
    "fimI": "Fimbrial Assembly",

    # Stress Resistance / Virulence
    "mgtC": "Stress Resistance ",
    "mgtB": "Stress Resistance ",
    "sodCI": "Oxidative stress resistance",

    # Siderophores and Iron Acquisition
    "fepC": "Siderophore Transport",
    "fepG": "Siderophore Transport",
    "entB": "Siderophore Biosynthesis",
    "entA": "Siderophore Biosynthesis",

    # Virulence factorss
    "misL": "Virulence factors",
    "lpfA": "Virulence factors",
    "lpfB": "Virulence factors",
    "lpfC": "Virulence factors",
    "lpfD": "Virulence factors",
    "lpfE": "Virulence factors",
    "ratB": "Virulence factors",
    "sinH": "Virulence factors",
    "sseI/srfH": "Virulence factors",
    "ompA": "Virulence factors",
      "avrA": "Virulence factors",
    "mig-14": "Virulence factors",
    "pipB2": "Virulence factors",
    "sspH2": "Virulence factors",
    "rck": "Virulence factors",

    # Pathogenicity Island (SPI-2)
    "ssaU": "Pathogenicity Island (SPI-2)",
    "ssaT": "Pathogenicity Island (SPI-2)",
    "ssaS": "Pathogenicity Island (SPI-2)",
    "ssaR": "Pathogenicity Island (SPI-2)",
    "ssaQ": "Pathogenicity Island (SPI-2)",
    "ssaP": "Pathogenicity Island (SPI-2)",
    "ssaO": "Pathogenicity Island (SPI-2)",
    "ssaN": "Pathogenicity Island (SPI-2)",
    "ssaV": "Pathogenicity Island (SPI-2)",
    "ssaM": "Pathogenicity Island (SPI-2)",
    "ssaL": "Pathogenicity Island (SPI-2)",
    "ssaK": "Pathogenicity Island (SPI-2)",
    "ssaJ": "Pathogenicity Island (SPI-2)",
    "ssaI": "Pathogenicity Island (SPI-2)",
    "ssaH": "Pathogenicity Island (SPI-2)",
    "ssaG": "Pathogenicity Island (SPI-2)",
    "sseG": "Pathogenicity Island (SPI-2)",
    "sseF": "Pathogenicity Island (SPI-2)",
    "sscB": "Pathogenicity Island (SPI-2)",
    "sseE": "Pathogenicity Island (SPI-2)",
    "sseD": "Pathogenicity Island (SPI-2)",
    "sseC": "Pathogenicity Island (SPI-2)",
    "sscA": "Pathogenicity Island (SPI-2)",
    "sseB": "Pathogenicity Island (SPI-2)",
    "sseA": "Pathogenicity Island (SPI-2)",
    "ssaE": "Pathogenicity Island (SPI-2)",
    "ssaD": "Pathogenicity Island (SPI-2)",
    "ssaC": "Pathogenicity Island (SPI-2)",

    # Curli assembly System
    "csgC": "Curli assembly",
    "csgA": "Curli assembly",
    "csgB": "Curli assembly",
    "csgD": "Curli assembly",
    "csgE": "Curli assembly",
    "csgF": "Curli assembly",
    "csgG": "Curli assembly",

    # Plasmid-Encoded Factors
    "pefD": "Plasmid-Encoded Factor",
    "pefC": "Plasmid-Encoded Factor",
    "pefA": "Plasmid-Encoded Factor",
    "pefB": "Plasmid-Encoded Factor",
    "spvC": "Plasmid-Encoded Factor",
    "spvB": "Plasmid-Encoded Factor",
    "spvR": "Plasmid-Encoded Factor",
    # Type III Secretion System (T3SS) Effectors
   
 
    # Fimbrial Adhesins
    "fimF": "Fimbrial Adhesin",
    "fimH": "Fimbrial Adhesin",
    "fimD": "Fimbrial Assembly",
    "fimC": "Fimbrial Assembly",
    "fimI": "Fimbrial Assembly",

    # Stress Resistance / Virulence
      "sodCI": "Oxidative stress resistance",
    "grvA": "Oxidative stress resistance",

    # Siderophores and Iron Acquisition
    "fepC": "Siderophore Transport",
    "fepG": "Siderophore Transport",
    "entB": "Siderophore Biosynthesis",
    "entA": "Siderophore Biosynthesis",
    "entE": "Siderophore Biosynthesis",
    "entD": "Siderophore Biosynthesis",
    "fes": "Siderophore Biosynthesis",
    "entS": "Siderophore Biosynthesis",
    "fepD": "Siderophore Transport",
    "fepB": "Siderophore Transport",
    "entC": "Siderophore Biosynthesis",
    "entF": "Siderophore Biosynthesis",

    # Virulence factorss
    "misL": "Virulence factors",
    "lpfA": "Virulence factors",
    "lpfB": "Virulence factors",
    "lpfC": "Virulence factors",
    "lpfD": "Virulence factors",
    "lpfE": "Virulence factors",
    "ratB": "Virulence factors",
    "sinH": "Virulence factors",
    "sseI/srfH": "Virulence factors",
     "avrA": "Virulence factors",
    "mig-14": "Virulence factors",
    "pipB2": "Virulence factors",
    "sspH2": "Virulence factors",
    "rck": "Virulence factors",
    "shdA": "Virulence factors",

   
      # Beta-lactam resistance
    "blaTEM-1B_1": "Beta-lactam resistance",
    "blaTEM-1D_1": "Beta-lactam resistance",
    "blaTEM-1A_1": "Beta-lactam resistance",
    "blaCTX-M-9_1": "Extended-spectrum Beta-lactamase (ESBL)",
    "blaCTX-M-15_1": "Extended-spectrum Beta-lactamase (ESBL)",
    "blaSHV-12_1": "Extended-spectrum Beta-lactamase (ESBL)",
    "blaOXA-1_1": "Carbapenem resistance",
    "blaCARB-2_1": "Carbapenem resistance",
       # Beta-lactam resistance
    "blaTEM-1B_1": "Beta-lactam resistance",
    "blaTEM-1D_1": "Beta-lactam resistance",
    "blaTEM-1A_1": "Beta-lactam resistance",
    "blaCTX-M-9_1": "Extended-spectrum Beta-lactamase (ESBL)",
    "blaCTX-M-15_1": "Extended-spectrum Beta-lactamase (ESBL)",
    "blaSHV-12_1": "Extended-spectrum Beta-lactamase (ESBL)",
    "blaCARB-2_1": "Carbapenem resistance",

    # Aminoglycoside Modifying Enzymes (AMEs)
    "aac(6')-Iaa_1": "Aminoglycoside resistance",
    "aac(3)-Id_1": "Aminoglycoside resistance",
    "aac(3)-IIa_1": "Aminoglycoside resistance",
    "aac(3)-VIa_2": "Aminoglycoside resistance",
    "aac(3)-IVa_1": "Aminoglycoside resistance",
    "aac(3)-IId_1": "Aminoglycoside resistance",
    "ant(3'')-Ia_1": "Aminoglycoside resistance",
    "ant(2'')-Ia_1": "Aminoglycoside resistance",
    "aph(3'')-Ib_1": "Aminoglycoside resistance",
    "aph(3'')-Ib_2": "Aminoglycoside resistance",
    "aph(3'')-Ib_5": "Aminoglycoside resistance",
    "aph(3')-Ia_3": "Aminoglycoside resistance",
    "aph(3')-Ia_7": "Aminoglycoside resistance",
    "aph(6)-Id_1": "Aminoglycoside resistance",
    "aph(4)-Ia_1": "Aminoglycoside resistance",

    # Tetracycline Resistance
    "tet(A)_6": "Tetracycline resistance",
    "tet(B)_2": "Tetracycline resistance",
    "tet(C)_3": "Tetracycline resistance",
    "tet(G)_2": "Tetracycline resistance",
    "tet(M)_4": "Tetracycline resistance",
    "tet(M)_8": "Tetracycline resistance",
    "tet(K)_1": "Tetracycline resistance",

    # Sulfonamide Resistance
    "sul1_5": "Sulfonamide resistance",
    "sul2_3": "Sulfonamide resistance",
    "sul3_2": "Sulfonamide resistance",

    # FluoroQuinolo resistance Resistance
    "qnrB19_1": "Quinolone resistance",
    "qnrA1_1": "Quinolone resistance",
    "qnrS1_1": "Quinolone resistance",
    "qnrD1_1": "Quinolone resistance",

    # Fosfomycin resistance
    "fosA7_1": "Fosfomycin resistance",

    # Macrolide resistance
    "msr(E)_1": "Macrolide resistance",
    "mph(E)_1": "Macrolide resistance",
    "mph(B)_1": "Macrolide resistance",
    "mef(B)_1": "Macrolide resistance",

    # Phenicol resistance
    
    # Trimethoprim resistance
    "dfrA1_10": "Trimethoprim resistance",
    "dfrA12_8": "Trimethoprim resistance",
    "dfrA16_2": "Trimethoprim resistance",
    "dfrA17_1": "Trimethoprim resistance",
    "dfrA15_2": "Trimethoprim resistance",
    "dfrA14_5": "Trimethoprim resistance",

    # Colistin Resistance
   

    # Other Resistance Genes
    "armA_1": "16S RRNA Methyltransferase",
    "lnu(F)_1": "Lincosamide resistance",
    "lnu(A)_1": "Lincosamide resistance",
    "rck": "Virulence-associated Gene",
    "golS": "Copper/Silver Resistance",
    "mdsA": "Copper/Silver Resistance",
    "mdsB": "Copper/Silver Resistance",
    "mdsC": "Copper/Silver Resistance",
     "sspH2": "Virulence factors",
    "pipB2": "Virulence factors",
    "avrA": "Virulence factors",
    "mig-14": "Virulence factors",
    "csgA": "Curli assembly",
    "csgB": "Curli assembly",
    "csgC": "Curli assembly",

   
 
    # Sulfonamide Resistance
  
    # Fluoroquinolone Resistance
     "qnrA1_1": "Quinolone resistance",
    "qnrS1_1": "Quinolone resistance",
    "qnrD1_1": "Quinolone resistance",

    # Fosfomycin resistance
    "fosA7_1": "Fosfomycin resistance",

    
    # Phenicol resistance
    "floR_2": "Phenicol resistance",
    "catA1_1": "Phenicol resistance",
    "catA2_1": "Phenicol resistance",
    "cmlA1_1": "Phenicol resistance",

    # Trimethoprim resistance
    "dfrA1_10": "Trimethoprim resistance",
    "dfrA12_8": "Trimethoprim resistance",
    "dfrA16_2": "Trimethoprim resistance",
    "dfrA17_1": "Trimethoprim resistance",
    "dfrA15_2": "Trimethoprim resistance",
    "dfrA14_5": "Trimethoprim resistance",

    # Colistin Resistance
    "mcr-1.1_1": "Colistin resistance",
    "mcr-9_1": "Colistin resistance",

    # Other Resistance Genes
       "lnu(F)_1": "Lincosamide resistance",
    "lnu(A)_1": "Lincosamide resistance",
    "rck": "Virulence-associated Gene",
    "golS": "Copper/Silver Resistance",
    "mdsA": "Copper/Silver Resistance",
    "mdsB": "Copper/Silver Resistance",
    "mdsC": "Copper/Silver Resistance",
     "sspH2": "Virulence factors",
    "pipB2": "Virulence factors",
    "avrA": "Virulence factors",
    "mig-14": "Virulence factors",
    "csgA": "Curli assembly",
    "csgB": "Curli assembly",
    "csgC": "Curli assembly",



 
   # Sulfonamide Resistance

   # Fluoroquinolone Resistance
       "qnrA1_1": "Quinolone resistance",
   "qnrS1_1": "Quinolone resistance",
   "qnrD1_1": "Quinolone resistance",

   # Fosfomycin resistance
   "fosA7_1": "Fosfomycin resistance",

   # Macrolide resistance
   "msr(E)_1": "Macrolide resistance",
   "mph(E)_1": "Macrolide resistance",
   "mph(B)_1": "Macrolide resistance",
   "mef(B)_1": "Macrolide resistance",

   # Phenicol resistance
   "floR_2": "Phenicol resistance",
   "catA1_1": "Phenicol resistance",
   "catA2_1": "Phenicol resistance",
   "cmlA1_1": "Phenicol resistance",

   # Trimethoprim resistance
   "dfrA1_10": "Trimethoprim resistance",
   "dfrA12_8": "Trimethoprim resistance",
   "dfrA16_2": "Trimethoprim resistance",
   "dfrA17_1": "Trimethoprim resistance",
   "dfrA15_2": "Trimethoprim resistance",
   "dfrA14_5": "Trimethoprim resistance",

   
   # Other Resistance Genes
   "mdf(A)_1": "Efflux pumps",
      "lnu(F)_1": "Lincosamide resistance",
   "lnu(A)_1": "Lincosamide resistance",
   "rck": "Virulence factors",
   "golS": "Copper/Silver Resistance",
   "mdsA": "Copper/Silver Resistance",
   "mdsB": "Copper/Silver Resistance",
   "mdsC": "Copper/Silver Resistance",
     "sspH2": "Virulence factors",
   "pipB2": "Virulence factors",
   "avrA": "Virulence factors",
   "mig-14": "Virulence factors",
   "csgA": "Curli assembly",
   "csgB": "Curli assembly",
   "csgC": "Curli assembly",

    # Siderophores and Iron Acquisition
    "fepC": "Siderophore Transport",
    "fepG": "Siderophore Transport",
    "entB": "Siderophore Biosynthesis",
    "entA": "Siderophore Biosynthesis",

    # Virulence factorss
    "misL": "Virulence factors",
    "lpfA": "Virulence factors",
    "lpfB": "Virulence factors",
    "lpfC": "Virulence factors",
    "lpfD": "Virulence factors",
    "lpfE": "Virulence factors",
    "ratB": "Virulence factors",
    "sinH": "Virulence factors",
    "sseI/srfH": "Virulence factors",
    "avrA": "Virulence factors",
    "mig-14": "Virulence factors",
    "pipB2": "Virulence factors",
    "sspH2": "Virulence factors",
    # Macrolídeos (Macrolide resistance)
    "ermA": "Aminoglycoside resistance",
    "ermB": "Aminoglycoside resistance",
    "ermC": "Aminoglycoside resistance",
    "msrA": "Macrolide resistance",
    "msrB": "Macrolide resistance",
    "mphA": "Macrolide resistance",
    "mphB": "Macrolide resistance",
    "mefA": "Macrolide resistance",
    "mefE": "Macrolide resistance",

    # Glicopeptídeos (Glycopeptides)
    "vanA": "Glycopeptide resistance",
    "vanB": "Glycopeptide resistance",
    "vanC": "Glycopeptide resistance",
    "vanD": "Glycopeptide resistance",
    "vanE": "Glycopeptide resistance",
    "vanG": "Glycopeptide resistance",
    "vanX": "Glycopeptide resistance",
    "vanY": "Glycopeptide resistance",
    "vanZ": "Glycopeptide resistance",

    # Quinolonas (Quinolones)
    "gyrA": "Quinolone resistance",
    "gyrB": "Quinolone resistance",
    "parC": "Quinolone resistance",
    "parE": "Quinolone resistance",
    "qnrA": "Quinolone resistance",
    "qnrB": "Quinolone resistance",
    "qnrC": "Quinolone resistance",
    "qnrD": "Quinolone resistance",
    "qnrS": "Quinolone resistance",
    "aac(6')-Ib-cr": "Quinolone resistance",
    "QNRA": "Quinolone resistance",


  
  # Specific Subgroups
   "IncHI1B(CIT)_1_pNDM-CIT_JX182975": "Large plasmid",


    # Incompatibility Groups (Inc)
   "IncX4_1__CP002895": "Large plasmid",
    "IncR_1__DQ449578": "Large plasmid",
    "IncX1_1__EU370913": "Large plasmid",
    "IncHI2_1__BX664015": "Large plasmid",
    "IncHI2A_1__BX664015": "Large plasmid",
    "IncFIB(AP001918)_1__AP001918": "Large plasmid",
    "IncFII_1__AY458016": "Large plasmid",
    "IncB/O/K/Z_4__FN868832": "Large plasmid",
    "IncFIC(FII)_1__AP001918": "Large plasmid",
    "IncY_1__K02380": "Large plasmid",
    "IncX3_1__JN247852": "Large plasmid",
    "IncN_1__AY046276": "Large plasmid",
    "IncQ1_1__M28829.1": "Large plasmid",
    "IncL/M(pOXA-48)_1_pOXA-48_JN626286": "Large plasmid",
    "IncP1_3__KX377410": "Large plasmid",
    "IncI1_1_Alpha_AP005147": "Large plasmid",
    "IncI_Gamma_1_AP011954": "Large plasmid",

    # Small plasmids (Colicin-type)
    "ColRNAI_1__DQ298019": "Small plasmid",
    "Col(MG828)_1__NC_008486": "Small plasmid",
    "Col156_1__NC_009781": "Small plasmid",
    "Col440I_1__CP023920.1": "Small plasmid",
    "Col440II_1__CP023921.1": "Small plasmid",
    "Col(KPHS6)_1__NC_016841": "Small plasmid",
    "Col(BS512)_1__NC_010656": "Small plasmid",
    "Col8282_1__DQ995352": "Small plasmid",
    "ColpVC_1__JX133088": "Small plasmid",
    
  # Large plasmids
    "IncHI1B(R27)_1_R27_AF250878": "Large plasmid",
    "IncFIA(HI1)_1_HI1_AF250878": "Large plasmid",
    "IncFIB(K)_1_Kpn3_JN233704": "Large plasmid",
    "IncFIB(pHCM2)_1_pHCM2_AL513384": "Large plasmid",
    # Small plasmids
    "RepA_1_pKPC-CAV1321_CP011611": "Small plasmid",
    "p0111_1__AP010962": "Small plasmid",
    "pENTAS02_1__CP003028": "Small plasmid",
   "IncI2_1_Delta_AP002527": "Large plasmid",
   "Col3M_1__JX514065": "Small plasmid",
   "IncFIB(pECLA)_1_pECLA_CP001919": "Large plasmid",
    # Tetraciclinas (Tetracycline resistance)
    "tetA": "Tetracycline resistance",
    "tetB": "Tetracycline resistance",
    "tetC": "Tetracycline resistance",
    "tetD": "Tetracycline resistance",
    "tetE": "Tetracycline resistance",
    "tetG": "Tetracycline resistance",
    "tetK": "Tetracycline resistance",
    "tetL": "Tetracycline resistance",
    "tetM": "Tetracycline resistance",
    "tetO": "Tetracycline resistance",
    "tetQ": "Tetracycline resistance",
    "tetS": "Tetracycline resistance",
    "tetX": "Tetracycline resistance",
    "TETK": "Tetracycline resistance",
   "TETG": "Tetracycline resistance",

    # Aminoglicosídeos (Aminoglycoside resistance)
    "aac(3')": "Aminoglycoside resistance",
    "aac(6')": "Aminoglycoside resistance",
    "aadA": "Aminoglycoside resistance",
    "ant(2'')": "Aminoglycoside resistance",
    "aph(3')": "Aminoglycoside resistance",
    "aph(4)": "Aminoglycoside resistance",
    "aph(6')": "Aminoglycoside resistance",
    "armA": "Aminoglycoside resistance",
    "rmtA": "Aminoglycoside resistance",
    "rmtB": "Aminoglycoside resistance",
    "rmtC": "Aminoglycoside resistance",

    # Sulfametoxazol/Trimetoprima (Sulfamethoxazole/Trimethoprim)
    "sul1": "Sulfonamide resistance",
    "sul2": "Sulfonamide resistance",
    "sul3": "Sulfonamide resistance",
    "dfrA": "Trimethoprim resistance",
    "dfrB": "Trimethoprim resistance",
    "dfrC": "Trimethoprim resistance",

    # Fenicóis (Phenicols)
    "catA": "Phenicol resistance",
    "catB": "Phenicol resistance",
    "cmlA": "Phenicol resistance",
    "floR": "Phenicol resistance",
    "cmr": "Phenicol resistance",
    "optrA": "Phenicol resistance",

    # Oxazolidinonas (Oxazolidinones)
    "cfr": "Oxazolidinone resistance",
    "poxtA": "Oxazolidinone resistance",

    # Rifamicinas (Rifamycins)
    "rpoB": "Rifamycin resistance",
    "arr": "Rifamycin resistance",

    # Lipopeptídeos (Lipopeptides)
    "mprF": "Lipopeptides resistance",
    "yycG": "Lipopeptides resistance",
    "yycF": "Lipopeptides resistance",

    # Polimixinas (Polymyxins)
    "mcr-1": "Colistin resistance",
    "mcr-2": "Phosphoethanolamine transferase",
    "mcr-3": "Phosphoethanolamine transferase",
    "mcr-4": "Phosphoethanolamine transferase",
    "mcr-5": "Phosphoethanolamine transferase",

    # Fosfomicina (Fosfomycin)
    "fosA": "Fosfomycin resistance",
    "fosB": "Fosfomycin resistance",
    "fosC": "Fosfomycin resistance",
    "fosX": "Fosfomycin resistance",
    "Klebsiella_pneumoniae_KpnG": "Type VI secretion system",
    "Klebsiella_pneumoniae_KpnH": "Type VI secretion system",
    "Klebsiella_pneumoniae_OmpK37": "Outer membrane protein",
    "Klebsiella_pneumoniae_KpnE": "Virulence factors",
    "Klebsiella_pneumoniae_KpnF": "Efflux pumps",
    "SHV-182": "Beta-lactam resistance",
    "QnrB17": "Quinolone resistance",
    "AAC(6')-Ib-cr": "Aminoglycoside resistance",
    "Klebsiella_pneumoniae_acrA": "Efflux pumps",
    "SHV-1": "Beta-lactam resistance",
    "FosA5": "Fosfomycin resistance",
    "SHV-107": "Beta-lactam resistance",
    "dfrA25": "Trimethoprim resistance",
    "AAC(6')-Ib9": "Aminoglycoside resistance",
    "rmtF": "Aminoglycoside resistance",
    "catB2": "",
    "dfrB1": "Trimethoprim resistance",
    "VIM-1": "Carbapenem resistance",
    "SHV-108": "Beta-lactam resistance",
    "arr-2": "Rifamycin resistance",
    "OKP-B-8": "Beta-lactam resistance",
    "OXA-9": "Beta-lactam resistance",
    "TEM-141": "Beta-lactam resistance",
    "KPC-3": "Carbapenem resistance",
    "EreA2": "Macrolide resistance",
    "dfrA19": "Trimethoprim resistance",
    "catB3": "Phenicol resistance",
    "rmtE": "Aminoglycoside resistance",
    "dfrA7": "Trimethoprim resistance",
    "TEM-104": "Beta-lactam resistance",
    "SHV-27": "Beta-lactam resistance",
    "NDM-1": "Carbapenem resistance",
    "APH(3')-VI": "Aminoglycoside resistance",
    "SCO-1": "Beta-lactam resistance",
    "SHV-62": "Beta-lactam resistance",
    "QnrB1": "Quinolone resistance",
    "APH(6)-Ic": "Aminoglycoside resistance",
    "APH(3')-IIa": "Aminoglycoside resistance",
       "LAP-2": "Beta-lactam resistance",
    "oqxA_1": "Efflux pumps",
     # Beta-lactâmicos
    "blaSHV-27_1": "Beta-lactam resistance",
    "blaSHV-187_1": "Beta-lactam resistance",
    "blaSHV-182_1": "Beta-lactam resistance",
    "blaDHA-17_1": "Beta-lactam resistance",
    "blaCTX-M-173_1": "Beta-lactam resistance",
    "blaOXA-48_1": "Beta-lactam resistance",
    "blaOXA-9_1": "Beta-lactam resistance",
    "blaTEM-141_1": "Beta-lactam resistance",
    "blaKPC-3_1": "Beta-lactam resistance",
    "blaSCO-1_1": "Beta-lactam resistance",
   
    # Resistência a quinolonas
    "aac(6')-Ib-cr_1": "Quinolone resistance",
    "qnrB4_1": "Quinolone resistance",
    "qnrB1_1": "Quinolone resistance",

    # Resistência a aminoglicosídeos
    "rmtE_1": "Aminoglycoside resistance",
    "rmtF_1": "Aminoglycoside resistance",
    "ARR-3_4": "Aminoglycoside resistance",
    "APH(3')-VI_1": "Aminoglycoside resistance",
    "aph(6)-Ic_1": "Aminoglycoside resistance",
    "aph(3')-IIa_2": "Aminoglycoside resistance",

        # Resistência a trimetoprim
    "dfrA19_1": "Trimethoprim resistance",
    "dfrB1_1": "Trimethoprim resistance",
    "dfrA7_5": "Trimethoprim resistance",

    # Resistência a tetraciclinas
  
    # Resistência a cloranfenicol
    "catB2_1": "Phenicol resistance",
    "catB3_2": "Phenicol resistance",

    # Resistência a fosfomicina
    "fosA5_1": "Fosfomycin resistance",
    "fosA6_1": "Fosfomycin resistance",
    "fosA_3": "Fosfomycin resistance",

    # Resistência a macrolídeos
    "mph(A)_2": "Macrolide resistance",
    "mefC": "Macrolide resistance",
      # Beta-lactâmicos
    "blaSHV-27_1": "Beta-lactam resistance",
    "blaSHV-187_1": "Beta-lactam resistance",
    "blaSHV-182_1": "Beta-lactam resistance",
    "blaDHA-17_1": "Beta-lactam resistance",
    "blaCTX-M-173_1": "Beta-lactam resistance",
    "blaOXA-48_1": "Beta-lactam resistance",
    "blaOXA-9_1": "Beta-lactam resistance",
    "blaTEM-141_1": "Beta-lactam resistance",
    "blaKPC-3_1": "Beta-lactam resistance",
    "blaSCO-1_1": "Beta-lactam resistance",
      "blaOKP-B-8_1": "Beta-lactam resistance",
    "blaTEM-1C_1": "Beta-lactam resistance",
    "blaSHV-110_1": "Beta-lactam resistance",

    # Resistência a quinolonas
    "aac(6')-Ib-cr_1": "Quinolone resistance",
    "qnrB4_1": "Quinolone resistance",
    "qnrB1_1": "Quinolone resistance",
      "aph(3')-VI_1": "Aminoglycoside resistance",
    # Resistência a aminoglicosídeos
    "rmtE_1": "Aminoglycoside resistance",
    "rmtF_1": "Aminoglycoside resistance",
    "ARR-3_4": "Aminoglycoside resistance",
    "APH(3')-VI_1": "Aminoglycoside resistance",
    "aph(6)-Ic_1": "Aminoglycoside resistance",
    "aph(3')-IIa_2": "Aminoglycoside resistance",

  

    # Resistência a trimetoprim
    "dfrA19_1": "Trimethoprim resistance",
    "dfrB1_1": "Trimethoprim resistance",
    "dfrA7_5": "Trimethoprim resistance",



    # Resistência a cloranfenicol
    "catB2_1": "Phenicol resistance",
    "catB3_2": "Phenicol resistance",

    # Resistência a fosfomicina
    "fosA5_1": "Fosfomycin resistance",
    "fosA6_1": "Fosfomycin resistance",
    "fosA_3": "Fosfomycin resistance",
    "fosA_6": "Fosfomycin resistance",

    # Resistência a macrolídeos
    "mph(A)_2": "Macrolide resistance",
    "mefC": "Macrolide resistance",

    # Resistência a carbapenêmicos (Metalo-Beta-lactam resistance)
    "blaVIM-1_1": "Carbapenem resistance",
    "blaNDM-1_1": "Carbapenem resistance",

    # Resistência a colistina


    # Fatores de virulência
    "iutA": "Virulence factors",
    "iucC": "Virulence factors",
    "iucB": "Virulence factors",
    "iucA": "Virulence factors",
    "iucD": "Virulence factors",
    "cseA": "Virulence factors",
    "hlyA": "Virulence factors",
    "hlyB": "Virulence factors",
    "hlyD": "Virulence factors",
    "sat": "Virulence factors",
    "papX": "Virulence factors",
    "papG": "Virulence factors",
    "papF": "Virulence factors",
    "papE": "Virulence factors",
    "papK": "Virulence factors",
    "papJ": "Virulence factors",
    "papD": "Virulence factors",
    "papC": "Virulence factors",
    "papH": "Virulence factors",
    "papA": "Virulence factors",
    "papB": "Virulence factors",
    "papI": "Virulence factors",
    "chuA": "Virulence factors",
    "chuS": "Virulence factors",
    "chuT": "Virulence factors",
    "chuW": "Virulence factors",
    "chuX": "Virulence factors",
    "chuY": "Virulence factors",
    "chuU": "Virulence factors",
    "chuV": "Virulence factors",
    "kpsD": "Virulence factors",
    "kpsM": "Virulence factors",
    "pic": "Virulence factors",
    "sfaX": "Virulence factors",
    "sfaY": "Virulence factors",
    "focH": "Virulence factors",
    "focG": "Virulence factors",
    "focF": "Virulence factors",
    "focD": "Virulence factors",
    "focC": "Virulence factors",
    "focA": "Virulence factors",
    "sfaD": "Virulence factors",
    "sfaB": "Virulence factors",
    "sfaC": "Virulence factors",
    "vat": "Virulence factors",
    "essA": "Virulence factors",
    "astA": "Virulence factors",

  
  
    
    # Resistência a colistina
    "mcr-4.2_1": "Colistin resistance",

    # Fatores de virulência
    "iutA": "Virulence factors",
    "iucC": "Virulence factors",
    "iucB": "Virulence factors",
    "iucA": "Virulence factors",
    "cseA": "Virulence factors",
    "hlyA": "Virulence factors",
    "hlyB": "Virulence factors",
    "hlyD": "Virulence factors",
    "sat": "Virulence factors",
    "papX": "Virulence factors",
    "chuA": "Virulence factors",
    "chuS": "Virulence factors",
    "chuT": "Virulence factors",
    "chuW": "Virulence factors",
    "chuX": "Virulence factors",

   

    # Genes diversos
    "blaOKP-B-8_1": "Beta-lactam resistance",
    "blaTEM-1C_1": "Beta-lactam resistance",
    "blaSHV-110_1": "Beta-lactam resistance",
    "oqxB_1": "Efflux pumps",
             # Beta-lactâmicos
    "blaSHV-27_1": "Beta-lactam resistance",
    "blaSHV-187_1": "Beta-lactam resistance",
     "blaDHA-17_1": "Beta-lactam resistance",
    "blaCTX-M-173_1": "Beta-lactam resistance",
     "blaOXA-9_1": "Beta-lactam resistance",
    "blaTEM-141_1": "Beta-lactam resistance",
    "blaKPC-3_1": "Beta-lactam resistance",
    "blaSCO-1_1": "Beta-lactam resistance",

    # Resistência a quinolonas
    "aac(6')-Ib-cr_1": "Quinolone resistance",
    "qnrB4_1": "Quinolone resistance",

    # Resistência a aminoglicosídeos
    "rmtE_1": "Aminoglycoside resistance",
    "rmtF_1": "Aminoglycoside resistance",
    "ARR-3_4": "Aminoglycoside resistance",
    "APH(3')-VI_1": "Aminoglycoside resistance",
    "aph(6)-Ic_1": "Aminoglycoside resistance",
    "aph(3')-IIa_2": "Aminoglycoside resistance",

    # Resistência a sulfonamidas
    "sul2_2": "Sulfonamide resistance",
    "sul1_15": "Sulfonamide resistance",

    # Resistência a trimetoprim
    "dfrA19_1": "Trimethoprim resistance",
    "dfrB1_1": "Trimethoprim resistance",
    "dfrA7_5": "Trimethoprim resistance",

    # Resistência a tetraciclinas
    "tet(D)_1": "Tetracycline resistance",
    "tet(A)_4": "Tetracycline resistance",


    # Resistência a fosfomicina
    "fosA5_1": "Fosfomycin resistance",
    "fosA6_1": "Fosfomycin resistance",
    "fosA_3": "Fosfomycin resistance",

    # Resistência a macrolídeos
    "mph(A)_2": "Macrolide resistance",
    "mefC": "Macrolide resistance",

         # Resistência a colistina
  
    # Fatores de virulência
    "iutA": "Virulence factors",
    "iucC": "Virulence factors",
    "iucB": "Virulence factors",
    "iucA": "Virulence factors",
    "cseA": "Virulence factors",
    "hlyA": "Virulence factors",
    "hlyB": "Virulence factors",
    "hlyD": "Virulence factors",
    "sat": "Virulence factors",
    "papX": "Virulence factors",
    "chuA": "Virulence factors",
    "chuS": "Virulence factors",
    "chuT": "Virulence factors",
    "chuW": "Virulence factors",
    "chuX": "Virulence factors",

    # Bombas de efluxo
    "eefA": "Efflux pumps",
    "ydeO": "Efflux pumps",
    "ygiW": "Efflux pumps",
    "yhcN": "Efflux pumps",

    # Genes diversos
    "blaOKP-B-8_1": "Beta-lactam resistance",
    "blaTEM-1C_1": "Beta-lactam resistance",
    "blaSHV-110_1": "Beta-lactam resistance",
     "blaDHA-1_1": "Beta-lactam resistance",
      "blaSHV-107_1": "Beta-lactam resistance",
    "blaSHV-145_1": "Beta-lactam resistance",
    "blaSHV-106_1": "Beta-lactam resistance",
      "blaSHV-11_1": "Beta-lactam resistance",
    "fosA_5": "Fosfomycin resistance",
     "blaCTX-M-14_1": "Beta-lactam resistance",
    "blaSHV-62_1": "Beta-lactam resistance",
    "aac(6')-Ib_1": "Aminoglycoside resistance",
    "dfrA25_1": "Trimethoprim resistance",
    "blaSHV-190_1": "Beta-lactam resistance",
       "mef(C)_1": "Macrolide resistance",
    "ere(A)_2": "Macrolide resistance",
      "arr-3": "Rifamycin resistance",
    "AAC(6')-Ib4": "Aminoglycoside resistance",
    "SHV-110": "Beta-lactam resistance",
    "SHV-11": "Beta-lactam resistance",
    "DHA-1": "Beta-lactam resistance",
    "QnrB4": "Quinolone resistance",
    "SHV-106": "Beta-lactam resistance",
    "TEM-122": "Beta-lactam resistance",
    "SHV-187": "Beta-lactam resistance",
    "AAC(6')-Ib10": "Aminoglycoside resistance",
    "DHA-17": "Beta-lactam resistance",
      # Efflux pumps
    "MDTA": "Efflux pumps",
    "MDTB": "Efflux pumps",
    "MDTC": "Efflux pumps",
    "BAES": "Efflux pumps",
    "BAER": "Efflux pumps",
    "ACRD": "Efflux pumps",
    "MDTI": "Efflux pumps",
    "MDTJ": "Efflux pumps",
    "MDTK": "Efflux pumps",
    "MARA": "Efflux pumps",
    "MARR": "Efflux pumps",

    # Beta-lactam resistance
    "PBP2": "Beta-lactam resistance",
    "CTX": "Beta-lactam resistance",
    "TEM": "Beta-lactam resistance",

    # Colistin resistance
    "PMRG": "Colistin resistance",

    # Macrolide resistance
    "EMRD": "Macrolide resistance",
    "EMRK": "Macrolide resistance",

    # Quinolone resistance
    "SOXS": "Quinolone resistance",
    "CPXAR": "Quinolone resistance",

    # Rifampicin resistance
    "ROBA": "Rifampicin resistance",

    # Aminoglycoside resistance
    "ASMA": "Aminoglycoside resistance",

    # Tetracycline resistance
    "TETA": "Tetracycline resistance",

    # Bacitracin resistance
    "BCR": "Bacitracin resistance",

    # Sulfonamide resistance
    "SULI": "Sulfonamide resistance",
    "CMLA": "Beta-lactam resistance",  # Novo
    "SHV": "Beta-lactam resistance",   # Novo
    "CARB": "Beta-lactam resistance",  # Novo

    # Colistin resistance
    "PMRG": "Colistin resistance",

    # Tetracycline resistance
    "TETA": "Tetracycline resistance",
    "TETD": "Tetracycline resistance", # Novo

    # Aminoglycoside resistance
    "APH3-DPRIME": "Aminoglycoside resistance",  # Novo
        "MDTG": "Efflux pumps",
    "MDTH": "Efflux pumps",
    "BACA": "Bacitracin resistance",
    "YOGI": "Efflux pumps",
    "MDTE": "Efflux pumps",
    "MDTF": "Efflux pumps",
    "GADX": "Efflux pumps",
    "MDFA": "Efflux pumps",
    "EPTA": "Colistin resistance",
    "EMRY": "Efflux pumps",
    "MPHB": "Macrolide resistance",
    "PMRF": "Biocide resistance",
    "KDPE": "Metal ion resistance",
    "SULII": "Sulfonamide resistance",
    "APH6": "Aminoglycoside resistance",
    "MSBA": "Efflux pumps",
    "AMPH": "Beta-lactam resistance",
    "ACRB": "Efflux pumps",
    "ACRA": "Efflux pumps",
    "HNS": "Efflux pumps",
    "MDTO": "Efflux pumps",
    "MDTN": "Efflux pumps",
    "ACRS": "Efflux pumps",
    "ACRE": "Efflux pumps",
    "ACRF": "Efflux pumps",
    "EMRR": "Efflux pumps",
    "EMRA": "Efflux pumps",
    "EMRB": "Efflux pumps",
    "UGD": "Polymyxin resistance",
    "SULIII": "Sulfonamide resistance",
    "QACL": "Disinfectant resistance",
    "TETB": "Tetracycline resistance",
    "QNRS": "Quinolone resistance",
    "LNUF": "Lincosamide resistance",
    "MCR": "Colistin resistance",
    "FOSA": "Fosfomycin resistance",
    "AAC3": "Aminoglycoside resistance",
    "APH4": "Aminoglycoside resistance",
    "FLOR": "Phenicol resistance",
    "MEFB": "Macrolide resistance",
    "AMPC": "Beta-lactam resistance",
    "NDM": "Carbapenem resistance",
    "BRP": "Biofilm formation",
    "TETM": "Tetracycline resistance",
    "MEFC": "Macrolide resistance",
    "ERMB": "Macrolide resistance",
    "MPHA": "Macrolide resistance",
    "APH3-PRIME": "Aminoglycoside resistance",
    "LNUG": "Lincosamide resistance",
    "QNRB": "Quinolone resistance",
    "MPHE": "Macrolide resistance",
    "MSRE": "Macrolide resistance",
    "ARMA": "Ribosome protection",
    "OQXA": "Quinolone resistance",
    "AAC2-PRIME": "Aminoglycoside resistance",
    "QNRD": "Quinolone resistance",
    "TETC": "Tetracycline resistance",
    "LAP": "Beta-lactam resistance",
    "KPN": "Efflux pumps",
    "MPHG": "Macrolide resistance",
    "SDIA": "Efflux pumps",
    "GESA": "Efflux pumps",
    "IncX4_1": "Large plasmid",
"IncR_1": "Large plasmid",
"IncX1_1": "Large plasmid",
"IncFIB(AP001918)_1": "Large plasmid",
"IncFII_1": "Large plasmid",
"IncB/O/K/Z_4": "Large plasmid",
"RepA_1_pKPC-CAV1321": "Large plasmid",
"IncHI2_1": "Large plasmid",
"IncHI2A_1": "Large plasmid",
"IncFIC(FII)_1": "Large plasmid",
"IncI1_1_Alpha": "Large plasmid",
"IncHI1B(R27)_1_R27": "Large plasmid",
"IncHI1A_1": "Large plasmid",
"IncFIA_1": "Large plasmid",
"IncX4_2": "Large plasmid",
"IncY_1": "Large plasmid",
"IncFII(pRSB107)_1_pRSB107": "Large plasmid",
"IncX3_1": "Large plasmid",
"Col156_1": "Small plasmid",
"ColRNAI_1": "Small plasmid",
"Col(MG828)_1": "Small plasmid",
"ColE10_1": "Small plasmid",
"Col8282_1": "Small plasmid",
"Col440II_1": "Small plasmid",
"Col(KPHS6)_1": "Small plasmid",
"Col440I_1": "Small plasmid",
"p0111_1": "Large plasmid",
"ColpVC_1": "Small plasmid",
"pENTAS02_1": "Large plasmid",
"RepA_1_pKPC-CAV1321": "Large plasmid",
   "IncQ1_1": "Large plasmid",
    "IncB/O/K/Z_2": "Large plasmid",
    "Col(BS512)_1": "Small plasmid",
    "IncHI1B(CIT)_1_pNDM-CIT": "Large plasmid",
    "IncI2_1": "Large plasmid",
    "IncX1_4": "Large plasmid",
    "IncFII(29)_1_pUTI89": "Large plasmid",
    "IncI_Gamma_1": "Large plasmid",
    "IncFIB(Mar)_1_pNDM-Mar": "Large plasmid",
    "IncHI1B_1_pNDM-MAR": "Large plasmid",
    "IncFII(pHN7A8)_1_pHN7A8": "Large plasmid",
    "IncI2_1_Delta": "Large plasmid",
    "Col(MP18)_1": "Small plasmid",
    "IncFIB(K)_1_Kpn3": "Large plasmid",
    "IncFIB(pHCM2)_1_pHCM2": "Large plasmid",
    "IncFIA(HI1)_1_HI1": "Large plasmid",
    "IncFII(pSE11)_1_pSE11": "Large plasmid",
    "IncN_1": "Large plasmid",
    "IncFIB(pECLA)_1_pECLA": "Large plasmid",
    "IncFII(pCoo)_1_pCoo": "Large plasmid",
    "IncL/M(pOXA-48)_1_pOXA-48": "Large plasmid",
    "Col3M_1": "Small plasmid",
    "Col(Ye4449)_1": "Small plasmid",
    "IncFII(pHN7A8)_1_pHN7A8": "Large plasmid",
    "IncI2_1_Delta": "Large plasmid",
    "IncX1_4": "Large plasmid",
    "IncFIA(HI1)_1_HI1": "Large plasmid",
    "IncFIB(K)_1_Kpn3": "Large plasmid",
    "IncI2_1": "Large plasmid",
    "IncFII(pCoo)_1_pCoo": "Large plasmid",
    "IncB/O/K/Z_2": "Large plasmid",
    "IncFII(29)_1_pUTI89": "Large plasmid",
    "IncI_Gamma_1": "Large plasmid",
    "IncFIB(pB171)_1_pB171": "Large plasmid",
    "Col(BS512)_1": "Small plasmid",
    "IncHI1B(CIT)_1_pNDM-CIT": "Large plasmid",
    "IncX2_1": "Large plasmid",
    "IncL/M(pOXA-48)_1_pOXA-48": "Large plasmid",
    "IncFII(Yp)_1_Yersenia": "Large plasmid",
    "Col(Ye4449)_1": "Small plasmid",
    "Col3M_1": "Small plasmid",
    "IncFII(SARC14)_1_SARC14": "Large plasmid",
    "IncFII(pAMA1167-NDM-5)_1_pAMA1167-NDM-5": "Large plasmid",
    "Col(VCM04)_1": "Small plasmid",
    "IncP1_3": "Large plasmid",
    "IncB/O/K/Z_3": "Large plasmid",
    "IncX2_1": "Large plasmid",
    "IncFII_1_pSFO": "Large plasmid",
   "IncFIB(pENTAS01)_1_pENTAS01": "Large plasmid",
   "IncFII(S)_1": "Large plasmid",
   "IncFIB(S)_1": "Large plasmid",
   "IncFIB(pCTU1)_1_pCTU1": "Large plasmid",
   "ASMA": "Aminoglycoside resistance",
    "APH3-DPRIME": "Aminoglycoside resistance",  # Novo
    "RMTF": "Aminoglycoside resistance",         # Novo

    # Tetracycline resistance
    "TETA": "Tetracycline resistance",
    "TETD": "Tetracycline resistance", # Novo

    # Bacitracin resistance
    "BCR": "Bacitracin resistance",
    "BACA": "Bacitracin resistance",  # Novo

    # Sulfonamide resistance
    "SULI": "Sulfonamide resistance",

    # Beta-lactam resistance
    "CMLA": "Beta-lactam resistance",  # Novo
    "PBP2": "Beta-lactam resistance",
    "CTX": "Beta-lactam resistance",
    "SHV": "Beta-lactam resistance",   # Novo
    "CARB": "Beta-lactam resistance",  # Novo
    "DFRB": "Trimethoprim resistance",  # Novo, trimethoprim

    # Colistin resistance
    "PMRG": "Colistin resistance",

    # Rifamycin resistance
    "ARR": "Rifamycin resistance",     # Novo

    # Carbapenem resistance
    "VIM": "Carbapenem resistance",    # Novo

    # Efflux pumps
    "MDTG": "Efflux pumps",
    "MDTH": "Efflux pumps",
    "MDTE": "Efflux pumps",
    "MDTF": "Efflux pumps",
    "MDFA": "Efflux pumps",
    "YOGI": "Efflux pumps",
    "GADX": "Efflux pumps",
    "KMRA": "Efflux pumps",             # Novo (MFS-type)
    "KDEA": "Efflux pumps",             # Novo
    "OQXB": "Efflux pumps",    
  # Porin alteration / Outer membrane protein (related to resistance)
    "OMP37": "Porin loss-associated resistance",  # Novo
    "KPNO": "Porin loss-associated resistance",   # Novo
"PHOR": "Efflux pumps",  # Novo
    "KPNE": "Virulence factors",  # Novo
    "KPNF": "Efflux pumps",  # Novo   
# Resistência a beta-lactâmicos
    "DHA": "Beta-lactam resistance",
    "KPC": "Beta-lactam resistance",

    # Resistência a macrolídeos
    "EREA": "Macrolide resistance",

    # Plasmídeos de grande porte associados à resistência
    "IncFII_1_pKP91": "Large plasmid",
    "FII(pBK30683)_1": "Large plasmid",
    "IncFII(pCRY)_1_pCRY": "Large plasmid",
    "IncFIB(pQil)_1_pQil": "Large plasmid",
    "IncFIB(pKPHS1)_1_pKPHS1": "Large plasmid",
# Resistência a beta-lactâmicos
    "SCO": "Beta-lactam resistance",

    # Plasmídeos de grande porte associados à resistência
    "IncFII(K)_1": "Large plasmid",
    "IncA/C2_1": "Large plasmid",
    "IncFII(pMET)_1_pMET1": "Large plasmid",
  "rep7a_16_repC(Cassette)": "Small plasmid",
"rep13_4_rep(pKH13)": "Small plasmid",
"rep21_1_rep(pWBG754)": "Small plasmid",
"rep21_17_rep(pLNU8)": "Small plasmid",
"rep21_3_rep(pS0385)": "Small plasmid",
"rep21_25_rep(pSHaeA)": "Small plasmid",
    



}

    
    



def extract_sample_name(file_name):
    """Extrai o nome da amostra do nome do arquivo, removendo:
    - Extensões (.tsv, .fasta)
    - Sufixo '_plasmidfinder' (se presente)
    """
    # Remove extensões (.tsv, .fasta, etc.)
    base_name = file_name.split('.')[0].replace('.fasta', '').replace('.tsv', '')
    # Remove '_plasmidfinder' apenas se estiver presente
    if '_plasmidfinder' in base_name:
        return base_name.split('_plasmidfinder')[0]
    return base_name

def process_plasmidial_file(plasmidial_path, gene_to_group):
    """Processa arquivos plasmidiais com tratamento robusto de colunas"""
    try:
        df = pd.read_csv(plasmidial_path, sep='\t')
        
        # Mapeamento completo de colunas (incluindo Source_File e variações)
        column_mapping = {
            'SAMPLE': ['File', 'Sample', 'SAMPLE', 'Filename', 'Source_File', 'FILE'],
            'GENE': ['GENE', 'Gene', 'gene'],
            'SEQUENCE': ['SEQUENCE', 'Sequence', 'sequence', 'SEQUENCE_ID'],
            'START': ['START', 'Start', 'start', 'Inicio'],
            'END': ['END', 'End', 'end', 'Fim']
        }

        # Padronização de colunas
        for standard_col, possible_cols in column_mapping.items():
            found = False
            for col in possible_cols:
                if col in df.columns:
                    df.rename(columns={col: standard_col}, inplace=True)
                    found = True
                    break
            
            if not found:
                if standard_col == 'SAMPLE':
                    raise ValueError(f"Coluna de amostra não encontrada em {plasmidial_path}. Colunas existentes: {df.columns.tolist()}")
                df[standard_col] = None

        # Fallback explícito para Source_File
        if 'SAMPLE' not in df.columns and 'Source_File' in df.columns:
            df['SAMPLE'] = df['Source_File']

        # Processamento final
        df['SAMPLE'] = df['SAMPLE'].apply(extract_sample_name)
        df['Group'] = df['GENE'].map(gene_to_group).fillna('Unclassified')
        
        return df.groupby(['SAMPLE', 'Group']).size().reset_index(name='Gene_Count')

    except Exception as e:
        print(f"\nERRO no processamento de {plasmidial_path}")
        print(f"Tipo: {type(e).__name__}")
        print(f"Detalhes: {str(e)}")
        return None

def process_file(file_path, gene_to_group):
    """Processa arquivos regulares com tratamento especial para BacMet"""
    try:
        data = pd.read_csv(file_path, sep='\t')
        
        # Mapeamento dinâmico de colunas
        column_mapping = {
            'SAMPLE': ['File', 'Sample', 'SAMPLE', 'Filename', 'Source_File'],
            'GENE': ['GENE', 'Gene', 'gene'],
            'SEQUENCE': ['SEQUENCE', 'Sequence', 'SEQUENCE_ID', 'Query'],
            'START': ['START', 'Start', 'start'],
            'END': ['END', 'End', 'end']
        }

        # Padronização
        for standard_col, possible_cols in column_mapping.items():
            found = False
            for col in possible_cols:
                if col in data.columns:
                    data.rename(columns={col: standard_col}, inplace=True)
                    found = True
                    break
            
            if not found:
                if standard_col in ['SAMPLE', 'GENE']:
                    raise ValueError(f"Coluna obrigatória {standard_col} não encontrada em {file_path}")
                data[standard_col] = None

        # Processamento
        data['SAMPLE'] = data['SAMPLE'].apply(extract_sample_name)
        data['Group'] = data['GENE'].map(gene_to_group).fillna('Unclassified')
        
        return data[['SAMPLE', 'Group', 'GENE', 'SEQUENCE', 'START', 'END']]

    except Exception as e:
        print(f"Erro ao processar {file_path}: {str(e)}")
        return pd.DataFrame()

def deduplicate_data(df):
    """Remove duplicatas considerando:
    - Case insensitive para genes
    - Mesmo grupo de resistência
    - Mesma amostra
    """
    dedup_data = []
    seen = {}
    
    for _, row in df.iterrows():
        # Normalização
        sample = row['SAMPLE']
        normalized_gene = str(row['GENE']).lower() if pd.notna(row['GENE']) else None
        group = row['Group']  # Já mapeado para "Efflux pumps", etc.
        
        # Chave única: Amostra + Gene (case insensitive) + Grupo
        primary_key = (sample, normalized_gene, group)
        
        # Chave secundária: Amostra + Posição genômica (opcional)
        secondary_key = (sample, row.get('SEQUENCE'), row.get('START'), row.get('END'), group) if all(col in df.columns for col in ['SEQUENCE', 'START', 'END']) else None
        
        # Verificação de duplicatas
        if primary_key not in seen and (secondary_key is None or secondary_key not in seen):
            dedup_data.append(row)
            seen[primary_key] = True
            if secondary_key is not None:
                seen[secondary_key] = True
                
    return pd.DataFrame(dedup_data)

# Processamento principal
for config in bacteria_config:
    print(f"\nProcessando {config['name']}...")
    
    base_dir = config['base_dir']
    files = [os.path.join(base_dir, f) for f in config['file_names'] if os.path.exists(os.path.join(base_dir, f))]
    plasmidial = os.path.join(base_dir, config['plasmidial_file'])
    
    try:
        # Processar arquivos principais
        all_data = pd.concat([process_file(f, gene_to_group) for f in files])
        dedup_data = deduplicate_data(all_data)
        
        # Processar plasmidial
        if os.path.exists(plasmidial):
            plasmid_counts = process_plasmidial_file(plasmidial, gene_to_group)
            if plasmid_counts is not None:
                plasmid_counts.to_csv(os.path.join(base_dir, 'plasmidial_counts.tsv'), sep='\t', index=False)
        
        # Gerar resultados
        group_counts = dedup_data.groupby(['SAMPLE', 'Group']).size().reset_index(name='Gene_Count')
        group_counts.to_csv(os.path.join(base_dir, 'group_counts_total.tsv'), sep='\t', index=False)
        
        unclassified = dedup_data[dedup_data['Group'] == 'Unclassified']
        unclassified.to_csv(os.path.join(base_dir, 'unclassified_genes.tsv'), sep='\t', index=False)
        
    except Exception as e:
        print(f"Erro no processamento de {config['name']}: {str(e)}")

print("\nProcessamento concluído para todas as bactérias!")