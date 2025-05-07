#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: rafael

Script otimizado para análises de correlação fago-genoma:
  1. phage-gene (Negative Binomial)
  2. gene-gene (Negative Binomial)
  3. phage-phage (Fisher)
  4. phage-combo (Negative Binomial)

 - Gera apenas arquivos .tsv filtrados por p < 0.05 e coef > 0
 - Processamento em paralelo com até 60 threads
 - Salva resultados organizados por bactéria, tipo de análise e grupo genômico
"""

import os
import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count
from scipy.stats import fisher_exact
import statsmodels.api as sm  # <- ESSENCIAL para usar sm.GLM()


# ─────────────────────────── CONFIGURAÇÃO ───────────────────────────
# ───────────────────────── CONFIGURAÇÃO ─────────────────────────
# Tamanhos de cluster como fração do maior genoma (ex: 0.2 para 5 grupos)
CLUSTER_FRAC = 0.2

blocks = [
    "Antibiotic Resistance",
    "Metal Resistance",
    "Toxins",
    "Stress Resistance",
    "Pathogenicity",
    "Virulence",
    "Plasmids",
    "Other"
]


# Dicionário para mapear grupos de genes (Group) a blocos funcionais
gene_group_dict = {
       "Aminoglycoside resistance": "Antibiotic Resistance",
       "Bacitracin resistance": "Antibiotic Resistance",
    "Antibiotic Resistance Regulator": "Antibiotic Resistance",
    "Polimyxin resistance": "Antibiotic Resistance",
    "Bacitracin": "Antibiotic Resistance",
    "Quinolones resistance": "Antibiotic Resistance",
    "Beta-Lactam Resistance": "Antibiotic Resistance",
    "Bleomycin Resistance": "Antibiotic Resistance",
    "Carbapenemases": "Antibiotic Resistance",
    "Chloramphenicol Resistance": "Antibiotic Resistance",
    "Colistin Resistance": "Antibiotic Resistance",
    "Extended-Spectrum Beta-Lactamase (Esbl)": "Antibiotic Resistance",
    "Fosfomycin resistance": "Antibiotic Resistance",
    "Lincosamide Resistance": "Antibiotic Resistance",
    "Macrolide Resistance": "Antibiotic Resistance",
    "Nitrofurantoin Resistance": "Antibiotic Resistance",
    "Phenicol resistance": "Antibiotic Resistance",
    "Polymyxins": "Antibiotic Resistance",
    "Quinolone resistance": "Antibiotic Resistance",
    "Quinolones": "Antibiotic Resistance",
    "Tellurite Resistance": "Antibiotic Resistance",
    "Tellurium Resistance": "Antibiotic Resistance",
    "Tetracycline Resistance": "Antibiotic Resistance",
      "Aminoglycosides": "Antibiotic Resistance",
"Aminoglycosides resistance": "Antibiotic Resistance",
"AmpC beta-lactamase resistance": "Antibiotic Resistance",
"Bleomycin resistance": "Antibiotic Resistance",
"Cecropin resistance": "Antibiotic Resistance",
"Colistin resistance": "Antibiotic Resistance",
"Extended-spectrum Beta-lactamase (ESBL)": "Antibiotic Resistance",
"Fosfomycin resistance": "Antibiotic Resistance",
"Leucocidin": "Antibiotic Resistance",
"Lipopeptides resistance": "Antibiotic Resistance",
"Methicillin resistance": "Antibiotic Resistance",
"Mupirocin resistance": "Antibiotic Resistance",
"Nitrofurantoin resistance": "Antibiotic Resistance",
"Oxazolidinone resistance": "Antibiotic Resistance",
"Beta-lactam resistance": "Antibiotic Resistance",
"Macrolide resistance": "Antibiotic Resistance",
"Lincosamide resistance": "Antibiotic Resistance",
"Leucocidin resistance": "Antibiotic Resistance",
"Rifampin resistance": "Antibiotic Resistance",
"Rifamycin resistance": "Antibiotic Resistance",
"Sulfonamide resistance": "Antibiotic Resistance",
"Tellurite resistance": "Antibiotic Resistance",
"Tetracycline resistance": "Antibiotic Resistance",
"Trimethoprim resistance": "Antibiotic Resistance",
"Carbapenem resistance": "Antibiotic Resistance",
"Folate pathway resistance": "Antibiotic Resistance",
"Folate synthesis inhibitors": "Antibiotic Resistance",
"Fusidic acid": "Antibiotic Resistance",
"Glycopeptide resistance": "Antibiotic Resistance",
"Streptothricin resistance": "Antibiotic Resistance",
"Rifampicin resistance": "Antibiotic Resistance",
"Polymyxin resistance": "Antibiotic Resistance",

"Large plasmid": "Plasmids",
"Small plasmid": "Plasmids",

    "Copper Resistance": "Metal Resistance",
    "Ion Transporter": "Metal Resistance",
    "Iron Acquisition": "Metal Resistance",
    "Iron Acquisition System": "Metal Resistance",
    "Iron And Zinc Efflux": "Metal Resistance",
    "Iron Siderophore System": "Metal Resistance",
    "Iron Transport": "Metal Resistance",
    "Iron Uptake": "Metal Resistance",
    "Magnesium Transport": "Metal Resistance",
    "Manganese And Iron Transport": "Metal Resistance",
    "Manganese Efflux": "Metal Resistance",
    "Mercury Resistance": "Metal Resistance",
    "Metal Resistance": "Metal Resistance",
    "Metal Resistance (Arsenic)": "Metal Resistance",
    "Zinc Transport": "Metal Resistance",
    "Zinc Transporter": "Metal Resistance",
    "Copper resistance":"Metal Resistance",
    "Copper transport": "Metal Resistance",
"Copper/Silver Resistance": "Metal Resistance",
"Zinc transporter": "Metal Resistance",
    
    "Pathogenicity": "Pathogenicity",
    "Pathogenicity Island (SPI-2)": "Pathogenicity",

    "Enteropathogenic E. Coli": "Toxins",
    "Enterotoxin": "Toxins",
    "Exotoxin": "Toxins",
    "Shiga Toxin": "Toxins",
    "Shigella toxin": "Toxins",
    "Toxin production": "Toxins",

    "Acid Resistance": "Stress Resistance",
    "Disinfectant resistance": "Stress Resistance",
    "Heat shock proteins": "Stress Resistance",
    "Stress response": "Stress Resistance",
    "Stress response Regulators": "Stress Resistance",
    "Efflux Pumps": "Stress Resistance",
    "Acid resistance": "Stress Resistance",
    "Acid stress response": "Stress Resistance",
    "Disinfectant resistant": "Stress Resistance",
    "Efflux pumps": "Stress Resistance",
    "Environmental stress sensor": "Stress Resistance",
    "Oxidative stress resistance": "Stress Resistance",
    "Stress Resistanc": "Stress Resistance",
    "Biocide resistance": "Stress Resistance",


    "Siderophore Biosynthesis": "Nutrient Transport",
    "Siderophore Transport": "Nutrient Transport",
    "Phosphate Transport": "Nutrient Transport",
    "Formate Transport": "Nutrient Transport",
    "Glycerol Transport": "Nutrient Transport",

    "Adhesion": "Virulence",
    "Biofilm Regulators": "Virulence",
    "Capsule Biosynthesis": "Virulence",
    "Cytotoxic necrotizing factor": "Virulence",
    "Cytotoxin": "Virulence",
    "Fimbrial Adhesin": "Virulence",
    "Fimbrial Assembly": "Virulence",
    "Hemolysin": "Virulence",
    "Invasion": "Virulence",
    "Toxin Production": "Virulence",
    "Flagella synthesis": "Virulence",
    "Motility": "Virulence",
    "Adhesion and biofilm formation": "Virulence",
    "Biofilm regulators": "Virulence",
    "Capsule biosynthesis": "Virulence",
    "Fimbriae formation and adhesion": "Virulence",
    "Virulence factors": "Virulence",
    "Type III Secretion System (T3SS)": "Virulence",
    "Type VI secretion system":"Virulence",
    "Secretion System": "Virulence",
    "Curly Assembly": "Virulence",
    "EscE Secretion": "Virulence",
    "Two-component systems": "Virulence",
    "Quorum sensing": "Virulence",
    "Curly assembly": "Virulence",
    "Porin loss-associated resistance": "Virulence",
    "Biofilm formation": "Virulence",
    

    "16S RRNA Methyltransferase": "Antibiotic Resistance",
    "Protein Disulfide Isomerase": "Other",
    "Two-Component Systems": "Virulence",
    "Outer Membrane Proteins": "Virulence",
    "Global Regulators": "Virulence",
    "Effector protein": "Virulence",
        "Plasmid-Encoding Factors": "Virulence",
    "Environmental Stress Sensor": "Stress Resistance",
    "Dna Repair": "Other"
}

# Limiares de tamanho de genoma bacteriano (bp)
GENOME_THRESHOLDS = {
    'Grupo1': 6500645,
    'Grupo2': 5850580,
    'Grupo3': 5200516,
    'Grupo4': 4550451,
    'Grupo5': 3900387
}

def assign_genome_group(size):
    """Classifica tamanho de genoma em Grupo1 a Grupo5 ou Menor60."""
    for group, thr in GENOME_THRESHOLDS.items():
        if size >= thr:
            return group
    return 'Menor60'

# Diretório base de saída
base_out = "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Results/Poisson/"
# Subdiretórios por tipo de análise
dirs = {
    'phage-phage': os.path.join(base_out, 'co_occurrence'),
    'phage-combo': os.path.join(base_out, 'co_occurrence'),
    'gene-gene': os.path.join(base_out, 'gene_vs_gene'),
    'phage-gene': os.path.join(base_out, 'phage_vs_gene')
}
for d in dirs.values():
    os.makedirs(d, exist_ok=True)

paths = {
    "E. coli": {
        "genes": "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Diversity_Results_Rarefied/Escherichia_coli_filtered_genes.tsv",
        "phages": "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Phages/E.coli/E_coli_count_total_phages_sizes.tsv"
    },
    "Klebsiella": {
        "genes": "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Diversity_Results_Rarefied/Klebsiella_pneumoniae_filtered_genes.tsv",
        "phages": "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Phages/Klebsiella_p/Klebsiella_count_total_phages_sizes.tsv"
    },
    "Salmonella enterica": {
        "genes": "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Diversity_Results_Rarefied/Salmonella_Enterica_filtered_genes.tsv",
        "phages": "/media/rafael/dbdad3a1-6923-4ffa-a1e4-30a98d3941a5/Taxonomia/WGS_BACTERIA/Phages/Salmonella_enterica/Salmonella_count_total_phages_sizes.tsv"
    }
}

ANALYSIS_TYPES = ['phage-gene', 'gene-gene', 'phage-phage', 'phage-combo']

# ────────────────── FUNÇÕES AUXILIARES ──────────────────

def map_gene_to_category(name):
    return gene_group_dict.get(name, np.nan)


def load_data(files):
    # Genes
    gdf = pd.read_csv(files['genes'], sep='\t')
    gdf['Category'] = gdf['Group'].map(map_gene_to_category)
    gdf = gdf.dropna(subset=['Category'])

    # Phages + genome size
    raw = pd.read_csv(files['phages'], sep='\t')
    raw['Count'] = pd.to_numeric(raw['Count'], errors='coerce').fillna(0)
    # captura genome size por sample (único valor)
    if 'Bacterial_Genome_Size' in raw.columns:
        raw['Bacterial_Genome_Size'] = pd.to_numeric(raw['Bacterial_Genome_Size'], errors='coerce')
        genome = raw.groupby('Sample')['Bacterial_Genome_Size'].first()
    else:
        genome = None
    grouped = raw.groupby('Sample').agg({'Genus': list, 'Count': list}).reset_index()
    records = []
    for _,row in grouped.iterrows():
        rec = {'Sample':row.Sample}
        for g,c in zip(row.Genus,row.Count): rec[g]=c
        records.append(rec)
    pdf = pd.DataFrame(records).set_index('Sample').fillna(0).apply(pd.to_numeric)

    # cria clusters dinâmicos baseado no maior genoma
    if genome is not None:
        max_size = genome.max()
        n_clusters = int(1/CLUSTER_FRAC)
        pdf['Bacterial_Genome_Size'] = genome
        def assign_cluster(sz):
            for i in range(n_clusters-1):
                if sz >= max_size*(1 - (i+1)*CLUSTER_FRAC):
                    return f"G{i+1}"  # G1 = top 0-CLUSTER_FRAC
            return f"G{n_clusters}"  # último grupo
        pdf['Genome_Cluster'] = pdf['Bacterial_Genome_Size'].apply(assign_cluster)
    return gdf, pdf


def prepare_block_data(gdf,block):
    sub = gdf[gdf['Category']==block]
    pivot = sub.groupby(['Sample','Group'])['Gene_Count'].sum().unstack().fillna(0)
    pivot.index.name='Sample'
    return pivot.apply(pd.to_numeric)


def run_nb(y,X,verbose=False):
    common = y.index.intersection(X.index)
    y, X = y.loc[common], X.loc[common]
    X = sm.add_constant(X.astype(float))
    y = pd.to_numeric(y,errors='coerce').fillna(0)
    try:
        m = sm.GLM(y,X,family=sm.families.NegativeBinomial()).fit(disp=0)
        disp = m.deviance/m.df_resid if m.df_resid>0 else np.nan
        if verbose and disp>1.5:
            print(f"⚠️ Superdispersão (disp={disp:.2f})")
        return m
    except Exception as e:
        print(f"⚠️ Erro NB: {e}")
        return None


def process_task(args):
    bact,atype,block,files = args
    gdf,pdf = load_data(files)
    out=[]
    try:
        if atype in ['phage-gene','gene-gene','phage-combo']:
            pivot = prepare_block_data(gdf,block)
            # phages or genes
            items = pivot.columns if atype!='gene-gene' else pivot.columns
            for gene in items:
                for ph in pdf.columns.drop(['Bacterial_Genome_Size','Genome_Cluster'],errors='ignore'):
                    y = pivot[gene]
                    # constrói X com ph, genome size e cluster
                    cols = [ph]
                    if 'Bacterial_Genome_Size' in pdf: cols.append('Bacterial_Genome_Size')
                    X = pdf[cols].astype(float)
                    # dummies cluster
                    if 'Genome_Cluster' in pdf:
                        d = pd.get_dummies(pdf['Genome_Cluster'],prefix='Cluster',drop_first=True)
                        X = pd.concat([X,d],axis=1)
                    if y.sum()==0 or X.values.sum()==0: continue
                    m = run_nb(y,X)
                    if m:
                        # interação para combo
                        var = 'interaction' if atype=='phage-combo' else ph
                        coef = m.params.get(var,np.nan)
                        pval = m.pvalues.get(var,1)
                        if coef>0 and pval<0.05:
                            out.append({'analysis_type':atype,'bacteria':bact,'block':block,
                                        'gene':gene,'phage':var,'coef':coef,'p_value':pval})
        elif atype=='phage-phage':
            for i in range(len(pdf.columns)):
                for j in range(i+1,len(pdf.columns)):
                    a = pdf[pdf.columns[i]]>0; b=pdf[pdf.columns[j]]>0
                    if (a&b).sum()<6: continue
                    t=pd.crosstab(a,b)
                    if t.shape==(2,2):
                        orr,p=fisher_exact(t)
                        if p<0.05 and not np.isinf(orr):
                            out.append({'analysis_type':atype,'bacteria':bact,
                                        'phage1':pdf.columns[i],'phage2':pdf.columns[j],
                                        'odds_ratio':orr,'p_value':p})
    except Exception as e:
        print(f"[ERRO] {bact}-{atype}: {e}")
    return out

if __name__=='__main__':
    tasks=[]
    for bact,f in paths.items():
        for at in ANALYSIS_TYPES:
            if at=='phage-phage': tasks.append((bact,at,None,f))
            else:
                for block in blocks: tasks.append((bact,at,block,f))
    with Pool(cpu_count()-1) as pool:
        results=pool.map(process_task,tasks)
    all_res=[r for sub in results for r in sub]
    df=pd.DataFrame(all_res)
    # salva por bactéria e tipo
    for bact,f in paths.items():
        outd=os.path.join(base_out,bact.replace(' ','_'))
        os.makedirs(outd,exist_ok=True)
        for at in ANALYSIS_TYPES:
            sub=df[(df.bacteria==bact)&(df.analysis_type==at)] if not df.empty else pd.DataFrame()
            sub.to_csv(os.path.join(outd,f"{at}.tsv"),sep='\t',index=False)
    print("✅ Análises finalizadas, resultados salvos por bactéria.")
