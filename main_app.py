# -*- coding: utf-8 -*-
"""
main.py
Script principal que orquestra el pipeline de dades bioinformàtiques i anàlisi.
"""
import time
import networkx as nx
from Bio import SeqIO

# Importacions de la nostra arquitectura modular
from parsing_utils import TF_RISet_parse
from network_metrics import degree_distribution
from motifs_analysis import compare_with_random_graphs

def run_pipeline():
    GENOME_FILE = 'U00096.3.gb'
    TF_RISET_FILE = 'TF-RISet.tsv'
    TF_SET_FILE = 'TFSet.tsv'

    print("==========================================================")
    print(" E. coli Transcriptional Regulatory Network (TRN) Builder ")
    print("==========================================================")
    start_time = time.time()

    try:
        print(f"\n[1/4] Parsejant anotacions genòmiques ({GENOME_FILE})...")
        genoma = SeqIO.read(GENOME_FILE, 'genbank')
    except Exception as e:
        print(f"Error crític: {e}. Comproveu que els arxius de dades hi són presents.")
        return

    print("[2/4] Construint les xarxes dirigides (grafs)...")
    # Generem xarxa a nivell de gens individuals
    G_gene = TF_RISet_parse(TF_RISET_FILE, TF_SET_FILE, False, 100, genoma)
    # Generem xarxa avaluant operons complets (distància < 100 pb)
    G_operon = TF_RISet_parse(TF_RISET_FILE, TF_SET_FILE, True, 100, genoma)

    print("\n--- Estadístiques Topològiques Bàsiques ---")
    print(f"Xarxa (Resolució Gen):    {G_gene.number_of_nodes()} nodes, {G_gene.number_of_edges()} interaccions")
    print(f"Xarxa (Resolució Operó):  {G_operon.number_of_nodes()} nodes, {G_operon.number_of_edges()} interaccions")

    # Guardar els grafs en format estàndard
    nx.write_graphml(G_gene, 'Ecoli_TRN.graphml')
    nx.write_graphml(G_operon, 'Ecoli_operon_TRN.graphml')

    print("\n[3/4] Identificant 'Master Regulators' (Hubs Centrals)...")
    # Extraiem el top 5 de factors de transcripció (out-degree)
    _, top_out = degree_distribution(G_operon, 'out', 'TF', 5)
    for nom, grau in top_out.items():
        print(f"  -> El Factor de Transcripció {nom} regula directament {grau} dianes.")

    print("\n[4/4] Validant Motius de Xarxa contra Model Nul...")
    compare_with_random_graphs(G_gene, n_random=10)

    print(f"\nPipeline finalitzat correctament en {time.time() - start_time:.2f} segons.")

if __name__ == "__main__":
    run_pipeline()