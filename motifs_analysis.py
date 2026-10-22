# -*- coding: utf-8 -*-
"""
motifs_analysis.py
Mòdul avançat per a la detecció de motius de xarxa (Feed-Forward Loops, Bi-fans)
i simulacions estadístiques comparatives contra models nuls (Random Graphs).
"""
import networkx as nx
import numpy as np
from typing import List, Tuple

def find_ffl_motifs(G: nx.DiGraph) -> List[Tuple[str, str, str]]:
    ffl = []
    for a in G.nodes():
        succ_a = set(G.successors(a))
        for b in succ_a:
            succ_b = set(G.successors(b))
            common = succ_a.intersection(succ_b)
            for c in common:
                if a != b and b != c and a != c: 
                    ffl.append((a, b, c))
    return ffl

def compare_with_random_graphs(G: nx.DiGraph, n_random: int):
    print(f"\n--- Analitzant Graf ({len(G)} nodes, {len(G.edges)} arestes) ---")
    ffl_orig = find_ffl_motifs(G)
    print(f"[Original] Motius FFL detectats: {len(ffl_orig)}")
    
    in_degs = [d for n, d in G.in_degree()]
    out_degs = [d for n, d in G.out_degree()]
    ffl_rand = []
    
    print(f"Generant {n_random} Grafs Aleatoris segons el model de configuració...")
    for i in range(n_random):
        R_multi = nx.directed_configuration_model(in_degs, out_degs, seed=i)
        R = nx.DiGraph(R_multi)
        R.remove_edges_from(nx.selfloop_edges(R))
        ffl_rand.append(len(find_ffl_motifs(R)))
        
    mitjana = np.mean(ffl_rand)
    desviacio = np.std(ffl_rand)
    print(f"[Model Nul] Mitjana FFL (Grafs Aleatoris): {mitjana:.2f} +/- {desviacio:.2f}")
    
    if len(ffl_orig) > mitjana + (2 * desviacio):
        print(">> RESULTAT: La presència de FFL és estadísticament significativa (per sobre del soroll aleatori).")