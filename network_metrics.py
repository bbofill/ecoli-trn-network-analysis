# -*- coding: utf-8 -*-
"""
network_metrics.py
Mòdul per a l'anàlisi topològica de grafs, incloent distribucions de grau,
càlcul de distàncies mitjanes i components connexes.
"""
import random
import networkx as nx

def largest_CC_graph(G: nx.Graph) -> nx.Graph:
    components = list(nx.weakly_connected_components(G)) if G.is_directed() else list(nx.connected_components(G))
    if not components: return nx.Graph()
    largest_comp = max(components, key=len)
    return G.subgraph(largest_comp).copy()

def degree_distribution(G: nx.DiGraph, degree_type: str, node_type: str, bestN: int) -> tuple:
    deg_dict = dict(G.in_degree()) if degree_type == 'in' else dict(G.out_degree()) if degree_type == 'out' else dict(G.degree())
    filtered = {}
    
    for node, deg in deg_dict.items():
        is_tf = G.out_degree(node) > 0
        is_tg = G.in_degree(node) > 0
        if node_type == 'all' or (node_type == 'TF' and is_tf) or (node_type == 'TG' and is_tg):
            filtered[node] = deg

    if not filtered: return [], {}
    
    max_deg = max(filtered.values())
    freq = [0] * (max_deg + 1)
    for deg in filtered.values(): freq[deg] += 1
    
    best = dict(sorted(filtered.items(), key=lambda x: x[1], reverse=True)[:bestN])
    return freq, best

def average_distance(G: nx.Graph, iterations: int) -> float:
    total_distance, valid_paths = 0, 0
    nodes = list(G.nodes())
    if len(nodes) < 2: return 0.0

    for _ in range(iterations):
        u, v = random.sample(nodes, 2)
        if nx.has_path(G, u, v):
            total_distance += nx.shortest_path_length(G, source=u, target=v)
            valid_paths += 1
            
    return (total_distance / valid_paths) if valid_paths > 0 else 0.0