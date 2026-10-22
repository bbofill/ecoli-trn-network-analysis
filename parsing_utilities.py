# -*- coding: utf-8 -*-
"""
parsing_utils.py
Mòdul encarregat de l'extracció de dades genòmiques, processament d'arxius TSV
i predicció d'operons a partir de dades GenBank.
"""
import csv
import networkx as nx
from Bio.SeqRecord import SeqRecord

def _trobar_index_caracteristica(locus_tag: str, genoma: SeqRecord) -> int:
    for i, carac in enumerate(genoma.features):
        if carac.type == 'gene' and 'locus_tag' in carac.qualifiers:
            if carac.qualifiers['locus_tag'][0] == locus_tag:
                return i
    return -1

def gene_qualifier(consulta: str, camp_consulta: str, camp_objectiu: str, genoma: SeqRecord) -> tuple:
    num_carac = 0
    valor_retorn = ''
    max_caracteristiques = len(genoma.features)
    
    while valor_retorn == '' and num_carac < max_caracteristiques:
        carac_actual = genoma.features[num_carac]
        if carac_actual.type == 'gene':
            if camp_consulta in carac_actual.qualifiers:
                if carac_actual.qualifiers[camp_consulta][0] == consulta:
                    if camp_objectiu in carac_actual.qualifiers:
                        valor_retorn = carac_actual.qualifiers[camp_objectiu][0]
        if valor_retorn == '':
            num_carac += 1
    return num_carac, valor_retorn

def operon(locus_tag: str, dist_intergenica_max: int, genoma: SeqRecord) -> list:
    llista_opero = []
    index_actual = _trobar_index_caracteristica(locus_tag, genoma)
    if index_actual == -1: return llista_opero

    llista_opero.append(locus_tag)
    carac_actual = genoma.features[index_actual]
    cadena_lider = carac_actual.location.strand
    i = index_actual + 1
    max_caracteristiques = len(genoma.features)
    
    while i < max_caracteristiques:
        carac_seguent = genoma.features[i]
        if carac_seguent.type != 'gene' or 'locus_tag' not in carac_seguent.qualifiers:
            i += 1; continue 
        if carac_seguent.location.strand != cadena_lider: break

        dist = (carac_seguent.location.start - carac_actual.location.end) if cadena_lider == 1 else (carac_actual.location.start - carac_seguent.location.end)
        
        if dist <= dist_intergenica_max:
            llista_opero.append(carac_seguent.qualifiers['locus_tag'][0])
            carac_actual = carac_seguent
            i += 1
        else: break
    return llista_opero

def TF_RISet_parse(tf_riset_file: str, tf_set_file: str, detectar_operons: bool, max_dist: int, genoma: SeqRecord) -> nx.DiGraph:
    dic_FT = {}
    try:
        with open(tf_set_file, 'r', encoding='utf-8') as f:
            lector = csv.reader(f, delimiter='\t'); next(lector, None)
            for fila in lector:
                if len(fila) > 4: dic_FT[fila[1]] = {'gene': fila[3], 'locus_tag': fila[4]}
    except Exception as e: 
        return nx.DiGraph()

    G = nx.DiGraph()
    try:
        with open(tf_riset_file, 'r', encoding='utf-8') as f:
            lector = csv.reader(f, delimiter='\t'); next(lector, None)
            for fila in lector:
                if len(fila) < 17: continue
                TG = fila[16]; FT = fila[3]
                if TG != '' and FT in dic_FT:
                    locus_tag_FT = dic_FT[FT]['locus_tag']
                    locus_tag_TG = gene_qualifier(TG, 'gene', 'locus_tag', genoma)[1]
                    if not locus_tag_TG: continue
                    
                    if not G.has_node(locus_tag_FT):
                        G.add_node(locus_tag_FT, ntype='TF', name=dic_FT[FT]['gene'])
                    
                    gens_opero = operon(locus_tag_TG, max_dist, genoma) if detectar_operons else [locus_tag_TG]
                    for OG in gens_opero:
                        if not OG: continue
                        nom_gen_OG = gene_qualifier(OG, 'locus_tag', 'gene', genoma)[1]
                        if not G.has_node(OG): G.add_node(OG, ntype='TG', name=nom_gen_OG)
                        G.add_edge(locus_tag_FT, OG)
    except: 
        pass
    return G