import os
os.environ["OMP_NUM_THREADS"] = "8"

import sys

from graph_tool.all import *
import pandas as pd
import numpy as np
import scipy as sp
import statsmodels.api as sm

if __name__ == '__main__':

    if snakemake.params.type == "layer": 
        data_dir = snakemake.input[0]
        input_list = []
        for file in os.listdir(data_dir):
            if file.endswith(".tsv"):
                input_list.append(file)
        input_names = list(map(lambda p: p[:p.rfind('.')], input_list))
        gene_expr = []
        for file in input_list:
            gene_expr_raw = pd.read_table(os.path.join(data_dir, file))
            gene_expr.append(gene_expr_raw.T)
    elif snakemake.params.type == "batch":
        gene_expr_raw = pd.read_table(snakemake.input[0])
        gene_expr = [gene_expr_raw.T]
    else: 
        sys.exit('Parameter "type" is not batch or layer')

    n_layers = len(gene_expr)

    X_centered = []
    for ge in gene_expr:
        aux = (ge - ge.mean()) / np.sqrt(ge.var())
        X_centered.append(aux)

    gene_names = []
    for ge in gene_expr:
        gene_names = set().union(gene_names, ge.columns)
    gene_names = list(set(gene_names))
    
    n_genes = len(gene_names)
    g = Graph(directed=False)
    g.add_vertex(n = n_genes)
    genes = g.new_vertex_property("string", np.array(gene_names, dtype = "str"))
    g.vertex_properties["genes"] = genes

    spearman = g.new_ep("double", 0)
    pval = g.new_ep("double", 0)
    layer = g.new_ep("int", 0)
    dataset = g.new_ep("string", "")

    for k in range(n_layers):
        ge = gene_expr[k]
        ge_scale = X_centered[k]
        n_gene_layer = ge.shape[1]
        for i in range(n_gene_layer):
            for j in range(i):
                v_i = gene_names.index(ge.columns[i])
                v_j = gene_names.index(ge.columns[j])
                spearman_r = sp.stats.spearmanr(ge_scale.iloc[:,i], ge_scale.iloc[:,j])
                e = g.add_edge(v_i, v_j)
                if snakemake.params.type == "layer":
                    layer[e] = k
                    dataset[e] = input_names[k]
                pval[e] = spearman_r[1]
                spearman[e] = spearman_r[0]

    g.edge_properties["pvalue"] = pval
    g.edge_properties["spearman"] = spearman
    g.edge_properties["z_s"] = g.new_edge_property("double", (2*np.arctanh(spearman.a)))
    if snakemake.params.type == "layer":
        g.edge_properties["layer"] = layer
        g.edge_properties["dataset"] = dataset

    g.save(snakemake.output[0])