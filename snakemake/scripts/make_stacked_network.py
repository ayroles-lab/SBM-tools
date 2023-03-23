import os
os.environ["OMP_NUM_THREADS"] = "8"

from graph_tool.all import *
import pandas as pd
import numpy as np
import scipy as sp
import statsmodels.api as sm

if __name__ == '__main__':

    gene_expr_raw = pd.read_table(snakemake.input[0])
    gene_expr = gene_expr_raw.T

    X_centered = (gene_expr - gene_expr.mean()) / np.sqrt(gene_expr.var())

    n_genes = gene_expr.shape[1]
    g = Graph(directed=False)
    g.add_vertex(n = n_genes)

    spearman = g.new_ep("double", 0)
    pval = g.new_ep("double", 0)
    genes = g.new_vertex_property("string", np.array(np.array(gene_expr.columns, dtype = "str")))
    g.vertex_properties["genes"] = genes

    for i in range(n_genes):
        for j in range(i):
            spearman_r = sp.stats.spearmanr(X_centered.iloc[:,i], X_centered.iloc[:,j])
            g.add_edge(i, j)
            e = g.edge(i, j)
            pval[e] = spearman_r[1]
            spearman[e] = spearman_r[0]

    g.edge_properties["pvalue"] = pval
    g.edge_properties["spearman"] = spearman

    g.save(snakemake.output[0])