import os
os.environ["OMP_NUM_THREADS"] = '8' # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = '8' # export OPENBLAS_NUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = '8' # export NUMEXPR_NUM_THREADS=6

import sys,os
import logging, traceback
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    )
def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(''.join(["Uncaught exception: ",
                         *traceback.format_exception(exc_type, exc_value, exc_traceback)
                         ])
                 )
# Install exception handler
sys.excepthook = handle_exception

from graph_tool.all import *
import pandas as pd
import numpy as np
import scipy as sp
import statsmodels.api as sm

if __name__ == '__main__':

    if snakemake.params.type == "layer": 
        logging.info("Building layered network...")
        data_dir = snakemake.input[0]
        input_list = []
        for file in os.listdir(data_dir):
            if file.endswith(".csv"):
                input_list.append(file)
        input_names = list(map(lambda p: p[:p.rfind('.')], input_list))
        gene_expr = []
        for file in input_list:
            gene_expr_raw = pd.read_table(os.path.join(data_dir, file))
            gene_expr.append(gene_expr_raw.T)
    elif snakemake.params.type == "batch":
        logging.info("Building network...")
        gene_expr_raw = pd.read_table(snakemake.input[0])
        gene_expr = [gene_expr_raw.T]
    else: 
        sys.exit('Parameter "type" is not batch or layer')
    logging.info("Read expression files...")


    n_layers = len(gene_expr)

    X_centered = []
    logging.info("Centering data matrix...")
    for ge in gene_expr:
        aux = (ge - ge.mean()) / np.sqrt(ge.var())
        X_centered.append(aux)

    gene_names = []
    number_of_genes = ""
    for ge in gene_expr:
        gene_names = set().union(gene_names, ge.columns)
        number_of_genes = number_of_genes + ", " + str(len(ge.columns)) + " genes"
    gene_names = list(set(gene_names))
    logging.info("Concatenating gene names" + number_of_genes)

    
    n_genes = len(gene_names)
    logging.info("Total number of genes: " + str(n_genes))

    logging.info("Creating graph...")
    g = Graph(directed=False)
    g.add_vertex(n = n_genes)
    genes = g.new_vertex_property("string", np.array(gene_names, dtype = "str"))
    g.vertex_properties["genes"] = genes

    spearman = g.new_ep("double", 0)
    pval = g.new_ep("double", 0)
    layer = g.new_ep("int", 0)
    dataset = g.new_ep("string", "")

    unit_corr = dict()
    for k in range(n_layers):
        logging.info("Starting layer " + str(k) + " - " + str(input_names[k]) + "...")
        ge = gene_expr[k]
        ge_scale = X_centered[k]
        n_gene_layer = ge.shape[1]
        for i in range(n_gene_layer):
            for j in range(i):
                v_i = find_vertex(g, g.vp.genes, ge.columns[i])[0]
                v_j = find_vertex(g, g.vp.genes, ge.columns[j])[0]
                spearman_r = sp.stats.spearmanr(ge_scale.iloc[:,i], ge_scale.iloc[:,j])
                e = g.add_edge(v_i, v_j)
                if snakemake.params.type == "layer":
                    layer[e] = k
                    dataset[e] = input_names[k]
                pval[e] = spearman_r[1]
                spearman[e] = spearman_r[0]
                if(abs(spearman_r[0]) > 0.999): 
                    key_i = str(g.vp.genes[v_i])
                    key_j = str(g.vp.genes[v_j])
                    logging.info("Warning: Spearman correlation between " + 
                                 key_i + " and "  + key_j + " is one. Dropping one of them.")
                    if not ((key_i in unit_corr) or (key_j in unit_corr)):
                        unit_corr[str(g.vp.genes[v_i])] = v_i
        logging.info("Finished layer " + str(k) + ".")

    logging.info("Setting edge properties...")
    g.edge_properties["pvalue"] = pval
    g.edge_properties["spearman"] = spearman    
    
    if snakemake.params.type == "layer":
        logging.info("Setting layer properties...")
        g.edge_properties["layer"] = layer
        g.edge_properties["dataset"] = dataset
    
    for gene in unit_corr.keys():
        logging.info("Removing vertex " + str(gene) + "!")
    g.remove_vertex(unit_corr.values()) 
    
    logging.info("Setting edge weights...")
    g.edge_properties["z_s"] = g.new_edge_property("double", (2*np.arctanh(g.ep.spearman.a)))
    
    logging.info("Writting graph...")
    g.save(snakemake.output[0])
    logging.info("Done!")
