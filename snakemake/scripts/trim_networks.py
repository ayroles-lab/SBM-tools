import os
os.environ["OMP_NUM_THREADS"] = "32"

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
from sklearn.covariance import LedoitWolf, OAS
import statsmodels.api as sm

from multipy.fdr import lsu

def filterByEdge(g, corr, cutOff, keepOnlyMain):
    # Filtering edges
    corr = g.edge_properties[corr]
    sign = g.new_ep("bool", True)
    sign.a = np.array(np.abs(corr.a) > cutOff)

    tv = GraphView(g, efilt=sign)

    # Keeping largest component
    if keepOnlyMain:
        comp, hist = label_components(tv)
        main_component = tv.new_vp("bool", (comp.a == np.where(hist == max(hist))[0][0]))
        tv.vertex_properties["main_component"] = main_component
        tv.set_vertex_filter(main_component)
    return tv

def filterByFDR(g, level, keepOnlyMain):
    # Filtering edges
    pvals = np.array(g.edge_properties["pvalue"].a)

    fdr_ep = g.new_ep("bool", True)
    fdr_ep.a = lsu(pvals, q=level)

    tv = GraphView(g, efilt=fdr_ep)

    # Keeping largest component
    if keepOnlyMain:
        comp, hist = label_components(tv)
        main_component = tv.new_vp("bool", (comp.a == np.where(hist == max(hist))[0][0]))
        tv.vertex_properties["main_component"] = main_component
        tv.set_vertex_filter(main_component)
    return tv

def saveGeneTable(g, gene_expr, output=None, transpose = False):
    gene_list = []
    for i in g.vertex_properties['genes']:
        gene_list.append(i)
    expr_df = gene_expr[gene_list]
    if transpose is True:
        expr_df = expr_df.T
    if output is not None:
        expr_df.to_csv(output)
    return expr_df

if __name__ == '__main__':

    logging.info("FDR level:" + str(snakemake.wildcards.fdr))
    fdr = float(snakemake.wildcards.fdr)

    logging.info("Reading full graph...")
    g = load_graph(snakemake.input[0])

    logging.info("Trimming...")

    tv = filterByFDR(g, fdr, True)
    gi = Graph(tv, prune = True)

    N = len(gi.get_vertices())
    logging.info(str(N) + " genes")
    Et = (N * N - N)/2
    E = len(gi.get_edges())
    density = E/Et
    logging.info("Density: " + str(density))
    logging.info("Min correlation: " + str(min(gi.edge_properties["spearman"].a)))
    logging.info("Max correlation: " + str(max(gi.edge_properties["spearman"].a)))
    
    logging.info("Adding z-transformed weights...")
    spearman = gi.edge_properties["spearman"]
    gi.edge_properties["z_s"] = gi.new_edge_property("double", (2*np.arctanh(spearman.a)))

    logging.info("Writting trimmed graph...")
    gi.save(snakemake.output[0])
    logging.info("Done!")
