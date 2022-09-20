import os
os.environ["OMP_NUM_THREADS"] = '16' # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = '16' # export OPENBLAS_NUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = '16' # export NUMEXPR_NUM_THREADS=6

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

import dill

from graph_tool.all import *
import numpy as np
import pandas as pd

def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)
    if not os.path.exists(file_path):
        os.makedirs(file_path)

def get_group(x, state):
    levels = state.get_levels()
    n_levels = 7#len(levels)
    r = np.zeros(n_levels)
    r[0] = levels[0].get_blocks()[x]
    for i in range(1, n_levels):
        r[i] = levels[i].get_blocks()[r[i-1]]
    r = r.astype(int)
    return r

def create_nestedBlock_df(g, corr, state):
    genes = g.vertex_properties["genes"]
    nested_block_df = pd.DataFrame(columns=('Gene', "Degree", "E_corr", 'B1', "B2", "B3", "B4", "B5", "B6", "B7"))
    for v in g.vertex_index:
        line = [genes[v]]
        line.append(g.get_total_degrees([v])[0])
        line.append(np.mean(np.abs(g.get_all_edges(v, [corr] )[:,2])))
        [line.append(i) for i in get_group(v, state)]
        nested_block_df.loc[v] = line
    return nested_block_df

if __name__ == '__main__':

    logging.info("Reading trimmed graph...")
    g = load_graph(snakemake.input.graph)

    logging.info("Arctan transform on correlations...")
    corr = g.edge_properties["spearman"]
    g.ep.z_s = g.new_edge_property("double", (2*np.arctanh(corr.a)))
    
    logging.info("Loading blocks...")
    with open (snakemake.input.blocks, "rb") as fh:
        bs = dill.load(fh)[0:7]

    logging.info("Creating nested block model...")
    state = minimize_nested_blockmodel_dl(g, init_bs=bs,
                                              state_args=dict(recs=[g.ep.z_s],
                                              rec_types=["real-normal"]))

    logging.info("Creating block DataFrame...")
    block_df = create_nestedBlock_df(g, corr, state)

    logging.info("Calculating block sizes...")
    blocks = [list(set(block_df[b])) for b in block_df.filter(like='B', axis=1)]
    block_sizes = [len(b) for b in blocks]
    block_sizes = -np.sort(-np.array(list(set(block_sizes))))
    block_sizes = [x for x in block_sizes if x >= 2]

    block_df["Gene"].to_csv(snakemake.params.blockDir + "/background.csv", header=False, index=False)

    logging.info("Creating gene lists...")
    n_levels = len(block_sizes)
    output_df = pd.DataFrame(columns=('Nested_Level', 'Block', 'File', 'N_genes', 'Internal_degree', 'Average_internal_degree', 'Total_degree', 'Average_total_degree' , 'Assortativity'))
    l = 0
    for i in range(n_levels):
        logging.info("At level: " + str(i+1))
        bl = blocks[i]
        for b in bl:

            line = [i+1]
            line.append(b)

            df = block_df[block_df['B' + str(i+1)]==b]
            genes = df["Gene"]
            file_name = "/" + '-'.join([str(num) for num in list(df.filter(like='B', axis=1).iloc[0,range(i, n_levels)])]) + ".csv"

            line.append(file_name)
            N_genes = genes.shape[0]
            line.append(N_genes)

            ensure_dir(snakemake.params.blockDir + "/Level_" + str(i+1))
            genes.to_csv(snakemake.params.blockDir + "/Level_" + str(i+1) + file_name, header=False, index=False )

            # Weighted
            ers = adjacency(state.levels[i].bg, weight=state.levels[i].mrs)
            B = len(bl)
            E = ers.sum()
            q_r = (B/E) * (ers[b,b] - (ers[b,:].sum()**2/E))
            if np.abs(q_r) > 1:
                sys.error("q_r is larger than one.")

            line.append(ers[b,b])
            line.append(ers[b,b] / N_genes)
            line.append(ers[b,:].sum())
            line.append(ers[b,:].sum() / N_genes)
            line.append(q_r)

            output_df.loc[l] = line
            l = l + 1

    logging.info("Outputing block summary...")
    output_df.to_csv(snakemake.output.blockSummary, index=False )

    logging.info("Outputing block DataFrame...")
    block_df.to_csv(snakemake.output.blockDF, index=False )

    logging.info("Done!")



