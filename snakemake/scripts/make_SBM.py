import os
os.environ["OMP_NUM_THREADS"] = "8"

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
import pandas as pd
import numpy as np

if __name__ == '__main__':

    logging.info("Reading trimmed graph...")
    g = load_graph(snakemake.input[0])

    logging.info("Arctan transform on correlations...")
    corr = g.edge_properties["spearman"]
    g.ep.z_s = g.new_edge_property("double", (2*np.arctanh(corr.a)))
    
    logging.info("Creating nested block model...")
    #g = GraphView(g, vfilt=lambda v: np.random.uniform(0, 1) > 0.99)
    state_min = minimize_nested_blockmodel_dl(g, init_bs=None,
                                              state_args=dict(recs=[g.ep.z_s],
                                              rec_types=["real-normal"]))
                                                  
    logging.info("Saving blockstate...")
    block_state = state_min.get_bs()
    with open(snakemake.output[0], 'wb') as fh:
        dill.dump(block_state, fh, recurse=True)
    logging.info("Done!")
