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

import dill

from graph_tool.all import *
import pandas as pd
import numpy as np

if __name__ == '__main__':

    logging.info("Reading trimmed graph...")
    g = load_graph(snakemake.input[0])
    
    logging.info("Creating nested block model...")
    if snakemake.params.type == "layer": 
        state_min = minimize_nested_blockmodel_dl(g, init_bs=None,
                                                  state_args=dict(base_type=LayeredBlockState,
                                                   ec=g.ep.layer, layers=True,
                                                   recs=[g.ep.z_s], 
                                                   rec_types=["real-normal"]))
    elif snakemake.params.type == "batch":
        state_min = minimize_nested_blockmodel_dl(g, init_bs=None,
                                                  state_args=dict
                                                  (recs=[g.ep.z_s],
                                                   rec_types=["real-normal"]))
    else: 
        sys.exit('Parameter "type" is not batch or layer')
                                                    
    logging.info("Saving blockstate...")
    block_state = state_min.get_bs()
    with open(snakemake.output[0], 'wb') as fh:
        dill.dump(block_state, fh, recurse=True)
    logging.info("Done!")
