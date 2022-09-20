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
import numpy as np

if __name__ == '__main__':

    logging.info("Reading trimmed graph...")
    g = load_graph(snakemake.input.graph)

    logging.info("Arctan transform on correlations...")
    corr = g.edge_properties["spearman"]
    g.ep.z_s = g.new_edge_property("double", (2*np.arctanh(corr.a)))
    
    logging.info("Loading blocks...")
    with open (snakemake.input.blocks, "rb") as fh:
        bs = dill.load(fh)

    logging.info("Creating nested block model...")
    state_min = minimize_nested_blockmodel_dl(g, init_bs=bs,
                                              state_args=dict(recs=[g.ep.z_s],
                                              rec_types=["real-normal"]))

    logging.info("Starting annealing...")
    S1 = state_min.entropy()
    mcmc_anneal(state_min, beta_range=(1, 10), niter=1000,
                mcmc_equilibrate_args=dict(force_niter=10), verbose=False)
    S2 = state_min.entropy()
    logging.info("Improvement from annealing: " + str(S2 - S1))
    logging.info("Final entropy after annealing: " + str(state_min.entropy()))

    logging.info("Saving blockstate...")
    block_state = state_min.get_bs()
    with open(snakemake.output.blocks, 'wb') as fh:
        dill.dump(block_state, fh, recurse=True)
    logging.info("Done!")