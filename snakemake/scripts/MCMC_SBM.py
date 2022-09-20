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

    bs = []
    h = [np.zeros(g.num_vertices() + 1) for s in state_min.get_levels()]

    def collect_partitions(s):
        bs.append(s.get_bs())
        for l, sl in enumerate(s.get_levels()):
            B = sl.get_nonempty_B()
            h[l][B] += 1

    logging.info("Starting MCMC...")
    # Now we collect the marginals for exactly 1000*10 sweeps
    S1 = state_min.entropy()
    mcmc_equilibrate(state_min, force_niter=1000, mcmc_args=dict(niter=10),
                        callback=collect_partitions, verbose=True)

    pmode = PartitionModeState(bs, nested=True, converge=True)
    pv = pmode.get_marginal(g)

    # Get consensus estimate
    bs_max = pmode.get_max_nested()

    state = state_min.copy(bs=bs_max)
    S2 = state.entropy()


    logging.info("Description lenght improvement in MCMC: " + str(S2 - S1))
    logging.info("Final entropy after MCMC: " + str(state.entropy()))

    logging.info("Saving blockstate...")
    block_state = state.get_bs()
    with open(snakemake.output.blocks, 'wb') as fh:
        dill.dump(block_state, fh, recurse=True)

    logging.info("Saving block number histogram...")
    with open(snakemake.output.hist, 'wb') as fh:
        dill.dump(h, fh, recurse=True)

    logging.info("Done!")
