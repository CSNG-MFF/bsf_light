import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from clean_BSF import calc_I_fiber
import pickle
from clean_utils import type_cast_paramdict

results = calc_I_fiber(type_cast_paramdict(snakemake.params.simulation))

with open(str(snakemake.output), 'wb') as handle:
    pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)
