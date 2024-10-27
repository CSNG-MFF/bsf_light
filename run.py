import sys
from BSF import calc_I_fiber
from load_save_utils import load_yaml, save_pickle
from time import time
start = time()
results = calc_I_fiber(load_yaml(str(sys.argv[1])))
end = time()
results['comp_time_s'] = end - start

save_pickle(results, str(sys.argv[2]))