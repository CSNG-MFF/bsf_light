import resource
import sys
from bsf_light import calc_I_fiber, load_yaml, save_pickle
from time import time
start = time()
results = calc_I_fiber(load_yaml(str(sys.argv[1])))
end = time()
results['comp_time_s'] = end - start
results['MaxRS_KB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

save_pickle(results, str(sys.argv[2]))
