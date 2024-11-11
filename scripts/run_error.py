import sys
from BSF.BSF import calc_I_fiber_reproduce_error
from BSF import load_yaml, save_pickle
from time import time
start = time()
results = calc_I_fiber_reproduce_error(
        load_yaml(str(sys.argv[1])),
        fix_pencil_rho_z=(500,750),
        ang_conv_fill_nan=True,
        decrease_disk_conv_volume_radially=100
)
end = time()
results['comp_time_s'] = end - start

save_pickle(results, str(sys.argv[2]))
