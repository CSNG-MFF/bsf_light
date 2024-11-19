## Beam-spread-function-based light model

Replication and improvement of original model introduced by Yona et al. 2016 [1]. Replication details and model published in [2].

## installation using pip, python 3.x:

* ´´´pip install ./BSF´´´

* install additional dependency: ´´´pip install pandas´´´

## Run simulations and reproduce figures:

* Run in bash console: ´´´bash run_all.sh´´´

## Use code for modeling light propagation in cortical tissue

* ´BSF/BSF.py´ contains main model implementation. 

* ´BSF/utils.py´, ´BSF/load_save_utils.py´ contain utilitly functions

* ´scripts/load_original.py´ contains functions to load data from Yona et al. [1]

* To run simulation with custom parameters, adapt parametrization in ´params/default.yml´ and save as new file, e.g., ´params/custom.yml´ and run: 
´´´python scripts/run.py params/custom.yml results/custom.pickle´´´


## References
[1] G. Yona, N. Meitav, I. Kahn, S. Shoham, Realistic Numerical and Analytical Modeling of Light Scattering in Brain Tissue for Optogenetic Applications. eNeuro 3, ENEURO.0059-15.2015 (2016).
[2] To be announced.
