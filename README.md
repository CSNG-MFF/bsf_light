# Beam-spread-function-based light model implemented according to Yona et al. 2016 [1]

## (Known) software requirements:

* python 3.x
* numpy
* pandas
* matplotlib
* scipy
* tqdm

Tested on Ubuntu 20.04 LTS and 22.04.4 LTS.

## Conda installation:

'''conda create -n BSF -c conda-forge numpy pandas matplotlib scipy tqdm'''

## Reproduce figures

'''bash run_all.sh'''
'''bash run_all_plots.sh'''

## Use code for modeling light propagation in cortical tissue

* 'BSF.py' contains main model implementation. 

* 'utils.py', 'load_original.py', 'load_save_utils.py' contain utilitly functions

* To run simulation with custom parameters, adapt parametrization in params/default.yml and save as under new file, e.g., params/custom.yml and run: '''python run.py params/custom.yml results/custom.pickle'''


## References
[1] G. Yona, N. Meitav, I. Kahn, S. Shoham, Realistic Numerical and Analytical Modeling of Light Scattering in Brain Tissue for Optogenetic Applications. eNeuro 3, ENEURO.0059-15.2015 (2016).
