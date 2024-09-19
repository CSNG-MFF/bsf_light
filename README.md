# Beam-spread-function-based light model implemented according to Yona et al. 2016 [1]

## (Known) software requirements:

* python 3.x
* numpy
* pandas
* snakemake >5.13
* matplotlib
* scipy
* tqdm

Tested on Ubuntu 20.04 LTS and 22.04.4 LTS.

## Conda installation:

'''conda create -n BSF -c conda-forge -c bioconda numpy pandas snakemake">=5.13" matplotlib scipy tqdm'''

## Reproduce figures

'''bash run_all.sh'''

## Use code for modeling light propagation in cortical tissue

* 'BSF.py' contains main model implementation. 

* 'utils.py' contains utilitly functions

* See example on how to use code for light simulations in Figure1.py.


## References
[1] G. Yona, N. Meitav, I. Kahn, S. Shoham, Realistic Numerical and Analytical Modeling of Light Scattering in Brain Tissue for Optogenetic Applications. eNeuro 3, ENEURO.0059-15.2015 (2016).
