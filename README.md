# Beams Simple and Fast LIGHT simulation (bsf_light) - light simulation tool using a beam-spread-function approach [![DOI](https://zenodo.org/badge/859378712.svg)](https://doi.org/10.5281/zenodo.15490385)

Replication and improvement of original model introduced by Yona et al. 2016 [1]. Replication details and model published in [2].

## software requirements:

* Requires python 3.13, dependencies will be installed automatically. See /bsf_light/setup.py file for more details. 

Tested on Ubuntu 20.04 LTS and 22.04.4 LTS.

## Pip installation:

* ```pip install bsf-light```

* install additional dependency: ```pip install pandas```

## Run simulations and reproduce figures:

* Run in bash console: ```bash run_all.sh```

## Use code for modeling light propagation in cortical tissue

1. Define simulation parameters, example can be found under `bsf_light/examples/default.yml` and under https://github.com/CSNG-MFF/bsf_light in `params/default.yml`

2. Use commandline to run simulation using the run-script provided in `bsf_light/examples/run.py` and under https://github.com/CSNG-MFF/bsf_light in `scripts/run.py`, providing the parameter-file defined in step 1. and the location where simulation output shall be written to: ```python scripts/run.py PARAMETER_FILE OUTPUT_LOCATION```, e.g., ```python scripts/run.py params/defaults.yml result.pickle``` will run the simulation using default parameters and write results to result.pickle.

## References

[1] G. Yona, N. Meitav, I. Kahn, S. Shoham, Realistic Numerical and Analytical Modeling of Light Scattering in Brain Tissue for Optogenetic Applications. eNeuro 3, ENEURO.0059-15.2015 (2016).

[2] D. Berling, J. Střeleček, T. Iser, J. Antolík, An open-source replication for fast and accessible light propagation modeling in brain tissue. In revision at PLOS ONE. 2025.
