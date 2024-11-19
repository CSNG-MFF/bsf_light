# Beam-spread-function-based light model

Replication and improvement of original model introduced by Yona et al. 2016 [1]. Replication details and model published in [2].

## software requirements:

* python 3.x
* numpy
* matplotlib
* scipy
* pyyaml

Tested on Ubuntu 20.04 LTS and 22.04.4 LTS.

## Pip installation:

'''pip install PATH_TO_BSF'''

## Use code for modeling light propagation in cortical tissue

1. Define simulation parameters, example can be found under ´´´examples/default.yml´´´ and under https://github.com/CSNG-MFF/BSF_light_model in ´´´params/default.yml´´´

2. Use commandline to run simulation using the run-script provided in ´´´examples/run.py´´´ and under https://github.com/CSNG-MFF/BSF_light_model in ´´´scripts/run.py´´´, providing the parameter-file defined in step 1. and the location where simulation output shall be written to:
	´´´python run.py PARAMETER_FILE OUTPUT_LOCATION´´´

## References

[1] G. Yona, N. Meitav, I. Kahn, S. Shoham, Realistic Numerical and Analytical Modeling of Light Scattering in Brain Tissue for Optogenetic Applications. eNeuro 3, ENEURO.0059-15.2015 (2016).
[2] To be announced.
