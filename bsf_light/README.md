# Beams Simple and Fast LIGHT simulation (bsf_light) - light simulation tool using a beam-spread-function approach

Replication and improvement of original model introduced by Yona et al. 2016 [1]. Replication details and model published in [2].

## software requirements:

* Requires python 3.13. Dependencies will be installed automatically. See bsf_light/setup.py for details.

Tested on Ubuntu 20.04 LTS and 22.04.4 LTS.

## Pip installation:

```pip install bsf-light```

## Use code for modeling light propagation in cortical tissue

1. Define simulation parameters, example can be found under `examples/default.yml` and under https://github.com/CSNG-MFF/bsf_light in `params/default.yml`

2. Use commandline to run simulation using the run-script provided in `examples/run.py` and under https://github.com/CSNG-MFF/bsf_light in `scripts/run.py`, providing the parameter-file defined in step 1. and the location where simulation output shall be written to:
    ```python run.py PARAMETER_FILE OUTPUT_LOCATION```

## References

[1] G. Yona, N. Meitav, I. Kahn, S. Shoham, Realistic Numerical and Analytical Modeling of Light Scattering in Brain Tissue for Optogenetic Applications. eNeuro 3, ENEURO.0059-15.2015 (2016).
[2] To be announced.

