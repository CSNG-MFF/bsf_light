import pytest
import numpy as np
from bsf_light import calc_I_fiber, load_yaml, load_pickle

def test_simulation_default_params():
    # Run the simulation with default parameters
    result = calc_I_fiber(load_yaml('tests/default_params.yml'))['final']['combined']
    
    # Load the default result (pre-saved output for default params)
    default_result = load_pickle("tests/default_result.pickle")['final']['combined']
    
    # Assert result matches within a tolerance
    np.testing.assert_allclose(result, default_result, rtol=1e-5, atol=1e-8)

