import numpy as np

def I_direct_cone(z, rho, params):
    """ 
    Intensity of direct light cone exiting infinitesimal point.
    
    Calculated analytically as angular convolution of the expression
    for direct light pencil beam component.
    
    params:
    -------
    rho: np.ndarray
        radial coordinate
    z: np.ndarray
        z-coordinate
    params: dict
        For values, see I_fiber.
    """
    assert np.all(rho >=0) and np.all(z >= 0), 'rho or z is negative.'
    on = np.arctan(rho/z) <= params['theta_div']
    R_spherical = np.sqrt(rho*rho+z*z)
    norm = (1 - np.cos(params['theta_div'])) * 2 * np.pi * R_spherical*R_spherical
    return np.exp(-(params['mu_a']+params['mu_s'])*R_spherical) * on / norm

