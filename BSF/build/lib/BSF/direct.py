import numpy as np

def I_direct_cone(z, rho, params):
    """ 
    Compute the intensity of a direct light cone emitted from an infinitesimal point.

    This function calculates the intensity of direct light emitted from an infinitesimal
    point source by performing an analytical angular convolution over the direct light
    component.

    Parameters
    ----------
    rho : np.ndarray
        Radial coordinate, representing the distance from the center along the
        xy-plane.
    z : np.ndarray
        z-coordinate, representing the height along the z-axis.
    params : dict
        Parameters used in the calculation, refer to the `calc_I_fiber` function in the
        `fiber.py` module.

    Returns
    -------
    intensity : np.ndarray
        Computed intensity values for given `rho` and `z`.
    """
    assert np.all(rho >=0) and np.all(z >= 0), 'rho or z is negative.'
    on = np.arctan(rho/z) <= params['theta_div']
    R_spherical = np.sqrt(rho*rho+z*z)
    norm = (1 - np.cos(params['theta_div'])) * 2 * np.pi * R_spherical*R_spherical
    return np.exp(-(params['mu_a']+params['mu_s'])*R_spherical) * on / norm

