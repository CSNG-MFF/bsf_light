import numpy as np
from scipy.special import gamma

def h(z, rho, tau, c):
    """
    Spatial-angular distribution of scattered photons.

    The result depends on the radial coordinate (`rho`), the z-coordinate (`z`), 
    the multipath time (`tau`), and the speed of light (`c`).
    The function diverges if rho=0 for tau->0.

    Parameters
    ----------
    z : float
        The z-coordinate, representing the distance along the depth.
    rho : float
        The radial coordinate, representing the distance from the center in the xy-plane.
    tau : float
        The multipath time, which affects how far scattered photons can travel.
    c : float
        The speed of light in the medium.

    Returns
    -------
    h: float

    """
    _3by4taucz = 3 / (4 * tau * c * z)
    return _3by4taucz / np.pi * np.exp(- _3by4taucz * rho*rho)

def time_dispersion_moments(z, g, mu_s, c, version: str):
    """
    Compute time dispersion moments thorugh diff. formulas as in McLean et al.

    The calculations for the 2nd moment of time dispersion assume that the 2nd moment of the
    scattering angle (`Theta^4`) is negligible.

    Parameters
    ----------
    z : float
        Distance along the z-axis, representing the propagation depth of the light.
    g : float
        Mean cosine of the scattering angle (as described by the anisotropy factor).
    mu_s : float
        The scattering coefficient, representing the number of scattering events per unit length.
    c : float
        The speed of light in the medium.
    version : str
        Formulas to use for calculation:
        - 'eq4': Uses the eq. 4 from McLean et al. for precise calculation of `mu` and approximates
          mu2_by_sigma2 according to approximations by Lutomirski et al. as in table 1 by McLean et
          al.
        - 'table1_Lutomirski': Uses approximations from Lutomirski et al. as in table 1 by McLean
        et al.
        - 'table1_vandeHulst': Uses approximations from van de Hulst and Kattawar as in table 1 by
        McLean et al.

    Returns
    -------
    mu : float
        1st moment of time dispersion

    mu2_by_sigma2 : float
        The ratio of the first and second moment of time dispersion squared.
    """
    if version == 'eq4':
        v = 1 - g # g = mean(cos(theta))
        bzv = mu_s * z * v # bz represents number of scattering lengths 
        mu = (z/c) * (1 - (1 - np.exp(-bzv)) / bzv)
        mu2_by_sigma2 = (1/4)**2 / (1/24)
    
    elif version == 'table1_Lutomirski':
        mean_Theta2 = 2 * (1 - g) # comment below table1
        mu = (z/c) * (1/4) * mu_s * z * mean_Theta2
        mu2_by_sigma2 = (1/4)**2 / (1/24)
    
    elif version == 'table1_vandeHulst':
        mean_Theta2 = 2 * (1 - g) # comment below table1
        mu = (z/c) * (1/12) * mu_s * z * mean_Theta2
        mu2_by_sigma2 = (1/12)**2 / (7/720)

    else:
        raise ValueError("version must be one out of ['eq4', 'table1_vandeHulst', 'table1_Lutomirski']")

    return mu, mu2_by_sigma2

def G(z, tau, g, mu_s, c, moments: str):
    """
    Time dispersion distribution.

    Diverges for tau -> 0 and not defined for z=0. Moments can be calculated in different
    ways, see ´time_dispersion_moments´ function.

    Parameters
    ----------
    z : float
        Distance along the z-axis, representing the propagation depth of the light.
    tau : float
        The multipath time, which affects how far scattered photons can travel.
    g : float
        Mean cosine of the scattering angle (as described by the anisotropy factor).
    mu_s : float
        The scattering coefficient, representing the number of scattering events per unit length.
    c : float
        The speed of light in the medium.
    moments : str
        Formulas to use for calculation of the moments of multipath time:
        - 'eq4': Uses the eq. 4 from McLean et al. for precise calculation of `mu` and approximates
          mu2_by_sigma2 according to approximations by Lutomirski et al. as in table 1 by McLean et
          al.
        - 'table1_Lutomirski': Uses approximations from Lutomirski et al. as in table 1 by McLean
        et al.
        - 'table1_vandeHulst': Uses approximations from van de Hulst and Kattawar as in table 1 by
        McLean et al.

    Returns
    -------
    G: float
        The time dispersion distribution value at the given parameters.
    """
    # moments of time dispersion
    mu, mu2_by_sigma2 = time_dispersion_moments(z, g, mu_s, c, version=moments)
    mu_by_sigma2 = mu2_by_sigma2 / mu
    mu_tau_by_sigma2 = mu_by_sigma2 * tau
    # G in three factors
    G1 = mu_by_sigma2/gamma(mu2_by_sigma2)
    G2 = mu_tau_by_sigma2**(mu2_by_sigma2 - 1)
    G3 = np.exp(-mu_tau_by_sigma2)
    return G1 * G2 * G3

def pencil_scattered(z, rho, tau, g, mu_s, mu_a, c, G_version: str):
    """
    Pencil beam of scattered photons before time integration.

    Parameters
    ----------
    z : float
        Distance along the z-axis, representing the propagation depth of the light.
    rho : float
        Distance from the center in the xy-plane.
    tau : float
        Multipath time, which affects how far scattered photons can travel.
    g : float
        Mean cosine of the scattering angle (as described by the anisotropy factor).
    mu_s : float
        Scattering coefficient, representing the number of scattering events per unit length.
    mu_a : float
        Absorption coefficient, representing the number of absorbed photons per unit length.
    c : float
        The speed of light in the medium.
    G_version: str
        Formulas to use for calculation of the moments of multipath time:
        - 'eq4': Uses the eq. 4 from McLean et al. for precise calculation of `mu` and approximates
          mu2_by_sigma2 according to approximations by Lutomirski et al. as in table 1 by McLean et
          al.
        - 'table1_Lutomirski': Uses approximations from Lutomirski et al. as in table 1 by McLean
        et al.
        - 'table1_vandeHulst': Uses approximations from van de Hulst and Kattawar as in table 1 by
        McLean et al.

    Returns
    -------
    pencil_scattered: float
        Scattered pencil beam component before time integration.
    """
    scattered = 1 - np.exp(-mu_s * z)
    not_absorbed = np.exp(-mu_a * (z + c * tau))
    return scattered*not_absorbed*G(z, tau, g, mu_s, c, G_version)*h(z, rho, tau, c)

def pencil_scattered_time_integrated(z, rho, tau, g, mu_s, mu_a, c, G_version: str):
    """
    Pencil beam of scattered photons (after integration over multipath time).

    Important: Function assumes multipath time (tau) to be on axis/dim 2 in arrays.

    Parameters
    ----------
    z : array_like, tau on dim 2
        Distance along the z-axis, representing the propagation depth of the light.
        Note functions assumes multipath time to be on dim 2 of the array.
    rho : array_like, tau on dim 2
        Distance from the center in the xy-plane.
        Note functions assumes multipath time to be on dim 2 of the array.
    tau : array_like, tau on dim 2
        Multipath time, which affects how far scattered photons can travel.
        Note functions assumes multipath time to be on dim 2 of the array.
    g : float
        Mean cosine of the scattering angle (as described by the anisotropy factor).
    mu_s : float
        Scattering coefficient, representing the number of scattering events per unit length.
    mu_a : float
        Absorption coefficient, representing the number of absorbed photons per unit length.
    c : float
        The speed of light in the medium.
    G_version: str
        Formulas to use for calculation of the moments of multipath time:
        - 'eq4': Uses the eq. 4 from McLean et al. for precise calculation of `mu` and approximates
          mu2_by_sigma2 according to approximations by Lutomirski et al. as in table 1 by McLean et
          al.
        - 'table1_Lutomirski': Uses approximations from Lutomirski et al. as in table 1 by McLean
        et al.
        - 'table1_vandeHulst': Uses approximations from van de Hulst and Kattawar as in table 1 by
        McLean et al.

    Returns
    -------
    pencil_scattered_time_integrated: float
        Scattered pencil beam component after time integration.
    """
    integral = np.sum(
        pencil_scattered(
            z[:,:,:-1], rho[:,:,:-1], tau[:,:,:-1], g, mu_s, mu_a, c, G_version
        ) * np.diff(tau, axis=2),
        axis=2
    )
    return integral
