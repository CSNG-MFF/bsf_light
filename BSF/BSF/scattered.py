import numpy as np
from scipy.special import gamma

def h(z, rho, tau, c):
    """
    spatial-angular distribution of scattered photons.
    Diverges if rho=0 for tau->0.
    """
    _3by4taucz = 3 / (4 * tau * c * z)
    return _3by4taucz / np.pi * np.exp(- _3by4taucz * rho*rho)

def first_moment_of_time_dispersion_mu_tau(
    z, g, mu_s, c,
    version: str
):
    if version == 'eq4':
        # use eq. 4 in McLean et al. to calculate mu precisely
        v = 1 - g # g = mean(cos(theta))
        bzv = mu_s * z * v # bz represents number of scattering lengths 
        mu = (z/c) * (1 - (1 - np.exp(-bzv)) / bzv)
        # use approximations by Lutomirski etc. in table1 to relate mu and sigma,
        # further assume mean(Theta**4) irrelevant
        mu2_by_sigma2 = (1/4)**2 / (1/24)
    
    elif version == 'table1_Lutomirski':
        # mu_tau in table 1 in McLean et al., Dolin, Ishimaru, Lutomirski et al.
        mean_Theta2 = 2 * (1 - g) # comment below table1
        mu = (z/c) * (1/4) * mu_s * z * mean_Theta2
        # use approximations by Lutomirski etc. in table1 to relate mu and sigma,
        # further assume mean(Theta**4) irrelevant
        mu2_by_sigma2 = (1/4)**2 / (1/24)
    
    elif version == 'table1_vandeHulst':
        # mu_tau in table 1 in McLean et al., van de Hulst and Kattawar
        mean_Theta2 = 2 * (1 - g) # comment below table1
        mu = (z/c) * (1/12) * mu_s * z * mean_Theta2
        mu2_by_sigma2 = (1/12)**2 / (7/720)
    
    else:
        raise ValueError("version must be one out of ['eq4', 'table1_vandeHulst', 'table1_Lutomirski']")
    
    return mu, mu2_by_sigma2
    

def G(z, tau, g, mu_s, c, first_moment: str):
    """
    Time dispersion distribution.
    Diverges for tau -> 0 and not defined for z=0.
    
    Moments can be calculated via 4 ways, which change the way the
    first moment (mu_tau) is calculated:
    'eq4', 'table1_vandeHulst', or 'table1_Lutomirski'
    """
    # moments of time dispersion
    mu, mu2_by_sigma2 = first_moment_of_time_dispersion_mu_tau(z, g, mu_s, c, version=first_moment)
    mu_by_sigma2 = mu2_by_sigma2 / mu
    mu_tau_by_sigma2 = mu_by_sigma2 * tau
    # G in three factors
    G1 = mu_by_sigma2/gamma(mu2_by_sigma2)
    G2 = mu_tau_by_sigma2**(mu2_by_sigma2 - 1)
    G3 = np.exp(-mu_tau_by_sigma2)
    return G1 * G2 * G3

def pencil_scattered(z, rho, tau, g, mu_s, mu_a, c, G_version: str):
    scattered = 1 - np.exp(-mu_s * z)
    not_absorbed = np.exp(-mu_a * (z + c * tau))
    return scattered*not_absorbed*G(z, tau, g, mu_s, c, G_version)*h(z, rho, tau, c)

def pencil_scattered_time_integrated(z, rho, tau, g, mu_s, mu_a, c, G_version: str):
    """
    Assume tau on axis/dim 2 in arrays.
    """
    integral = np.sum(
        pencil_scattered(
            z[:,:,:-1], rho[:,:,:-1], tau[:,:,:-1], g, mu_s, mu_a, c, G_version
        ) * np.diff(tau, axis=2),
        axis=2
    )
    return integral

