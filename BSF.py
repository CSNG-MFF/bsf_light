from scipy.special import gamma
import numpy as np
from utils import rotate_cyl_coords_2angles_return_rho_z
from utils import log_smplng, calc_pencil_rho_z_max
from utils import Interpolator, calc_dependent_params, disk_conv_numpy

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
        # assume approximations by Lutomirski etc. in table1 for mu and sigma,
        # further assume mean(Theta**4) irrelevant
        mu2_by_sigma2 = (1/4)**2 / (1/24)
    
    elif version == 'table1_Lutomirski':
        # mu_tau in table 1 in McLean et al., Dolin, Ishimaru, Lutomirski et al.
        mean_Theta2 = 2 * (1 - g) # comment below table1
        mu = (z/c) * (1/4) * mu_s * z * mean_Theta2
        # assume approximations by Lutomirski etc. in table1 for mu and sigma,
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

def ang_conv(rho, z, func_rho_z, params):
    """
    Numerical angular convolution to obtain scattered light of
    light cone exiting infinitesimal point in emitter surface.

    Numerical convolution over angles theta
    and phi of function depending on rho and z.
    """
    # uniform sampling
    thetas = np.linspace(0, params['theta_div'], params['nstepstheta'])
    dtheta = np.diff(thetas)[0]
    phis = np.arange(0, 2*np.pi, 2*np.pi/params['nstepsphi'])
    dphi = np.diff(phis)[0]
    thetas, phis = np.meshgrid(thetas, phis, indexing='ij')
    thetas = thetas.flatten()[np.newaxis,np.newaxis,:]
    phis = phis.flatten()[np.newaxis,np.newaxis,:]

    rho_ = rho[:,:,np.newaxis]
    z_ = z[:,:,np.newaxis]

    # rotate the coordinates, rescale their value such that after 
    # rounding it corresponds to some index in intensity_prof
    rho_r, z_r = rotate_cyl_coords_2angles_return_rho_z(
        rho=rho_, phi=0, z=z_, 
        alpha=phis, 
        beta=thetas
    )
    # get intensity from interpolation of pencil beam:
    ang_conv = np.sum(
        func_rho_z(np.abs(rho_r), z_r)\
        * np.sin(thetas) * dtheta * dphi * (rho_*rho_ + z_*z_),
        axis=2
    )
    
    norm = 2 * np.pi * (rho*rho + z*z)
    return ang_conv/norm

def calc_I_fiber(params):

    params = calc_dependent_params(params)
    # Calculate pencil beam of scattered
    ## Define space
    rho_max_pen, z_max_pen = calc_pencil_rho_z_max(
        params['theta_div'],
        params['xymax'],
        params['zmax']
    )
    
    if params['rho_exp_smpl']:
        # use exp.-sampling with 0 at beginning for rho
        rhos = log_smplng(
            min_=params['rhoexpmin'], 
            max_=rho_max_pen+params['rhoexpmin'], 
            n_smpls=params['n_rhosmpls']
        )
        rhos -= params['rhoexpmin']
    else:
        # use uniform sampling
        rhos = np.arange(
            0, 
            rho_max_pen+params['rhostep'],
            params['rhostep']
        )
    zs = np.arange(1, z_max_pen+params['dz'], params['dz'])

    ## Define multi-path time
    if params['tau_exp_smpl']:
        taus = log_smplng(
            min_=params['taumin'],
            max_=params['taumax'],
            n_smpls=params['n_tausmpls']
        )
    else:
        taus = np.arange(
            params['taumin'],
            params['taumax']+params['taustep'],
            params['taustep']
        )

    rho3_pen, z3_pen, tau3_pen = np.meshgrid(rhos, zs, taus, indexing='ij')
    pencil_beam_scattered = pencil_scattered_time_integrated(
        z3_pen, rho3_pen, tau3_pen, 
        g=params['g'], 
        mu_s=params['mu_s'], 
        mu_a=params['mu_a'], 
        c=params['c'], 
        G_version=params['mu_tau']
    )   

    # create interpolator for pencil_beam:
    interpolator_pencil_beam = Interpolator(
        rho3_pen[:,:,0], z3_pen[:,:,0], pencil_beam_scattered, fill_value=0
    )

    # Calculate scattered light
    ## define space
    rhos = np.arange(
        0, params['xymax']+params['dxy'],params['dxy']
    )
    zs = np.arange(
        1, params['zmax']+params['dz'], params['dz'])

    rho2_cone, z2_cone = np.meshgrid(rhos, zs, indexing='ij')

    # cone
    cone_scattered = ang_conv(
        rho2_cone, z2_cone, interpolator_pencil_beam.calc, 
        params = params
    )
    # create Interpolator for cone_scattered
    interp_cone_scattered = Interpolator(rho2_cone, z2_cone, cone_scattered)
    # disk
    disk_scattered = disk_conv_numpy(
        rho=rho2_cone, 
        z=z2_cone, 
        I_rho_z=interp_cone_scattered.calc, 
        opt_radius=params['opt_radius'], 
        dxy=params['dxy_scattered_disk']
    )
    # Calculate direct light
    def I_direct_cone_fixed_params(rho, z):
        return I_direct_cone(z, rho, params)
    disk_direct = disk_conv_numpy(
        rho=rho2_cone, 
        z=z2_cone, 
        I_rho_z=I_direct_cone_fixed_params, 
        opt_radius=params['opt_radius'], 
        dxy=params['dxy_direct_disk']
    )
    
    # results
    results = dict(
        pencil = dict(rho=rho3_pen[:,:,0], z=z3_pen[:,:,0], scattered=pencil_beam_scattered),
        cone = dict(rho=rho2_cone, z=z2_cone, scattered=cone_scattered),
        final = dict(rho=rho2_cone, z=z2_cone, scattered=disk_scattered, direct=disk_direct, combined=disk_direct+disk_scattered)
    )
    return results

def calc_I_fiber_reproduce_error(
        params, 
        fix_pencil_rho_z, 
        ang_conv_fill_nan,
        decrease_disk_conv_volume_radially
        ):

    params = calc_dependent_params(params)
    # Calculate pencil beam of scattered
    ## Define space
    rho_max_pen, z_max_pen = calc_pencil_rho_z_max(
        params['theta_div'],
        params['xymax'],
        params['zmax']
    )
    if fix_pencil_rho_z:
        rho_max_pen = fix_pencil_rho_z[0]
        z_max_pen = fix_pencil_rho_z[1]
    
    if params['rho_exp_smpl']:
        # use exp.-sampling with 0 at beginning for rho
        rhos = np.insert(log_smplng(
            min_=params['rhoexpmin'], 
            max_=rho_max_pen, 
            n_smpls=params['n_rhosmpls']-1
        ), 0, 0)
    else:
        # use uniform sampling
        rhos = np.arange(
            0, 
            rho_max_pen+params['rhostep'],
            params['rhostep']
        )
    zs = np.arange(1, z_max_pen+params['dz'], params['dz'])

    ## Define multi-path time
    if params['tau_exp_smpl']:
        taus = log_smplng(
            min_=params['taumin'],
            max_=params['taumax'],
            n_smpls=params['n_tausmpls']
        )
    else:
        taus = np.arange(
            params['taumin'],
            params['taumax']+params['taustep'],
            params['taustep']
        )

    rho3_pen, z3_pen, tau3_pen = np.meshgrid(rhos, zs, taus, indexing='ij')
    pencil_beam_scattered = pencil_scattered_time_integrated(
        z3_pen, rho3_pen, tau3_pen, 
        g=params['g'], 
        mu_s=params['mu_s'], 
        mu_a=params['mu_a'], 
        c=params['c'], 
        G_version=params['mu_tau']
    )   

    # create interpolator for pencil_beam:
    if ang_conv_fill_nan == True:
        interpolator_pencil_beam = Interpolator(
            rho3_pen[:,:,0], z3_pen[:,:,0], pencil_beam_scattered, fill_value=np.nan
        )
    else:
        interpolator_pencil_beam = Interpolator(
            rho3_pen[:,:,0], z3_pen[:,:,0], pencil_beam_scattered, fill_value=0
        )

    # Calculate scattered light
    ## define space
    rhos = np.arange(
        0, params['xymax']+params['dxy'],params['dxy']
    )
    zs = np.arange(
        1, params['zmax']+params['dz'], params['dz'])

    rho2_cone, z2_cone = np.meshgrid(rhos, zs, indexing='ij')

    # cone
    cone_scattered = ang_conv(
        rho2_cone, z2_cone, interpolator_pencil_beam.calc, 
        params = params
    )
    if ang_conv_fill_nan == True:
        cone_scattered[np.isnan(cone_scattered)] = 0
    # create Interpolator for cone_scattered
    if decrease_disk_conv_volume_radially:
        rho_idx_rm = int(decrease_disk_conv_volume_radially / params['dxy'])
        interp_cone_scattered = Interpolator(
            rho2_cone[:-rho_idx_rm,:], z2_cone[:-rho_idx_rm,:], cone_scattered[:-rho_idx_rm,:], fill_value=0
        )
    else:
        interp_cone_scattered = Interpolator(rho2_cone, z2_cone, cone_scattered)
    # disk
    disk_scattered = disk_conv_numpy(
        rho=rho2_cone, 
        z=z2_cone, 
        I_rho_z=interp_cone_scattered.calc, 
        opt_radius=params['opt_radius'], 
        dxy=params['dxy_scattered_disk']
    )
    if decrease_disk_conv_volume_radially == True:
        disk_scattered[np.isnan(disk_scattered)] = 0
    # Calculate direct light
    def I_direct_cone_fixed_params(rho, z):
        return I_direct_cone(z, rho, params)
    disk_direct = disk_conv_numpy(
        rho=rho2_cone, 
        z=z2_cone, 
        I_rho_z=I_direct_cone_fixed_params, 
        opt_radius=params['opt_radius'], 
        dxy=params['dxy_direct_disk']
    )
    
    # results
    results = dict(
        pencil = dict(rho=rho3_pen[:,:,0], z=z3_pen[:,:,0], scattered=pencil_beam_scattered),
        cone = dict(rho=rho2_cone, z=z2_cone, scattered=cone_scattered),
        final = dict(rho=rho2_cone, z=z2_cone, scattered=disk_scattered, direct=disk_direct, combined=disk_direct+disk_scattered)
    )
    return results
