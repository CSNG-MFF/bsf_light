import numpy as np
from bsf_light.utils import log_smplng, calc_pencil_rho_z_max, Interpolator, calc_dependent_params
from bsf_light.convolutions import ang_conv, disk_conv
from bsf_light.direct import I_direct_cone
from bsf_light.scattered import pencil_scattered_time_integrated

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
    disk_scattered = disk_conv(
        rho=rho2_cone, 
        z=z2_cone, 
        func_rho_z=interp_cone_scattered.calc, 
        opt_radius=params['opt_radius'], 
        dxy=params['dxy_scattered_disk']
    )
    # Calculate direct light
    def I_direct_cone_fixed_params(rho, z):
        return I_direct_cone(z, rho, params)
    disk_direct = disk_conv(
        rho=rho2_cone, 
        z=z2_cone, 
        func_rho_z=I_direct_cone_fixed_params, 
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
    disk_scattered = disk_conv(
        rho=rho2_cone, 
        z=z2_cone, 
        func_rho_z=interp_cone_scattered.calc, 
        opt_radius=params['opt_radius'], 
        dxy=params['dxy_scattered_disk']
    )
    if decrease_disk_conv_volume_radially == True:
        disk_scattered[np.isnan(disk_scattered)] = 0
    # Calculate direct light
    def I_direct_cone_fixed_params(rho, z):
        return I_direct_cone(z, rho, params)
    disk_direct = disk_conv(
        rho=rho2_cone, 
        z=z2_cone, 
        func_rho_z=I_direct_cone_fixed_params, 
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
