import numpy as np
from bsf_light.utils import rotate_cyl_coords_2angles_return_rho_z

def ang_conv(rho, z, func_rho_z, params):
    """
    Numerical angular convolution. 

    Used to obtain angular convolution of scattered light to
    calculate the scattered component of a light cone emitted
    from an infinitesimal point in emitter surface.

    Performs convolution over angles theta and phi of function
    depending on rho and z.

    Parameters
    ----------
    z : array_like
        Distance along the z-axis, representing the propagation depth of the light.
    rho : array_like
        Distance from the center in the xy-plane.
    func_rho_z: function
        Function depending on rho and z (1st, 2nd arguments). Angular convolution
        will be performed over this function.
    params : dict
        Parameters used in the calculation, refer to the `calc_I_fiber` function in the
        `fiber.py` module.

    Returns
    -------
    I_res: array_like
        Intensities resulting from angular convolution in the shape of rho, z.

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
    
    I_res = ang_conv / (2 * np.pi * (rho*rho + z*z))
    return I_res

def disk_conv(rho, z, func_rho_z, opt_radius: float, dxy: float):
    """
    Numerical disk convolution. 

    Used to generalize from a light cone to the light emitted from a
    circular surface.

    Warning: Makes use of symmetry along y-axis. Instead of calculating
    contribution from all 4 x-y-quadrants, calculates only quadrants
    with positive y and multiplies by 2.

    Parameters
    ----------
    rho : array_like
        Distance from the center in the xy-plane.
    z : array_like
        Distance along the z-axis, representing the propagation depth of the light.
    func_rho_z: function
        Function depending on rho and z (1st, 2nd arguments). Disk convolution
        will be performed over this function.
    opt_radius:
        Radius of the optical fiber as given in params. For details refer to the params
        dictionary explained in docstring of `calc_I_fiber` function in the `fiber.py`
        module.
    dxy:
        Step size of the convolution in x/y-direction. For details refer to the params 
        dictionary explained in docstring of `calc_I_fiber` function in the `fiber.py`
        module.

    Returns
    -------
    I_res: array_like
        Intensities resulting from disk convolution in the shape of rho, z.
    """

    x_shift = np.arange(-1 * opt_radius, opt_radius + dxy, dxy)
    y_shift = np.arange(0, opt_radius + dxy, dxy)
    xx_shift, yy_shift = np.meshgrid(x_shift, y_shift, indexing='ij')
    
    # Apply disk constraint to the shifts
    within_disk = (xx_shift ** 2 + yy_shift ** 2) <= opt_radius ** 2
    xx_shift = xx_shift[within_disk].flatten()
    yy_shift = yy_shift[within_disk].flatten()

    # Calculate shifted coordinates
    rho_shifted = np.sqrt((rho[:, :, np.newaxis] - xx_shift) ** 2 + yy_shift ** 2)
    
    z_shifted = z[:, :, np.newaxis] + np.zeros(xx_shift.shape)
    # Interpolate over the shifted coordinates in a vectorized manner
    I_res = np.sum(func_rho_z(rho_shifted, z_shifted) * 2 * dxy**2, axis=2)

    return I_res
