import numpy as np
from BSF.utils import rotate_cyl_coords_2angles_return_rho_z

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

def disk_conv(rho, z, I_rho_z, opt_radius: float, dxy: float):
    """
    Use disk convolution to generalize from a light cone existing an 
    infinitesimal point to the light emitted from a circular surface.

    Warning: Makes use of symmetry along y-axis. Instead of calculating
    contribution from all 4 x-y-quadrants, calculates only quadrants
    with positive y and multiplies by 2.
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
    I_res = np.sum(I_rho_z(rho_shifted, z_shifted) * 2 * dxy**2, axis=2)

    return I_res
