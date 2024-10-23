import numpy as np
from scipy.interpolate import griddata
from tqdm import tqdm

def conv_bool(inp):
    if (inp == 'True') or (inp == True):
        return True
    elif (inp == 'False') or (inp == False):
        return False
    else:
        raise ValueError(f"Value must be 'True' or 'False' but is {inp}")

def type_cast_paramdict(params):
    """
    Enforce data types for values in params.
    """
    params['c0']= float(params['c0'])
    params['ntissue']= float(params['ntissue'])
    params['mu_a']= float(params['mu_a'])
    params['mu_s']= float(params['mu_s'])
    params['g']= float(params['g'])
    params['NA']= float(params['NA'])
    params['opt_radius']= int(float(params['opt_radius']))

    params['xymax'] = float(params['xymax'])
    params['dxy'] = float(params['dxy'])
    params['zmax'] = float(params['zmax'])
    params['dz'] = float(params['dz'])

    params['rho_exp_smpl'] = conv_bool(params['rho_exp_smpl'])
    params['rhoexpmin'] = float(params['rhoexpmin'])
    params['n_rhosmpls'] = int(float(params['n_rhosmpls']))
    params['rhostep'] = float(params['rhostep'])

    params['tau_exp_smpl'] = conv_bool(params['tau_exp_smpl'])
    params['taumin']= float(params['taumin'])
    params['taumax']= float(params['taumax'])
    params['n_tausmpls'] = int(float(params['n_tausmpls']))
    params['taustep']= float(params['taustep'])

    params['mu_tau'] = str(params['mu_tau'])

    params['nstepstheta']= int(float(params['nstepstheta']))
    params['nstepsphi']=int(float(params['nstepsphi']))

    return params

def disk_conv_nproll(intensity, opt_radius: float, dxy: float):
    """
    Use disk convolution to get full 3D light profile from a light cone 
    exiting infinitesimal point mapped to a circular emitter surface.
    """
    
    # set all np.nan to zero
    intensity[np.isnan(intensity)] = 0
    # generate sampling points for light in disk
    nxy = 2 * opt_radius / dxy + 1

    assert nxy%1==0, 'choose optical radius to be multiple of dxy'
    nxy = int(nxy)
    ix_disk_from_center = np.linspace(-opt_radius/dxy,opt_radius/dxy,nxy)
    iy_disk_from_center = np.linspace(-opt_radius/dxy,opt_radius/dxy,nxy)
    ixx_disk_from_center, iyy_disk_from_center = np.meshgrid(ix_disk_from_center, iy_disk_from_center, indexing='ij')
    within_disk = ((ixx_disk_from_center*dxy)**2 + (iyy_disk_from_center*dxy)**2) <= opt_radius**2
    ixx_disk_from_center_flat = ixx_disk_from_center[within_disk].flatten()
    iyy_disk_from_center_flat = iyy_disk_from_center[within_disk].flatten()
    assert np.sum(ixx_disk_from_center_flat % 1) == 0,\
    'opt radius is not integer multiple of dx or nx is not chosen to return integer steps'
    assert np.sum(iyy_disk_from_center_flat % 1) == 0,\
    'opt radius is not integer multiple of dy or ny is not chosen to return integer steps'
    ixx_disk_from_center_flat = ixx_disk_from_center_flat.astype(int)
    iyy_disk_from_center_flat = iyy_disk_from_center_flat.astype(int)

    # disk conv
    res = np.zeros(np.shape(intensity))
    for x_ishift, y_ishift in tqdm(list(zip(ixx_disk_from_center_flat, iyy_disk_from_center_flat))):
        res += np.roll(
            np.roll(
                intensity, shift=x_ishift, axis=0
            ), shift=y_ishift, axis=1
        ) * dxy * dxy
    return res

def disk_conv_numpy(rho, z, I_rho_z, opt_radius: float, dxy: float):
    """
    Use disk convolution to generalize from a light cone existing an 
    infinitesimal point to the light emitted from a circular surface.
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


def calc_pencil_rho_z_max(theta, xmax, zmax):
    if np.abs(theta) > np.pi/4:
        theta = np.pi/4
    z_pencil = np.sqrt(xmax**2 + zmax**2)
    rho_pencil = (xmax * np.cos(theta)+ zmax * np.sin(theta)) 
    return rho_pencil, z_pencil

def log_smplng(min_, max_, n_smpls):
    """
    Log-distributed samples.
    """
    return np.exp(np.linspace(np.log(min_), np.log(max_), n_smpls))

def calc_dependent_params(params):
    # speed of light in medium
    params['c'] = params['c0'] / params['ntissue']
    # divergence of light emitted from optical fiber
    params['theta_div'] = np.arcsin(params['NA'] / params['ntissue'])
    return params

def rotate_cyl_coords_2angles_return_rho_z(rho, phi, z, alpha, beta):
    """
    Takes cylindrical coordinates (rho, phi, z) of the coord-
    system in which the Riemann sum representing the angular
    convolution takes place and returns cylindrical coords
    (rho_, z_) in which the light beam is oriented
    along the z_-axis.
    """
    # Ensure input shapes are compatible
    rho = np.asarray(rho)
    z = np.asarray(z)

    # Convert cylindrical to Cartesian coordinates
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)

    # Rotate (x,y,z) by alpha around z-axis
    cosb = np.cos(beta)
    cosa = np.cos(alpha)
    sinb = np.sin(beta)
    sina = np.sin(alpha)

    # Apply rotation
    x_ = cosa * cosb * x + sina * cosb * y + sinb * z
    y_ = -1 * sina * x + cosa * y
    z_ = -1 * cosa * sinb * x - sina * sinb * y + cosb * z

    # Convert back to cylindrical coordinates
    rho_ = np.sqrt(x_**2 + y_**2)

    # Handle phi_ calculation properly
    #phi_ = np.arctan2(y_, x_)  # Use arctan2 for correct quadrant handling

    # Return both rho_ and z_ (and optionally phi_)
    #return rho_, phi_, z_
    return rho_, z_

class Interpolator:
    def __init__(self, rr, zz, data, fill_value=np.nan):
        self.type = 'griddata'
        self.rr = rr
        self.zz = zz
        self.data = data
        self.fill_value=fill_value
    def calc(self,rr,zz):
        flat_interp_data = griddata(
            points=np.array([self.rr.flatten(), self.zz.flatten()]).T,
            values=self.data.flatten(),
            xi=np.array([rr.flatten(), zz.flatten()]).T,
            fill_value=self.fill_value
        )
        return flat_interp_data.reshape(rr.shape) 
    
def round_x(data):
    d = data.copy()
    d.x/=0.05
    d.x = d.x.round(decimals=0)
    d.x *= 0.05
    return d.sort_values(by=['x','y'])

def reformat_yerr(y_vals, y_errs):
    return [y - y_vals[int(idx/2)] for idx, y in enumerate(y_errs)]

def reformat_error_data(data, data_err):
    """Use to sort errors from experimental depth data"""
    d_err = round_x(data_err)
    errs = reformat_yerr(y_vals=data.y, y_errs=d_err.y)
    return np.abs(np.array(errs).reshape(int(len(errs)/2),2).T)

def mirror_x_axis(arr, make_neg=False):
    """
    Mirrors a 2D array along the x-axis (the first dimension), ignoring the first row (x = 0).
    
    Parameters:
    arr (numpy.ndarray): The input 2D array of shape (a, b), where a is the number of rows.
    
    Returns:
    numpy.ndarray: The mirrored array of shape (2a-1, b).
    """
    # Get the rows excluding the first row (x = 0)
    arr_positive_x = arr[1:, :]  # Shape (a-1, b)
    
    # Mirror the array along the x-axis by flipping along the first axis
    arr_mirrored = np.flip(arr_positive_x, axis=0)  # Shape (a-1, b)
    if make_neg:
        # Concatenate the original array with the mirrored array
        result = np.vstack((-1*arr_mirrored, arr))  # Shape (2a-1, b)
    else:
        # Concatenate the original array with the mirrored array
        result = np.vstack((arr_mirrored, arr))  # Shape (2a-1, b)
    return result
