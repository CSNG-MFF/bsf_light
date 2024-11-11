import numpy as np
from scipy.interpolate import griddata

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
