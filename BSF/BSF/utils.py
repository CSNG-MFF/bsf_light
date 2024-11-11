import numpy as np
from scipy.interpolate import griddata

def calc_pencil_rho_z_max(theta, xmax, zmax):
    """
    Calculate dimensions of pencil beam volume to be sufficient during convolutions.

    Parameters
    ----------
    theta : float
        Opening angle of optical fiber in rad.
    xmax : float
        Max. length in x/y direction.
    zmax : float
        Max. length in depth direction.

    Returns
    -------
    rho, z: float, float
        Dimensions of pencil beam.
    """
    if np.abs(theta) > np.pi/4:
        theta = np.pi/4
    z_pencil = np.sqrt(xmax**2 + zmax**2)
    rho_pencil = (xmax * np.cos(theta)+ zmax * np.sin(theta)) 
    return rho_pencil, z_pencil

def log_smplng(min_, max_, n_smpls):
    """
    Generate log-distributed samples using natural logarithm.

    Parameters
    ----------
    min_ : float
        Start value of interval.
    max_ : float
        End of interval.
    n_smpls : int
        Number of samples

    Returns
    -------
    samples : array_like
        Log-distributed samples.
    """
    return np.exp(np.linspace(np.log(min_), np.log(max_), n_smpls))

def calc_dependent_params(params):
    """
    Calculate dependent parameters of param dictionary.

    For details on params, see function ´calc_I_fiber´ in fiber.py module.
    Calculates and adds speed of light in medium and opening angle of fiber
    to params.

    Parameters
    ----------
    params: dict
        Params of simulation, see function ´calc_I_fiber´ in fiber.py for details.
    Returns
    -------
    params: dict
        Params of simulation, see function ´calc_I_fiber´ in fiber.py for details.
    """
    # speed of light in medium
    params['c'] = params['c0'] / params['ntissue']
    # divergence of light emitted from optical fiber
    params['theta_div'] = np.arcsin(params['NA'] / params['ntissue'])
    return params

def rotate_cyl_coords_2angles_return_rho_z(rho, phi, z, alpha, beta):
    """
    Rotation in cylinder coordinates.

    Takes cylindrical coordinates (rho, phi, z) of the coord-
    system in which the Riemann sum representing the angular
    convolution takes place and returns cylindrical coords
    (rho_, z_) in which the light beam is oriented
    along the z_-axis.

    Parameters
    ----------
    rho : array_like
        Distance from the center in the xy-plane.
    phi : array_like
        Azimuthal angle to point in cylinder coorinates.
    z : array_like
        Distance along the z-axis, representing the propagation depth of the light.
    alpha : float
        Angle of rotation around z-axis in rad.
    beta : float
        Angle of rotation around y-axis in rad.

    Returns
    -------
    rho_, z_: array_like
        Rotated coordinates in shapes as rho, z

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

    return rho_, z_

class Interpolator:
    """
    Class used for interpolation during beam calculation.

    Ensures that arrays of arbitrary shaped can be interpolated and
    are returned in according shapes.

    Can be initialized with data with 2D dependency on rr, zz:
    interp = Interpolator(rr, zz, data) using optional fill_value.

    Interpolation of points using interp.calc(r,z).

    Parameters
    ----------
    rr : array_like
        Radial coordinate of arbitrary shape.
    zz : array_like
        Depth coordinate of arbitrary shape.
    data: array_like
        Corresponding reference data of shape as rr, zz.
    fill_value: float
        Which value to fill (float) or not fill (np.nan).

    Returns
    -------
    Interpolator : obj
        Interpolation object, use Interpolator.calc(r,z) to interpolate.
    """
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
    
    Parameters
    ----------
    arr : np.ndarray
        The input 2D array of shape (a, b), where a is the number of rows.
    make_neg : bool
        Whether to multiply the mirrored part with (-1).
    
    Returns
    -------
    result : np.ndarray
        The mirrored array of shape (2a-1, b).
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
