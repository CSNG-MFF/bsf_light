from scipy.io import loadmat
import numpy as np
import pandas as pd

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

def load_published_exp_data():
    """
    Returns 2 dictionaries with depth and radial data.

    Spatial units all [um], light transmission decimals.
    """
    ## import scraped exp data from Fig 5a (decay over depth)
    depth_data = pd.read_csv('data_from_paper/2205_David_a_points.txt', names=['x','y'], delimiter=';')
    depth_data_err = pd.read_csv('data_from_paper/2205_David_a_err.txt', names=['x','y'], delimiter=';')
    depth = dict()
    # sort errors and substract datapoint to adapt for plt.errorbar arguments
    depth['z'] = depth_data.x.values*1000 # convert from mm to um
    depth['transmission'] = depth_data.y.values/100 # convert from percent to decimal
    depth['transmission_err'] = reformat_error_data(depth_data, depth_data_err)/100 # convert from percent to decimal
    
    ## import scraped exp data from Fig 5b (radial profiles at z=300um, 600 um)
    radial_data_z600 = pd.read_csv('data_from_paper/2205_David_b_points_red.txt', names=['x','y'], delimiter=';')
    radial_data_z300 = pd.read_csv('data_from_paper/2205_David_b_points_blue.txt', names=['x','y'], delimiter=';')
    radial = dict()
    radial['x_z300'] = radial_data_z300.x.values*1000 # convert from mm to um
    radial['transmission_z300'] = radial_data_z300.y.values/100 # convert from percent to decimal
    radial['x_z600'] = radial_data_z600.x.values*1000 # convert from mm to um
    radial['transmission_z600'] = radial_data_z600.y.values/100 # convert from percent to decimal

    return depth, radial

def load_published_model_data():
    """
    Returns 2 dictionaries with depth and radial data.

    Spatial units all [um], light transmission decimals.
    """
    depth = dict()
    
    ## import model prediction from Fig 5b
    radial_curve_z300 = np.load('data_from_paper/BSF_model_radial_curve_300.npy')
    radial_curve_z600 = np.load('data_from_paper/BSF_model_radial_curve_600.npy')
    radial = dict()
    radial['x_z300'] = radial_curve_z300[:,0]*1000 # convert from mm to um
    radial['transmission_z300'] = radial_curve_z300[:,1]/100 # convert from percent to decimal
    radial['x_z600'] = radial_curve_z600[:,0]*1000 # convert from mm to um
    radial['transmission_z600'] = radial_curve_z600[:,1]/100 # convert from percent to decimal

    return depth, radial

def load_matlab_model_data():
    matlab_data = loadmat('matlab/output3D_26_08_2024.mat')['out']
    matlab_data_xz = matlab_data[int(matlab_data.shape[0]/2):,int(matlab_data.shape[1]/2),:]
    matlab_dx = 5
    data = dict()
    data['data'] = matlab_data_xz
    matlab_x = np.arange(0,81*matlab_dx,matlab_dx)
    matlab_z = np.arange(0,140*5,5)
    data['xx'], data['zz'] = np.meshgrid(matlab_x, matlab_z, indexing='ij')
    data['z300'] = int(300/matlab_dx)
    data['z600'] = int(600/matlab_dx)
    return data