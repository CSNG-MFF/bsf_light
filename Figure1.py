from BSF import calc_I_fiber, I_direct_cone
from utils import reformat_error_data
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import gridspec
from scipy.io import loadmat
import numpy as np
import pandas as pd
from time import time

params = {
    # natural constants & optical params
    "c0": 0.299792458, #(* um/fs *),
    # tissue properties
    "ntissue": 1.36, # refractive index of cortical tissue
    "mu_a": 0.00006, # um**-1
    "mu_s": 0.0211, # um**-1
    "g": 0.86, # anisotropy
    # optical fiber
    "NA": 0.37, # num aperture
    "opt_radius": 100, # fiber radius in um
    
    
    # final volume, NOTE: choose xymax to be 100um larger than your desired 
    #                     volume to avoid artifacts from disk convolution.
    #                     Check docstring of disk_conv_nproll in utils.py for
    #                     more information.
    'xymax':  700,
    'dxy':   5,
    'zmax':   700,
    'dz':    5,

    
    # Calculation of scattered pencil beam
    # exp.-sampling of rho
    "rho_exp_smpl" : True,
    "rhoexpmin": 1, # exp(linspace(log(rhomin), log(rhomax), n_rhosmpls)) if exp-smpl True
    "n_rhosmpls" : 20, # (line above)
    "rhostep": 2, # stepsize if exp-smpl False
    
    # multipath-time integral
    "tau_exp_smpl": True,
    "taumin": 5,  # fs 
    "taumax": 10000, # fs
    "n_tausmpls": 100, # if exp-smpl True, see definition of exp-sampling for rho above
    "taustep": 500, # fs stepsize if exp-smpl False
    
    "mu_tau": "eq4", # which equations to use to calculate the first moment 
                                   # of time dispersion, mu. Should be one of
                                   # 'eq4'               > McLean eq. 4
                                   # 'table1_vandeHulst' 
                                   #             -> McLean table1 van de Hulst & Kattawar
                                   # 'table1_Lutomirski' 
                                   #  -> McLean table1 Dolin, Ishimaru, Lutomirski et al.
    
    # angular convolution
    "nstepstheta": 24,  # ang conv steps
    "nstepsphi"  : 24,
}

#################################################### Calculation

st = time()
results = calc_I_fiber(params)
comp_time = time() - st
print(f"Computing took {comp_time} seconds.")

#################################################### Plots

matlab_data = loadmat('matlab/output3D_26_08_2024.mat')['out']
matlab_data_xz = matlab_data[:,int(matlab_data.shape[1]/2),:]
matlab_dx = 5
matlab_x = np.arange(-80*matlab_dx,81*matlab_dx,matlab_dx)
matlab_z = np.arange(0,140*5,5)
matlab_xx, matlab_zz = np.meshgrid(matlab_x, matlab_z, indexing='ij')

this_model = results['final']['combined']
this_model = this_model/this_model.max()
this_model[this_model==0] = 1e-30
shape = this_model.shape
y_cnt = int(shape[1]/2)
z_300 = int(300/params['dz'])
z_600 = int(600/params['dz'])
x = results['final']['x']
z = results['final']['z']

## import scraped data from Fig 5a (decay over depth)
data = pd.read_csv('data_from_paper/2205_David_a_points.txt', names=['x','y'], delimiter=';')
data_err = pd.read_csv('data_from_paper/2205_David_a_err.txt', names=['x','y'], delimiter=';')
# sort errors and substract datapoint to adapt for plt.errorbar arguments
yerrs = reformat_error_data(data, data_err)
## import scraped data from Fig 5b (radial profiles at z=300um, 600 um)
data_600 = pd.read_csv('data_from_paper/2205_David_b_points_red.txt', names=['x','y'], delimiter=';')
data_300 = pd.read_csv('data_from_paper/2205_David_b_points_blue.txt', names=['x','y'], delimiter=';')
curve_600 = np.load('data_from_paper/BSF_model_radial_curve_600.npy')
curve_300 = np.load('data_from_paper/BSF_model_radial_curve_300.npy')

# use data point at x=300um to normalize
norm = data.y[2]/100 / this_model[y_cnt, y_cnt, int(300/params['dz'])]

fs = 8
A4_w, A4_h = 8.27, 11.69 # inch
fig = plt.figure(figsize=(A4_w*0.8, A4_h*0.4))
axs = []
gs_row1 = gridspec.GridSpec(2, 4, width_ratios=[1,1,0.1,0.1])
gs_row2 = gridspec.GridSpec(2, 3)
axs.append(fig.add_subplot(gs_row1[0,0]))
axs.append(fig.add_subplot(gs_row1[0,1]))
axs.append(fig.add_subplot(gs_row1[0,2]))
axs.append(fig.add_subplot(gs_row2[1,0]))
axs.append(fig.add_subplot(gs_row2[1,1]))
axs.append(fig.add_subplot(gs_row2[1,2]))

# 2D profiles
logminmax = (1e-4,1)
cmap = 'Blues'
axs[0].set_title('Original (matlab app)', fontsize=fs)
mappable = axs[0].pcolormesh(matlab_xx, matlab_zz, matlab_data_xz, norm=LogNorm(*logminmax), shading='nearest', cmap=cmap)
axs[1].set_title('Replication', fontsize=fs)
mappable = axs[1].pcolormesh(x[:,y_cnt,:], z[:,y_cnt,:], this_model[:,y_cnt,:]*0.07, norm=LogNorm(*logminmax), shading='nearest', cmap=cmap)
for ax in axs[0:2]:
    ax.set_xlim(-400,400)
    ax.set_ylim(0,700)
cbar = plt.colorbar(mappable, cax=axs[2])
cbar.set_label('Normalized intensity', fontsize=fs)

for ax in axs[:2]:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])

    
# depth profile
# axs[3].set_title('Decay into depth', fontsize=fs)
axs[3].errorbar(data.x*1000,
            data.y/100,
            yerr=yerrs/100,
            label='experiment', 
            marker='x',
            c='black', 
            linestyle='',
            capsize=5)
axs[3].plot(z[y_cnt,y_cnt, :], this_model[y_cnt, y_cnt, :]*norm,
            ls='solid', color='orange', label='replication')
axs[3].plot(matlab_z, matlab_data_xz[int(matlab_data.shape[0]/2), :], 
            ls='solid', color='blue', label='original model')
#axs[3].plot(z[y_cnt,y_cnt, :], I_direct_cone(z=z[y_cnt,y_cnt, :], rho=0, params=params), color='gray', ls='dashed')
axs[3].set_xlim(0,700)
axs[3].set_ylim(0,0.3)
axs[3].set_ylabel('Normalized Intensity', fontsize=fs)
axs[3].set_xlabel('Depth [µm]', fontsize=fs)


# radial profiles
axs[4].set_title('depth = 300 µm', fontsize=fs)
axs[4].plot(data_300.x*1000,data_300.y/100, label='experiment', 
            marker='x',c='black', linestyle='')
axs[4].plot(x[:,y_cnt, 0], this_model[:, y_cnt, z_300]/this_model[y_cnt, y_cnt, z_300],
            ls='solid', color='orange', label='replication')
axs[4].plot(matlab_x, matlab_data_xz[:, int(300/5)]/matlab_data_xz[int(matlab_data.shape[1]/2), int(300/5)],
            ls='solid', color='blue', label='original model')
axs[4].plot(curve_300[:,0]*1000, curve_300[:,1]/100,
            ls='--', color='blue', label='original model')


axs[5].set_title('depth = 600 µm', fontsize=fs)
axs[5].plot(data_600.x*1000,data_600.y/100, label='experiment', 
            marker='x',c='black', linestyle='')
axs[5].plot(x[:,y_cnt, 0], this_model[:, y_cnt, z_600]/this_model[y_cnt, y_cnt, z_600],
            ls='solid', color='orange', label='replication')
axs[5].plot(matlab_x, matlab_data_xz[:, int(600/5)]/matlab_data_xz[int(matlab_data.shape[1]/2), int(600/5)],
            ls='solid', color='blue', label='original\n(matlab app)')
axs[5].plot(curve_600[:,0]*1000, curve_600[:,1]/100,
            ls='--', color='blue', label='original\n(publication)')

for ax in axs[4:]:
    ax.set_xlim(0,500)
    # ax.set_ylim(0,1)
    ax.set_xlabel('Radial distance [µm]', fontsize=fs)

for ax in axs[3:]:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

for ax in axs:
    ax.tick_params(axis='both', which='major', labelsize=fs)
axs[-1].legend(fontsize=fs, bbox_to_anchor=(0.6, 0.42), frameon=False)

# insert letters:
for ax, letter in zip(axs, ['a','b',None,'c','d','e']):
    if letter != None:
        ax.text(-0.1, 1.05, letter, fontsize=10, fontweight='bold', color='black', transform=ax.transAxes)
fig.savefig('figures/figure1.png', dpi=300)
plt.close()
