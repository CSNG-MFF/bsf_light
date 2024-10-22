from BSF import calc_I_fiber, I_direct_cone
from utils import reformat_error_data
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import gridspec
import scipy
import numpy as np
import pandas as pd
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
    
    
    # final volume
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
    "taustep": 5, # fs stepsize if exp-smpl False
    
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

    # disk convolution
    "dxy_direct_disk"   : 5,
    "dxy_scattered_disk": 5,
}
results = dict()
calc_types = ['table1_Lutomirski', 'table1_vandeHulst', 'eq4']
labels = ['approximation by\nLutomirski et al.', 'approximation by\nvan de Hulst & Kattawar', 'prec. calc. of µ and\nrelation based on approx.\nby Lutomirski et al.']
colors = ['tab:orange','tab:blue', 'tab:green']
for ct in calc_types:
    params['mu_tau'] = ct
    results[ct] = calc_I_fiber(params)
labels = ['approximation by\nLutomirski et al.', 'approximation by\nvan de Hulst & Kattawar', f'prec. calc. of $\mu$ + relating $\sigma$ to $\mu$\nbased on approx.by Lutomirski et al.']
## import scraped data from Fig 5a (decay over depth)
data = pd.read_csv('data_from_paper/2205_David_a_points.txt', names=['x','y'], delimiter=';')
data_err = pd.read_csv('data_from_paper/2205_David_a_err.txt', names=['x','y'], delimiter=';')
# sort errors and substract datapoint to adapt for plt.errorbar arguments
yerrs = reformat_error_data(data, data_err)
## import scraped data from Fig 5b (radial profiles at z=300um, 600 um)
data_600 = pd.read_csv('data_from_paper/2205_David_b_points_red.txt', names=['x','y'], delimiter=';')
data_300 = pd.read_csv('data_from_paper/2205_David_b_points_blue.txt', names=['x','y'], delimiter=';')
norm_z400 = data.iloc[3]['y']/100

fs = 8
A4_w, A4_h = 8.27, 11.69 # inch
fig = plt.figure(figsize=(A4_w*0.8, A4_h*0.2))
axs = []
gs = gridspec.GridSpec(1, 3)
axs.append(fig.add_subplot(gs[0,0]))
axs.append(fig.add_subplot(gs[0,1]))
axs.append(fig.add_subplot(gs[0,2]))

# depth profile
axs[0].errorbar(data.x*1000,
            data.y/100,
            yerr=yerrs/100,
            label='experiment', 
            marker='x',
            c='black', 
            linestyle='',
            capsize=5)

for ct, label, c in zip(calc_types, labels, colors):
    this_model = results[ct]['final']['combined']
    z_300 = int(300/params['dz'])
    z_600 = int(600/params['dz'])

    this_model[this_model==0] = 1e-30
    
    x = results[ct]['final']['rho']
    z = results[ct]['final']['z']
    axs[0].plot(z[0, :], this_model[0, :],
                ls='solid', label=label, color=c)
    axs[1].plot(x[:, 0], this_model[:, z_300]/this_model[0, z_300],
                ls='solid', label=label, color=c)
    axs[2].plot(x[:, 0], this_model[:, z_600]/this_model[0, z_600],
                ls='solid', label=label, color=c)
    
axs[0].set_xlim(0,700)
axs[0].set_ylim(0,0.3)
axs[0].set_ylabel('Normalized Intensity', fontsize=fs)
axs[0].set_xlabel('Depth [µm]', fontsize=fs)


# radial profiles
axs[1].set_title('depth = 300 µm', fontsize=fs)
axs[1].plot(data_300.x*1000,data_300.y/100, label='experiment', 
            marker='x',c='black', linestyle='')


axs[2].set_title('depth = 600 µm', fontsize=fs)
axs[2].plot(data_600.x*1000,data_600.y/100, label='experiment', 
            marker='x',c='black', linestyle='')

for ax in axs[1:]:
    ax.set_xlim(0,500)
    ax.set_ylim(0,1)
    ax.set_xlabel('Radial distance [µm]', fontsize=fs)

for ax in axs:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=fs)
axs[-1].legend(fontsize=fs, bbox_to_anchor=(1., 0.42), frameon=False)

# insert letters:
fig.text(0.073, 0.95, 'a', fontsize=10, fontweight='bold', color='black')
fig.text(0.357, 0.95, 'b', fontsize=10, fontweight='bold', color='black')
fig.text(0.63, 0.95, 'c', fontsize=10, fontweight='bold', color='black')
fig.savefig('figures/figure4.png', dpi=300, bbox_inches='tight')
plt.close()
