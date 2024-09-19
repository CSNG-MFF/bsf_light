from BSF import pencil_scattered, G, h, I_direct_cone
from utils import calc_dependent_params, log_smplng, reformat_error_data
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd

fs=8
font = {'size'   : fs}
matplotlib.rc('font', **font)
plt.rcParams.update({
    "font.family": "Arial"
})
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
    "tau_exp_smpl": False,
    "taumin": 5,  # fs 
    "taumax": 10000, # fs
    "n_tausmpls": 100, # if exp-smpl True, see definition of exp-sampling for rho above
    "taustep": 50, # fs stepsize if exp-smpl False
    
    "mu_tau": "eq4", # which equations to use to calculate the first moment 
                                   # of time dispersion, mu. Should be one of
                                   # 'eq4'               > McLean eq. 4
                                   # 'table1_vandeHulst' 
                                   #             -> McLean table1 van de Hulst & Kattawar
                                   # 'table1_Dolin' 
                                   #  -> McLean table1 Dolin, Ishimaru, Lutomirski et al.
                                   # 'table1_Stotts' 
                                   #             -> McLean table1 Stotts

    
    # angular convolution
    "nstepstheta": 24,  # ang conv steps
    "nstepsphi"  : 24,
}
params = calc_dependent_params(params)
zs = np.arange(1, 700, 5)
t_exp = log_smplng(
    min_=params['taumin'],
    max_=params['taumax'],
    n_smpls=params['n_tausmpls']
)
t_u5 = np.arange(
    params['taumin'],
    params['taumax']+5,
    5
)
t_u50 = np.arange(
    params['taumin'],
    params['taumax']+50,
    50
)
t_u500 = np.arange(
    params['taumin'],
    params['taumax']+500,
    500
)
zz_exp, tt_exp = np.meshgrid(zs, t_exp, indexing='ij')
zz_u5, tt_u5 = np.meshgrid(zs, t_u5, indexing='ij')
zz_u50, tt_u50 = np.meshgrid(zs, t_u50, indexing='ij')
zz_u500, tt_u500 = np.meshgrid(zs, t_u500, indexing='ij')

rho = 0
zz = [zz_exp, zz_u5, zz_u50, zz_u500]
tt = [tt_exp, tt_u5, tt_u50, tt_u500]
labels = ['log', 'uni@5fs', 'uni@50fs', 'uni@500fs']

pencils = [pencil_scattered(
    z, rho, t, g=params['g'], mu_s=params['mu_s'], mu_a=params['mu_a'], c=params['c'], G_version=params['mu_tau']
) for z, t in zip(zz,tt)
         ]
pencils_integrated = [
    np.sum(pencil[:,:-1] * np.diff(tt, axis=1), axis=1) for pencil,tt in zip(pencils, tt)
]
pencil_direct = I_direct_cone(zs, rho, params)


fs=8
A4_w, A4_h = 8.27, 11.69 # inch
fig, axs = plt.subplots(ncols=2, nrows=3, figsize=(A4_w*0.8, A4_h*0.6))
axs = axs.flatten()
for ax in axs:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

# pencil beam scattered before tau-integration for z=101 um
z_idx = 20
for tt_, zz_, pen_, label in zip(tt, zz, pencils, labels):
    axs[0].plot(tt_[z_idx,:], pen_[z_idx,:], label=label, marker='.')
axs[0].text(x=0.2, y=0.2, s="",transform=ax.transAxes, fontsize=fs)

for pi, label in zip(pencils_integrated, labels):
    axs[1].plot(zs, pi, label=label)
    axs[2].plot(zs, pencil_direct+pi, label=label)

# final beam data
final_exp = np.load('results/c0~0.299792458/ntissue~1.36/mu_a~6e-05/mu_s~0.0211/g~0.86/NA~0.37/opt_radius~100/xymax~700/dxy~5/zmax~700/dz~5/rho_exp_smpl~True/rhoexpmin~1.0/n_rhosmpls~20/rhostep~2/tau_exp_smpl~True/taumin~5.0/taumax~10000/n_tausmpls~100/taustep~5/mu_tau~eq4/nstepstheta~24/nstepsphi~24.pickle', allow_pickle=True)
final_u5 = np.load('results/c0~0.299792458/ntissue~1.36/mu_a~6e-05/mu_s~0.0211/g~0.86/NA~0.37/opt_radius~100/xymax~700/dxy~5/zmax~700/dz~5/rho_exp_smpl~True/rhoexpmin~1.0/n_rhosmpls~20/rhostep~2/tau_exp_smpl~False/taumin~5.0/taumax~10000/n_tausmpls~100/taustep~5/mu_tau~eq4/nstepstheta~24/nstepsphi~24.pickle', allow_pickle=True)
final_u50 = np.load('results/c0~0.299792458/ntissue~1.36/mu_a~6e-05/mu_s~0.0211/g~0.86/NA~0.37/opt_radius~100/xymax~700/dxy~5/zmax~700/dz~5/rho_exp_smpl~True/rhoexpmin~1.0/n_rhosmpls~20/rhostep~2/tau_exp_smpl~False/taumin~5.0/taumax~10000/n_tausmpls~100/taustep~50/mu_tau~eq4/nstepstheta~24/nstepsphi~24.pickle', allow_pickle=True)
final_u500 = np.load('results/c0~0.299792458/ntissue~1.36/mu_a~6e-05/mu_s~0.0211/g~0.86/NA~0.37/opt_radius~100/xymax~700/dxy~5/zmax~700/dz~5/rho_exp_smpl~True/rhoexpmin~1.0/n_rhosmpls~20/rhostep~2/tau_exp_smpl~False/taumin~5.0/taumax~10000/n_tausmpls~100/taustep~500/mu_tau~eq4/nstepstheta~24/nstepsphi~24.pickle', allow_pickle=True)
final = [final_exp, final_u5,final_u50, final_u500] 
# normalize by maximum of uniform-sampled beam with roughest sampling
norm = np.max(final_u500['final']['combined'])
for final_, label in zip(final, labels):
    z_f = final_['final']['z']
    x_f = final_['final']['x']
    I_f = final_['final']['combined'] / norm
    x_0, y_0 = (np.array(z_f.shape)[:2]/2).astype(int)
    axs[3].plot(z_f[x_0, y_0, :], I_f[x_0, y_0, :], label=label)
    axs[4].plot(x_f[:, y_0, int(300/params['dz'])], I_f[:, y_0, int(300/params['dz'])]/I_f[x_0, y_0, int(300/params['dz'])], 
                label=label)
    axs[5].plot(x_f[:, y_0, int(600/params['dz'])], I_f[:, y_0, int(600/params['dz'])]/I_f[x_0, y_0, int(600/params['dz'])], 
                label=label)

## import scraped data from Fig 5a (decay over depth)
data = pd.read_csv('data_from_paper/2205_David_a_points.txt', names=['x','y'], delimiter=';')
data_err = pd.read_csv('data_from_paper/2205_David_a_err.txt', names=['x','y'], delimiter=';')
# sort errors and substract datapoint to adapt for plt.errorbar arguments
yerrs = reformat_error_data(data, data_err)
axs[3].errorbar(data.x*1000,
            data.y/100,
            yerr=yerrs/100,
            label='experiment', 
            marker='x',
            c='black', 
            linestyle='',
            capsize=5)

## import scraped data from Fig 5b (radial profiles at z=300um, 600 um)
data_600 = pd.read_csv('data_from_paper/2205_David_b_points_red.txt', names=['x','y'], delimiter=';')
data_300 = pd.read_csv('data_from_paper/2205_David_b_points_blue.txt', names=['x','y'], delimiter=';')
axs[4].plot(data_300.x*1000,data_300.y/100, label='experiment', 
            marker='x',c='black', linestyle='')
axs[5].plot(data_600.x*1000,data_600.y/100, label='experiment', 
            marker='x',c='black', linestyle='')


# set limits
axs[0].set_xlim(0,500)
for ax in axs[1:]:
    ax.set_xlim(0,700)
for ax in axs:
    ax.set_ylim(0,None)
axs[2].set_ylim(0,0.5)
axs[3].set_ylim(0,0.5)
    
# set ax ticks, labels
for ax in axs[:4]:
    ax.set_ylabel('Transmission', fontsize=fs)
    ax.set_xlabel('Depth z [µm]', fontsize=fs)
for ax in axs[3:]:
    ax.set_ylabel('Normalized transmission')
    ax.set_xlabel('Depth z [µm]', fontsize=fs)
axs[0].set_yticks([0,0.5e-5,1e-5,1.5e-5])
axs[0].set_xlabel('Multipath time [fs]', fontsize=fs)
    
axs[0].set_title("Pencil beam (scattered) @ z=101µm\n before time integration", fontsize=fs, y=0.95)
axs[1].set_title("Pencil beam (scattered) after\nintegration over multipath time", fontsize=fs, y=0.95)
axs[2].set_title('Pencil (direct+scattered)', fontsize=fs, y=0.95)
axs[3].set_title('Final beam spread along depth', fontsize=fs, y=0.95)
axs[4].set_title('Final beam spread along radial direction\ndepth = 300 µm', fontsize=fs, y=0.92)
axs[5].set_title('Final beam spread along radial direction\ndepth = 600 µm', fontsize=fs, y=0.92)

axs[-1].legend(fontsize=fs, frameon=False)

plt.subplots_adjust(hspace=0.4)
for ax, letter in zip(axs, ['a','b','c','d','e','f']):
    ax.text(-0.17, 1.05, letter, fontsize=10, fontweight='bold', color='black', transform=ax.transAxes)
fig.savefig('figures/figure3.png', dpi=300)
plt.close()
