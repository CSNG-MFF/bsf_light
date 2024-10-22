import os, sys
from utils import type_cast_paramdict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter
from utils import reformat_error_data
import matplotlib
from matplotlib import gridspec

fs=8
font = {'size'   : fs}
matplotlib.rc('font', **font)
plt.rcParams.update({
    "font.family": "Arial"
})

def get_params_from_file(file):
    varis = [fn for fn in file.split('/') if '~' in fn]
    params = dict()
    key_val_pairs = [var.split('~') for var in varis]
    for key, val in key_val_pairs:
        params[key] = val
    params['nstepsphi'] = params['nstepsphi'][:-7]
    return type_cast_paramdict(params)
    
def get_results_from_dir(dire):
    files = []
    for dirpath, dirnames, filenames in os.walk('./results/'):
        for filename in filenames:
            # Construct the full path to the file
            files.append(os.path.join(dirpath, filename))
    paramdicts = [get_params_from_file(file) for file in files]
    return files, paramdicts

def plot_results(files, labelparams, legtitle=None, styles=['solid', 'dashed'], labels=None):
    if legtitle == None:
        legtitle = ' '.join(labelparams)
    if styles == None or len(styles)!=len(files):
        styles = ['solid'] * len(files)
    ## import scraped data from Fig 5a (decay over depth)
    data = pd.read_csv('data_from_paper/2205_David_a_points.txt', names=['x','y'], delimiter=';')
    data_err = pd.read_csv('data_from_paper/2205_David_a_err.txt', names=['x','y'], delimiter=';')
    # sort errors and substract datapoint to adapt for plt.errorbar arguments
    yerrs = reformat_error_data(data, data_err)
    
    A4_w, A4_h = 8.27, 11.69 # inch
    fig = plt.figure(figsize=(A4_w*0.8, A4_h*0.2))
    axs = []
    gs_row1 = gridspec.GridSpec(1, 3, width_ratios=[1,1,1])
    axs.append(fig.add_subplot(gs_row1[0,0]))
    axs.append(fig.add_subplot(gs_row1[0,1]))
    axs.append(fig.add_subplot(gs_row1[0,2]))
    
    for i, (file, ls) in enumerate(zip(files, styles)):
        params = get_params_from_file(file)
        if labels == None:
            label = ' '.join(['='.join([key,str(params[key])]) for key in labelparams])
            label = ' '.join([str(params[key]) for key in labelparams])
        else:
            label = labels[i]
        res = np.load(file, allow_pickle=True)
        x = res['final']['rho']
        z = res['final']['z']
        I = res['final']['combined']
        xcnt = 0
        axs[0].plot(z[xcnt, :], I[xcnt, :], label=label, ls=ls)
        axs[1].plot(x[:,int(300/params['dz'])], I[:, int(300/params['dz'])]/I[xcnt, int(300/params['dz'])], label=label, ls=ls)
        axs[2].plot(x[:,int(600/params['dz'])], I[:, int(600/params['dz'])]/I[xcnt, int(600/params['dz'])], label=label, ls=ls)
        
    axs[0].set_xlabel('Depth [µm]')
    axs[0].set_ylabel('Transmission (not normalized)')
    axs[0].set_ylim(0,None)
    axs[1].set_ylabel('Transmission (normalized)')
    for ax in axs[1:]:
        ax.set_xlabel('Radial distance [µm]')
        ax.set_ylim(0,1.05)
    for ax in axs:    
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
    axs[-1].legend(title=legtitle, bbox_to_anchor=(0.6, 0.42), frameon=False)
    axs[1].set_title('Depth = 300 µm', fontsize=fs, y=0.95)
    axs[2].set_title('Depth = 600 µm', fontsize=fs, y=0.95)
    return fig, axs
files, paramdicts = get_results_from_dir('results/')
params_df = pd.DataFrame(paramdicts)
index_order = ['rhoexpmin', 'n_rhosmpls', 'tau_exp_smpl', 'taumin', 'taumax', 'n_tausmpls', 'taustep', 'nstepstheta', 'nstepsphi']
params_df = params_df.set_index(index_order, append=True)
def get_idxs(rhoexpmin, n_rhosmpls, tau_exp_smpl, taumin, taumax, n_tausmpls, taustep, nanglesteps):
    return list(
        params_df.loc[:,rhoexpmin, n_rhosmpls, tau_exp_smpl, taumin, taumax, n_tausmpls, taustep, nanglesteps, nanglesteps].index
    )
idxs = get_idxs(
    rhoexpmin=1,
    n_rhosmpls=20,
    tau_exp_smpl=True,
    taumin=5,
    taumax=[100,1000,10000,100000],
    n_tausmpls=100,
    taustep=5,
    nanglesteps=24
)
idxs = [idx[0] for idx in idxs]
fig, ax = plot_results(
    files=itemgetter(*idxs)(files), 
    labelparams=['taumax'], 
    legtitle='tau_max [fs]',
    styles=['solid','solid','solid','dashed']
)
plt.tight_layout()
for ax, letter in zip(ax, ['d','e','f']):
    ax.text(-0.17, 1.05, letter, fontsize=10, fontweight='bold', color='black', transform=ax.transAxes)
fig.suptitle("Variation of upper integration limit", fontsize=fs, y=0.99)
fig.savefig('figures/figure2def.png', dpi=300)
plt.show()

idxs = get_idxs(
    rhoexpmin=1,
    n_rhosmpls=20,
    tau_exp_smpl=True,
    taumin=[0.01, 1, 5, 25],
    taumax=10000,
    n_tausmpls=100,
    taustep=5,
    nanglesteps=24
)
idxs = [idx[0] for idx in idxs]
fig, ax = plot_results(
    files=itemgetter(*idxs)(files), 
    labelparams=['taumin'], 
    legtitle='tau_min [fs]', 
)
plt.tight_layout()
for ax, letter in zip(ax, ['a','b','c']):
    ax.text(-0.17, 1.05, letter, fontsize=10, fontweight='bold', color='black', transform=ax.transAxes)
fig.suptitle("Variation of lower integration limit", fontsize=fs, y=0.99)
fig.savefig('figures/figure2abc.png', dpi=300)
plt.show()
idxs = get_idxs(
    rhoexpmin=1,
    n_rhosmpls=20,
    tau_exp_smpl=True,
    taumin=5,
    taumax=10000,
    n_tausmpls=100,
    taustep=5,
    nanglesteps=[24,48]
)
idxs = [idx[0] for idx in idxs]
fig, ax = plot_results(
    files=itemgetter(*idxs)(files), 
    labelparams=['nstepstheta'], 
    legtitle='steps'
)
for ax, letter in zip(ax, ['a','b','c']):
    ax.text(-0.17, 1.05, letter, fontsize=10, fontweight='bold', color='black', transform=ax.transAxes)
plt.tight_layout()
fig.savefig('figures/figure6.png', dpi=300)
idxs = get_idxs(
    rhoexpmin=1,
    n_rhosmpls=[20,50],
    tau_exp_smpl=True,
    taumin=5,
    taumax=10000,
    n_tausmpls=100,
    taustep=5,
    nanglesteps=24
)
idxs = [idx[0] for idx in idxs]
fig, ax = plot_results(
    files=itemgetter(*idxs)(files), 
    labelparams=['n_rhosmpls'], 
    legtitle='samples'
)
plt.tight_layout()
for ax, letter in zip(ax, ['d','e','f']):
    ax.text(-0.17, 1.05, letter, fontsize=10, fontweight='bold', color='black', transform=ax.transAxes)
fig.savefig('figures/figure5def.png', dpi=300)
idxs = get_idxs(
    rhoexpmin=[0.01, 1],
    n_rhosmpls=20,
    tau_exp_smpl=True,
    taumin=5,
    taumax=10000,
    n_tausmpls=100,
    taustep=5,
    nanglesteps=24
)
idxs = [idx[0] for idx in idxs]
fig, ax = plot_results(
    files=itemgetter(*idxs)(files), 
    labelparams=['rhoexpmin'], 
    legtitle='min of\nlog-sampling [µm]'
)
plt.tight_layout()
for ax, letter in zip(ax, ['a','b','c']):
    ax.text(-0.17, 1.05, letter, fontsize=10, fontweight='bold', color='black', transform=ax.transAxes)
fig.savefig('figures/figure5abc.png', dpi=300)
