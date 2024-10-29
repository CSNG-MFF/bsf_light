from load_save_utils import load_pickle, load_yaml
from load_original import load_matlab_model_data, load_published_model_data
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
from utils import mirror_x_axis
import numpy as np

# n_samples for simulations along radial or depth (z) direction:
comp_time_nsmps = dict()
comp_time_nsmps['laptop'] = [400/5, 700/5, 1000/5]
comp_time_nsmps['cluster'] = [400/5, 700/5, 1000/5, 2000/5, 4000/5]
# dictionary with computation times:
comp_time_rad = dict()
comp_time_rad['laptop'] = list()
comp_time_rad['cluster'] = list()
comp_time_z = dict()
comp_time_z['laptop'] = list()
comp_time_z['cluster'] = list()
for name in ['low_vol_radial', 'default', 'high_vol_radial']:
    comp_time_rad['laptop'].append(
            load_pickle('results/'+name+'.pickle')['comp_time_s']
    )
for name in ['low_vol_z', 'normal_vol_z', 'high_vol_z']:
    comp_time_z['laptop'].append(
            load_pickle('results/'+name+'.pickle')['comp_time_s']
    )
for name in ['low_vol_radial', 'default', 'high_vol_radial', 'very_high_vol_radial', 'very_very_high_vol_radial']:
    comp_time_rad['cluster'].append(
            load_pickle('results/'+name+'_cluster.pickle')['comp_time_s']
    )
for name in ['low_vol_z', 'normal_vol_z', 'high_vol_z', 'very_high_vol_z', 'very_very_high_vol_z']:
    comp_time_z['cluster'].append(
            load_pickle('results/'+name+'_cluster.pickle')['comp_time_s']
    )

# manually stopped matlab computation time:
comp_time_matlab = 28 # sec.
orig_matlab = load_matlab_model_data()
n_smps_matlab_radial, n_smps_matlab_z = np.shape(orig_matlab['zz'])

fs = 8
A4_w, A4_h = 8.27, 11.69 # inch
fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(A4_w*0.8, A4_h*0.2))
axs = axs.flatten()

for handle, ls in zip(['laptop','cluster'], ['solid', 'dotted']):
    axs[0].plot(
        comp_time_nsmps[handle],
        comp_time_rad[handle],
        ls=ls,
        marker='.',
        label=handle
    )
    axs[1].plot(
        comp_time_nsmps[handle],
        comp_time_z[handle],
        ls=ls,
        marker='.',
        label=handle
    )
axs[0].plot(
        [n_smps_matlab_radial],
        [comp_time_matlab],
        marker='x',
        ls=None,
        color='tab:red'
)
axs[0].annotate(
        'Original model', xy=(n_smps_matlab_radial, comp_time_matlab), 
        xytext=(n_smps_matlab_radial + 10, comp_time_matlab + 1),
        color='black',
        fontsize=fs
) 
axs[1].plot(
        [n_smps_matlab_z],
        [comp_time_matlab],
        marker='x',
        ls=None,
        color='tab:red'
)
axs[1].annotate(
        'Original model', xy=(n_smps_matlab_z, comp_time_matlab), 
        xytext=(n_smps_matlab_z + 10, comp_time_matlab + 1),
        color='black',
        fontsize=fs
) 
axs[1].legend(fontsize=fs, loc='upper right', borderaxespad=0., frameon=False)
axs[0].set_ylabel('Computing time [s]', fontsize=fs)
axs[0].set_xlabel('# radial samples', fontsize=fs)
axs[1].set_xlabel('# depth samples', fontsize=fs)

for ax in axs:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

# insert letters:
for ax, letter in zip(axs[:], ['a', 'b']):
    ax.text(-0.1, 1.07, letter, fontsize=10, fontweight='bold', color='black', transform=ax.transAxes)
plt.tight_layout()
fig.savefig('figures/figure4.png', dpi=300)
plt.close()
