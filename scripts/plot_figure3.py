from BSF import load_pickle, load_yaml
from load_original import load_matlab_model_data, load_published_model_data
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
from BSF.utils import mirror_x_axis
import numpy as np

# volumes for simulations along radial or depth (z) direction:
radial_samples = np.array([0.4, 0.7, 1., 2., 4.]) * 0.7 # volume in mm3
depth_samples = np.array([0.4, 0.7, 1., 2., 4.]) * 0.4 # vol in mm3

# dictionary with computation times:
comp_time_rad = []
comp_time_z = []
MaxRS_rad = []
MaxRS_z = []
for name in ['low_vol_radial', 'default', 'high_vol_radial', 'very_high_vol_radial', 'very_very_high_vol_radial']:
    comp_time_rad.append(
            load_pickle('results/'+name+'.pickle')['comp_time_s']
    )
    MaxRS_rad.append(
            load_pickle('results/'+name+'.pickle')['MaxRS_KB']
    )
for name in ['low_vol_z', 'normal_vol_z', 'high_vol_z', 'very_high_vol_z', 'very_very_high_vol_z']:
    comp_time_z.append(
            load_pickle('results/'+name+'.pickle')['comp_time_s']
    )
    MaxRS_z.append(
            load_pickle('results/'+name+'.pickle')['MaxRS_KB']
    )

# manually stopped matlab computation time:
comp_time_matlab = 28 # sec.
orig_matlab = load_matlab_model_data()
matlab_radial_samples, matlab_depth_samples = np.shape(orig_matlab['zz'])
matlab_radial_samples = (matlab_radial_samples * 0.005)**2 * 0.7
matlab_depth_samples = (matlab_depth_samples * 0.005) * 0.4**2


fs = 8
A4_w, A4_h = 8.27, 11.69 # inch
fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(A4_w*0.7, A4_h*0.22))
axs = axs.flatten()

axs[0].plot(
    radial_samples,
    comp_time_rad,
    ls='solid',
    marker='.',
    label='radial'
)
axs[0].plot(
    depth_samples,
    comp_time_z,
    ls='dashed',
    marker='.',
    label='depth'
)
axs[1].plot(
    radial_samples,
    np.array(MaxRS_rad)/1024**2,
    ls='solid',
    marker='.',
    label='radial'
)
axs[1].plot(
    depth_samples,
    np.array(MaxRS_z)/1024**2,
    ls='dashed',
    marker='.',
    label='depth'
)
axs[0].plot(
        [matlab_radial_samples],
        [comp_time_matlab],
        marker='x',
        ls=None,
        color='tab:red'
)
axs[0].annotate(
        'Original model', xy=(matlab_radial_samples, comp_time_matlab), 
        xytext=(matlab_radial_samples + 10, comp_time_matlab + 1),
        color='black',
        fontsize=fs
) 
axs[0].text(0.7, 0.3, r'$\Delta=5\mu$m', transform=axs[0].transAxes, fontsize=fs)
axs[1].legend(fontsize=fs, loc='lower right', borderaxespad=0., frameon=False)
axs[0].set_ylabel('Computing time [s]', fontsize=fs)
axs[1].set_ylabel('Physical memory consumption [GB]', fontsize=fs)
for ax in axs:
    ax.set_xlabel(r'Simulated $\rho$-z-plane [mm²]', fontsize=fs)

for ax in axs:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

# insert letters:
for ax, letter in zip(axs[:], ['a', 'b']):
    ax.text(-0.1, 1.07, letter, fontsize=10, fontweight='bold', color='black', transform=ax.transAxes)
plt.subplots_adjust(bottom=0.2)  # Increase bottom margin if needed
fig.savefig('figures/figure3.png', dpi=300)
plt.close()
