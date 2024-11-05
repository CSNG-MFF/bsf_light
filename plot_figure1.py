from load_save_utils import load_pickle, load_yaml
from load_original import load_published_exp_data, load_published_model_data
import matplotlib.pyplot as plt

replication_params = load_yaml('params/default.yml')
replication = load_pickle('results/default.pickle')['final']
orig_exp_depth, orig_exp_radial = load_published_exp_data()
orig_mod_depth, orig_mod_radial = load_published_model_data()

fs = 8
A4_w, A4_h = 8.27, 11.69 # inch
fig, axs = plt.subplots(ncols=3, figsize=(A4_w*0.7, A4_h*0.2))

   
# tranmission over depth
axs[0].errorbar(
    orig_exp_depth['z'],
    orig_exp_depth['transmission'],
    yerr=orig_exp_depth['transmission_err'],
    label='Experiment', 
    marker='x',
    color='black', 
    linestyle='',
    capsize=5
)
axs[0].plot(
    replication['z'][0, :],
    replication['combined'][0, :],
    ls='solid', 
    color='orange', 
    label='Replication'
)
axs[0].set_xlim(0,700)
axs[0].set_ylim(0,0.33)
axs[0].set_ylabel('Normalized intensity', fontsize=fs)
axs[0].set_xlabel('Depth [µm]', fontsize=fs)


# transmission over radial distance
axs[1].set_title('Depth = 300 µm', fontsize=fs)
axs[1].plot(
    orig_exp_radial['x_z300'],
    orig_exp_radial['transmission_z300'],
    marker='x',c='black', linestyle='',
    label='Experiment', 
)
axs[1].plot(
    orig_mod_radial['x_z300'],
    orig_mod_radial['transmission_z300'],
    ls='--', color='blue', 
    label='Original model'
)
z300 = int(300/replication_params['dz'])
axs[1].plot(
    replication['rho'][:, z300], 
    replication['combined'][:, z300]/replication['combined'][0, z300],
    ls='solid', color='orange', 
    label='Replication'
)
axs[2].set_title('depth = 600 µm', fontsize=fs)
axs[2].plot(
    orig_exp_radial['x_z600'],
    orig_exp_radial['transmission_z600'],
    marker='x',c='black', linestyle='',
    label='Experiment', 
)
axs[2].plot(
    orig_mod_radial['x_z600'],
    orig_mod_radial['transmission_z600'],
    ls='--', color='blue', 
    label='Original model'
)
z600 = int(600/replication_params['dz'])
axs[2].plot(
    replication['rho'][:, z600], 
    replication['combined'][:, z600]/replication['combined'][0, z600],
    ls='solid', color='orange', 
    label='Replication'
)

for ax in axs:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

for ax in axs[1:]:
    ax.set_xlim(0,500)
    ax.set_ylim(0,1.1)
    ax.set_xlabel('Radial distance [µm]', fontsize=fs)


for ax in axs:
    ax.tick_params(axis='both', which='major', labelsize=fs)
axs[-1].legend(fontsize=fs, bbox_to_anchor=(0.8, 0.82), loc='center',frameon=False)

# insert letters:
for ax, letter in zip(axs, ['a','b', 'c']):
    ax.text(-0.1, 1.05, letter, fontsize=10, fontweight='bold', color='black', transform=ax.transAxes)
plt.tight_layout()
fig.savefig('figures/figure1.png', dpi=300)
plt.close()
