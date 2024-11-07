from load_save_utils import load_pickle, load_yaml
from load_original import load_published_exp_data, load_published_model_data
import matplotlib.pyplot as plt

labels = ['Approximation by\nLutomirski et al.', 
          'Approximation by\nvan de Hulst and\nKattawar', 
          'Prec. calc. of µ\nand approximate\nrelation of σ to µ']
linestyles = ['dotted', 'dashed', 'solid']
colors = ['tab:orange','tab:blue', 'tab:green']
replication_paramss = [
    load_yaml(f) for f in ['params/default_Lutomirski.yml','params/default_vandeHulst.yml','params/default.yml']
]
replications = [
    load_pickle(f)['final'] for f in [
        'results/default_Lutomirski.pickle','results/default_vandeHulst.pickle','results/default.pickle']
]

orig_exp_depth, orig_exp_radial = load_published_exp_data()
orig_mod_depth, orig_mod_radial = load_published_model_data()

fs = 8
A4_w, A4_h = 8.27, 11.69 # inch
fig, axs = plt.subplots(ncols=2, nrows=2, figsize=(A4_w*0.7, A4_h*0.4))
axs = axs.flatten()
axs[1].axis('off')
axs = [axs[0], axs[2], axs[3]]
#axs[3].axis('off')
# model
# tranmission over depth
for replication, replication_params, label, ls, color in zip(replications, replication_paramss, labels, linestyles, colors):
    axs[0].plot(
        replication['z'][1, :],
        replication['combined'][1, :],
        ls=ls, 
        color=color, 
    )
    # transmission over radial distance
    z300 = int(300/replication_params['dz'])
    axs[1].plot(
        replication['rho'][:, z300], 
        replication['combined'][:, z300]/replication['combined'][0, z300],
        ls=ls, color=color, 
    )
    z600 = int(600/replication_params['dz'])
    axs[2].plot(
        replication['rho'][:, z600], 
        replication['combined'][:, z600]/replication['combined'][0, z600],
        ls=ls, color=color, 
        label=label
    )

# experimental
# tranmission over depth
axs[0].errorbar(
    orig_exp_depth['z'],
    orig_exp_depth['transmission'],
    yerr=orig_exp_depth['transmission_err'],
    marker='x',
    color='black', 
    linestyle='',
    capsize=5
)
# transmission over radial distance
axs[1].plot(
    orig_exp_radial['x_z300'],
    orig_exp_radial['transmission_z300'],
    marker='x',c='black', linestyle='',
)
axs[1].plot(
    orig_mod_radial['x_z300'],
    orig_mod_radial['transmission_z300'],
    ls='solid', color='black', 
)
axs[2].plot(
    orig_exp_radial['x_z600'],
    orig_exp_radial['transmission_z600'],
    marker='x',c='black', linestyle='',
    label='Experiment', 
)
axs[2].plot(
    orig_mod_radial['x_z600'],
    orig_mod_radial['transmission_z600'],
    ls='solid', color='black',
    label='Original model'
)

# ax properties
axs[0].set_xlim(0,700)
axs[0].set_ylim(0,0.33)
axs[0].set_ylabel('Transmission', fontsize=fs)
axs[0].set_xlabel('Depth [µm]', fontsize=fs)

axs[1].set_title('Depth = 300 µm', fontsize=fs)
axs[2].set_title('Depth = 600 µm', fontsize=fs)
for ax in axs:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

axs[1].set_ylabel('Normalized transmission', fontsize=fs)
for ax in axs[1:]:
    ax.set_xlim(0,500)
    ax.set_ylim(0,1.1)
    ax.set_xlabel('Radial distance [µm]', fontsize=fs)


for ax in axs:
    ax.tick_params(axis='both', which='major', labelsize=fs)
fig.legend(loc="upper right", bbox_to_anchor=(0.9, .95), frameon=False, fontsize=fs)


# insert letters:
for ax, letter in zip(axs, ['a','b', 'c']):
    ax.text(-0.2, 1.05, letter, fontsize=10, fontweight='bold', color='black', transform=ax.transAxes)
plt.tight_layout()
fig.savefig('figures/figure1.png', dpi=300)
plt.close()
