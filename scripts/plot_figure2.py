from BSF import load_pickle, load_yaml
from BSF.load_original import load_matlab_model_data, load_published_model_data
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
from BSF.utils import mirror_x_axis

replication_params = load_yaml('params/default.yml')
replication = load_pickle('results/default.pickle')['final']
replication_err_params = load_yaml('params/error.yml')
replication_err = load_pickle('results/error.pickle')['final']

orig_mod_depth, orig_mod_radial = load_published_model_data()
orig_matlab = load_matlab_model_data()

fs = 8
A4_w, A4_h = 8.27, 11.69 # inch
fig, axs = plt.subplots(ncols=3, nrows=2, figsize=(A4_w*0.7, A4_h*0.32))
axs = axs.flatten()
   
# tranmission over depth
axs[0].plot(
    orig_mod_depth['z'],
    orig_mod_depth['transmission'],
    ls='--', color='blue', 
    label='Original\n(publication)'
)
axs[0].plot(
    orig_matlab['zz'][1, :],
    orig_matlab['data'][1, :],
    ls='solid', 
    color='blue', 
    label='Original\n(matlab app)'
)
axs[0].plot(
    replication['z'][1, :],
    replication['combined'][1, :],
    ls='solid', 
    color='orange', 
    label='Replication'
)
axs[0].plot(
    replication_err['z'][1, :],
    replication_err['combined'][1, :],
    ls='dotted', color='red', 
    label='Replication with\nsmall-volume\nconvolutions'
)
axs[0].set_xlim(0,700)
axs[0].set_ylim(0,0.33)
axs[0].set_ylabel('Transmission', fontsize=fs)
axs[0].set_xlabel('Depth [µm]', fontsize=fs)


# transmission over radial distance
axs[1].set_title('Depth = 300 µm', fontsize=fs)
axs[1].plot(
    orig_mod_radial['x_z300'],
    orig_mod_radial['transmission_z300'],
    ls='--', color='blue', 
    label='Original\n(publication)'
)
axs[1].plot(
    orig_matlab['xx'][:,orig_matlab['z300']],
    orig_matlab['data'][:,orig_matlab['z300']]/orig_matlab['data'][0,orig_matlab['z300']],
    ls='solid', 
    color='blue', 
    label='Original\n(matlab app)'
)
z300 = int(300/replication_params['dz'])
axs[1].plot(
    replication['rho'][:, z300], 
    replication['combined'][:, z300]/replication['combined'][0, z300],
    ls='solid', color='orange', 
    label='Replication'
)
z300 = int(300/replication_err_params['dz'])
axs[1].plot(
    replication_err['rho'][:, z300], 
    replication_err['combined'][:, z300]/replication_err['combined'][0, z300],
    ls='dotted', color='red', 
    label='Replication with\nsmall-volume\nconvolutions'
)
axs[2].set_title('Depth = 600 µm', fontsize=fs)
axs[2].plot(
    orig_mod_radial['x_z600'],
    orig_mod_radial['transmission_z600'],
    ls='--', color='blue', 
    label='Original\n(publication)'
)
axs[2].plot(
    orig_matlab['xx'][:,orig_matlab['z600']],
    orig_matlab['data'][:,orig_matlab['z600']]/orig_matlab['data'][0,orig_matlab['z600']],
    ls='solid', 
    color='blue', 
    label='Original\n(matlab app)'
)
z600 = int(600/replication_params['dz'])
axs[2].plot(
    replication['rho'][:, z600], 
    replication['combined'][:, z600]/replication['combined'][0, z600],
    ls='solid', color='orange', 
    label='Replication'
)
z600 = int(600/replication_err_params['dz'])
axs[2].plot(
    replication_err['rho'][:, z600], 
    replication_err['combined'][:, z600]/replication_err['combined'][0, z600],
    ls='dotted', color='red', 
    label='Replication with\nsmall-volume\nconvolutions'
)

for ax in axs[:3]:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

axs[1].set_ylabel('Normalized transmission', fontsize=fs)
for ax in axs[1:3]:
    ax.set_xlim(0,500)
    ax.set_ylim(0,1.1)
    ax.set_xlabel('Radial distance [µm]', fontsize=fs)

for ax in axs[:3]:
    ax.tick_params(axis='both', which='major', labelsize=fs)
axs[2].legend(fontsize=fs, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., frameon=False)


# 2D profiles
def plot_2D(x, z, data, ax, **kws):
    xm = mirror_x_axis(x, make_neg=True)
    zm = mirror_x_axis(z, make_neg=False)
    datam = mirror_x_axis(data, make_neg=False)
    return ax.pcolormesh(xm,zm,datam, norm=LogNorm(1e-4,1), shading='nearest', cmap='Blues')

axs[3].set_title('Original (matlab app)', fontsize=fs)
plot_2D(orig_matlab['xx'], orig_matlab['zz'], orig_matlab['data'], axs[3])
axs[4].set_title('Replication with\nsmall-volume convolutions', fontsize=fs)
mapp = plot_2D(replication_err['rho'], replication_err['z'], replication_err['combined'], axs[4])
axs[5].set_title('Replication', fontsize=fs)
plot_2D(replication['rho'], replication['z'], replication['combined'], axs[5])


for ax in axs[3:]:
    # Position for arrow 1
    x,y = (-330, 680)
    ax.annotate(' ', xy=(x,y), xytext=(x + 100, y - 100),
                arrowprops=dict(arrowstyle="->", lw=1, color="black"),
                fontsize=fs, color="black")
    
    # Position for arrow 2
    x,y = (-330, 180) 
    ax.annotate(' ', xy=(x,y), xytext=(x + 100, y + 100),
                arrowprops=dict(arrowstyle="->", lw=1, color="black"),
                fontsize=fs, color="black")

    # Create a scale bar
    scalebar = AnchoredSizeBar(ax.transData,
                               100,                # Length of the scale bar in data units
                               '100 µm',           # Label for the scale bar
                               'lower right',      # Location of the scale bar
                               pad=0.5,
                               color='black',      # Color of the scale bar
                               frameon=False,
                               size_vertical=2,    # Thickness of the scale bar
                               fontproperties=fm.FontProperties(size=fs)
    )
    # Add the scale bar to the axis
    ax.add_artist(scalebar)

cax = fig.add_axes([0.84, 0.04, 0.01, 0.32])
cbar = plt.colorbar(mapp, cax=cax)
cbar.set_label('Transmission', fontsize=fs)

for ax in axs[3:]:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('equal')
    ax.set_xlim(-400,400)

# insert letters:
for ax, letter in zip(axs[:3], ['a', 'b', 'c']):
    ax.text(-0.3, 1.07, letter, fontsize=10, fontweight='bold', color='black', transform=ax.transAxes)
for ax, letter in zip(axs[3:],['d', 'e', 'f']):
    ax.text(-0.3, 1.17, letter, fontsize=10, fontweight='bold', color='black', transform=ax.transAxes)
plt.tight_layout()
fig.savefig('figures/figure2.png', dpi=300)
plt.close()
