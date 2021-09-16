
from functions.formatting import *
from functions.base import fitresults_dir, settings_dict
from functions.logging import write_log
from pyhdx.fileIO import csv_to_hdxm, csv_to_protein
from pyhdx.plot import plot_peptides

#%%
output = 'save'
output_dir = current_dir / 'figures' / 'Supplement'
output_dir.mkdir(parents=True, exist_ok=True)
fname_reference = 'Fig_SI_5_reference'
fname = 'Fig_SI_5_reduced_redundancy'

fit_dir = 'fits_ecsecb_reduced'
figure_width = page_width  # Width in mm
figure_width /= 25.4

#%%

f = 'ecSecB_r1'
fit_kwargs = settings_dict['ecsecb_tetramer_dimer']

fpath_protein = fitresults_dir / f / f"r1_{fit_kwargs['r1']}" / 'fit_result.csv'
protein = csv_to_protein(fpath_protein)
hdxm = csv_to_hdxm(fitresults_dir / f / f"r1_{fit_kwargs['r1']}" / 'HDXMeasurement.csv')

fig, axes = pplt.subplots(ncols=1, nrows=2, aspect=3.5, hspace=0.1,
                          sharey=False, axwidth=49.3/25.4)

plot_peptides(hdxm[1], axes[0])
cmap = mpl.cm.get_cmap('jet')
cbar = axes[0].colorbar(cmap, width=cbar_width, ticks=[0, 1], space=0)
cbar.ax.set_yticklabels([0, 1])
cbar.set_label('RFU', labelpad=-2)


yvals = protein['deltaG'] * 1e-3
rgba_colors = rgb_cmap(rgb_norm(yvals), bytes=True)
hex_colors = rgb_to_hex(rgba_colors)

axes[1].errorbar(protein.index, yvals, yerr=protein['covariance'] * 1e-3, **errorbar_kwargs, zorder=-1)
axes[1].scatter(protein.index, yvals, c=hex_colors, **scatter_kwargs)
axes[1].format(ylabel='', ylim=(45, 5), xlabel=r_xlabel, xlim=(0, 160))

axes.format(toplabels=[f"Np: {hdxm.Np}"])

add_colorbar(axes[1], axes[1], rgb_cmap, rgb_norm, tick_labels=[20, 40], label=dG_ylabel)

if output == 'show':
    plt.show()
elif output == 'save':
    plt.savefig(output_dir / f'{fname_reference}.png')
    plt.savefig(output_dir / f'{fname_reference}.pdf')


fig, axes = pplt.subplots(ncols=3, nrows=6, aspect=3.5, hspace=(0, None, 0, None, 0),
                          sharey=False, width=figure_width)
fracs = [0.75, 0.5, 0.25]
Np = 63
for i, ax in enumerate(axes):
    col = i % 3
    row = i // 3
    repeat = row // 2
    input_dir = fitresults_dir / fit_dir / f'frac_{fracs[col]}_{repeat}'

    hdxm = csv_to_hdxm(input_dir / 'hdxm.csv')
    print(hdxm.Np)
    protein = csv_to_protein(input_dir / 'fit_result.csv')
    if row % 2 == 0:
        plot_peptides(hdxm[1], ax)
    else:
        yvals = protein['deltaG'] * 1e-3
        rgba_colors = rgb_cmap(rgb_norm(yvals), bytes=True)
        hex_colors = rgb_to_hex(rgba_colors)

        ax.errorbar(protein.index, yvals, yerr=protein['covariance'] * 1e-3, **errorbar_kwargs, zorder=-1)
        ax.scatter(protein.index, yvals, c=hex_colors, **scatter_kwargs)
        if col == 0:
            ax.format(ylabel=dG_ylabel, ylim=(45, 5), xlabel=r_xlabel, xlim=(0, 160))
        else:
            ax.format(ylabel='', ylim=(45, 5), xlabel=r_xlabel, xlim=(0, 160))


toplabels = [f"{'abc'[col]}) Np: {int(np.round(fracs[col]*Np))}" for col in range(3)]
numerals = ['I'*(i+1) for i in range(3)]
fill = ['']*3
rightlabels = [val for pair in zip(numerals, fill) for val in pair]
axes.format(toplabels=toplabels, rightlabels=rightlabels)

if output == 'show':
    plt.show()
elif output == 'save':
    #bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    #width, height = bbox.width, bbox.height
    #print(width, height)  # get width for reference plot
    plt.savefig(output_dir / f'{fname}.png')
    plt.savefig(output_dir / f'{fname}.pdf')

