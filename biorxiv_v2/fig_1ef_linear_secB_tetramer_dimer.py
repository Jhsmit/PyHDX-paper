import pandas as pd
from matplotlib.axes import Axes
from pyhdx.fileIO import csv_to_protein

from functions.base import settings_dict
from functions.formatting import *
from functions.logging import write_log

write_log(__file__)


output = 'save'
output_dir = current_dir / 'figures' / 'Fig_1'
output_dir.mkdir(parents=True, exist_ok=True)
plt.savefig(output_dir / 'Fig_ecsecb_deltaG.pdf')
plt.savefig(output_dir / 'Fig_ecsecb_deltaG.png')

fname_tetramer = 'Fig_1e_ecsecb_tetramer'
fname_tetramer_dimer = 'Fig_1f_ecsecb_tetramer_dimer'

fit_kwargs = settings_dict['ecsecb_tetramer_dimer']

gray_value = 235/255
facecolor = [gray_value]*3

figure_width = 75
figure_width /= 25.4

hspace = 0.1
aspect = 4


def make_graph(axes, protein, cmap, norm, ylabel, tick_labels, **ax_format):
    yvals = protein['deltaG'] * 1e-3
    rgba_colors = cmap(norm(yvals), bytes=True)
    hex_colors = rgb_to_hex(rgba_colors)

    ax = axes[0]
    ax.errorbar(protein.index, yvals, yerr=protein['covariance'] * 1e-3, **errorbar_kwargs, zorder=-1)
    ax.scatter(protein.index, yvals, c=hex_colors, **scatter_kwargs)

    ax.format(ylabel='', xlabel=r_xlabel, **ax_format)
    #@ax.format(xlim=(0, None), ylim=(45, 5))
    add_colorbar(ax, ax, cmap, norm, tick_labels=tick_labels, label=ylabel)

    # Add light grey background such that white datapoints can be seen
    xlim = axes[0].get_xlim()
    ylim = axes[0].get_ylim()
    rect = mpl.patches.Rectangle((xlim[0], ylim[0]), xlim[1] - xlim[0], ylim[1] - ylim[0],
                                 facecolor=facecolor, edgecolor='none', zorder=-5)
    axes[0].add_patch(rect)

    # Add grey background for regions of no coverage beyond covered range
    xlim = axes[1].get_xlim()
    ylim = axes[1].get_ylim()
    rect = mpl.patches.Rectangle((ylim[0], xlim[0]), xlim[1]-xlim[0], ylim[1]-ylim[0],
                                 facecolor=no_coverage, edgecolor='none', zorder=-5)
    axes[1].add_patch(rect)

    img = np.expand_dims(yvals, 0)
    extent = [protein.index.min()-0.5, protein.index.max()+0.5, 0, 1]
    axes[1].format(yticks=[])
    Axes.imshow(axes[1], norm(img), aspect='auto', cmap=cmap, vmin=0, vmax=1, interpolation='None', extent=extent)



#-------------------------------------------------------------------------#
# ecSecB tetramer Single
#-------------------------------------------------------------------------#

state = 'ecSecB'
fpath_protein = current_dir / 'fits' / state / f"r1_{fit_kwargs['r1']}" / 'fit_result.csv'
protein = csv_to_protein(fpath_protein)
fig, axes = pplt.subplots(nrows=2, axwidth=figure_width, aspect=aspect, hratios=[5, 1], hspace=hspace)
make_graph(axes, protein, rgb_cmap, rgb_norm, dG_ylabel, tick_labels=[10, 20, 30, 40], xlim=(0, None), ylim=(45, 5))

if output == 'show':
    plt.show()
elif output == 'save':
    plt.savefig(output_dir / f'{fname_tetramer}.pdf')
    plt.savefig(output_dir / f'{fname_tetramer}.png')

#-------------------------------------------------------------------------#
# ecSecB tetramer/dimer delta
#-------------------------------------------------------------------------#

#
names = ['ecSecB', 'ecSecB_dimer']
states = ['SecB WT apo', 'SecB his dimer apo']

fit_output_file = current_dir / 'batch_fits' / 'SecB_tetramer_dimer' / f"r2_{fit_kwargs['r2']}" / 'fit_result.csv'
protein = csv_to_protein(fit_output_file)

p1, p2 = protein[states[0]], protein[states[1]]
dG = p1['deltaG'] - p2['deltaG']
cov = np.sqrt(p1['covariance'] ** 2 + p2['covariance'])

protein = pd.concat([dG, cov], axis=1, keys=['deltaG', 'covariance'])

fig, axes = pplt.subplots(nrows=2, axwidth=figure_width, aspect=aspect, hratios=[5, 1], hspace=hspace)
make_graph(axes, protein, diff_cmap, diff_norm, ddG_ylabel, tick_labels=[-10, 0, 10], xlim=(0, None), ylim=(15, -15))

if output == 'show':
    plt.show()
elif output == 'save':
    plt.savefig(output_dir / f'{fname_tetramer_dimer}.png')
    plt.savefig(output_dir / f'{fname_tetramer_dimer}.pdf')
