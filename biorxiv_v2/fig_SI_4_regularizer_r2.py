import re

from pyhdx.fileIO import csv_to_protein

from functions.base import *
from functions.formatting import *
from functions.logging import write_log

write_log(__file__)


output = 'save'
output_dir = current_dir / 'figures' / 'Supplement'
output_dir.mkdir(parents=True, exist_ok=True)
fname = 'Fig_SI_4_influence_r2'

fit_kwargs = settings_dict['ecsecb_tetramer_dimer']
fit_folder = 'SecB_tetramer_dimer'

states = ['SecB WT apo', 'SecB his dimer apo']

gray_value = 220/255
facecolor = [gray_value]*3
figure_width = 183  # Width in mm
figure_width /= 25.4

fig, axes = pplt.subplots(nrows=4, ncols=2, width=figure_width, aspect=3, sharex=False)
for ax, r2 in zip(axes, [0, 0.05, 0.1, 0.5, 1, 2, 5, 10]):
    fit_output_file = current_dir / 'batch_fits' / fit_folder / f"r2_{r2}" / 'fit_result.csv'
    protein = csv_to_protein(fit_output_file)

    p1, p2 = protein[states[0]], protein[states[1]]
    dG = p1['deltaG'] - p2['deltaG']
    yvals = dG *1e-3
    cov = np.sqrt(p1['covariance'] ** 2 + p2['covariance'])

    rgba_colors = diff_cmap(diff_norm(yvals), bytes=True)
    hex_colors = rgb_to_hex(rgba_colors)

    ax.errorbar(protein.index, yvals, yerr=cov * 1e-3, **errorbar_kwargs, zorder=-1)
    ax.scatter(protein.index, yvals, c=hex_colors, **scatter_kwargs)

    logfile =  current_dir / 'batch_fits' / fit_folder / f"r2_{r2}" / 'log.txt'
    line = logfile.read_text().split('\n')[2]
    result = re.search('\((.*?)\%\)', line)
    value = result[0][1:-1]

    ax.format(xlim=(0, None), ylim=(-15, 15), title=f'λ$_{2}$ = {r2}, Regularizer loss: {value}')
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    rect = mpl.patches.Rectangle((xlim[0], ylim[0]), xlim[1] - xlim[0], ylim[1] - ylim[0],
                                 facecolor=facecolor, edgecolor='none', zorder=-5)
    ax.add_patch(rect)

ticks = np.array([-10, 0, 10])
cbar = fig.colorbar(diff_cmap, width=cbar_width, ticks=diff_norm(ticks))
cbar.set_label(ddG_ylabel)
cbar.ax.set_yticklabels(ticks)

axes.format(xlabel=r_xlabel, ylabel=ddG_ylabel)


if output == 'show':
    plt.show()
elif output == 'save':
    plt.savefig(output_dir / f'{fname}.png')
    plt.savefig(output_dir / f'{fname}.pdf')

