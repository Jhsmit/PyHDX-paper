import re

from pyhdx.fileIO import csv_to_protein

from functions.formatting import *
from functions.logging import write_log

write_log(__file__)

output = 'save'
output_dir = current_dir / 'figures' / 'Supplement'
output_dir.mkdir(parents=True, exist_ok=True)
fname = 'Fig_SI_3_influence_r1'

state = 'ecSecB'
fit_dir = 'fits'

figure_width = 183  # Width in mm
figure_width /= 25.4

fig, axes = pplt.subplots(nrows=4, ncols=2, width=figure_width, aspect=3, sharex=False)
for ax, r1 in zip(axes, [0, 0.05, 0.1, 0.5, 1, 2, 5, 10]):
    reg_dir = f"r1_{r1}"
    protein = csv_to_protein(current_dir / fit_dir / state / reg_dir / 'fit_result.csv')

    yvals = protein['deltaG'] * 1e-3
    rgba_colors = rgb_cmap(rgb_norm(yvals), bytes=True)
    hex_colors = rgb_to_hex(rgba_colors)

    ax.errorbar(protein.index, yvals, yerr=protein['covariance'] * 1e-3, **errorbar_kwargs, zorder=-1)
    ax.scatter(protein.index, yvals, c=hex_colors, **scatter_kwargs)
    ax.format(ylabel='', ylim=(40, 0))

    logfile = current_dir / 'fits' / state / f'r1_{r1}' / 'log.txt'
    line = logfile.read_text().split('\n')[2]
    result = re.search('\((.*?)\%\)', line)
    value = result[0][1:-1]

    ax.format(xlim=(0, None), ylim=(5, 50), title=f'λ$_{1}$ = {r1}, Regularizer loss: {value}')

ticks = pplt.arange(10, 40, 5)
cbar = fig.colorbar(rgb_cmap, width=cbar_width, ticks=rgb_norm(ticks))
cbar.set_label(dG_ylabel)
cbar.ax.set_yticklabels(ticks)

axes.format(xlabel=r_xlabel, ylabel=dG_ylabel)

if output == 'show':
    plt.show()
elif output == 'save':
    plt.savefig(output_dir / f'{fname}.png')
    plt.savefig(output_dir / f'{fname}.pdf')



