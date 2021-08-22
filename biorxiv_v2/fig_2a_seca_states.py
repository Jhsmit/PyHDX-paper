import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import proplot as pplt
from matplotlib.axes import Axes
from pyhdx.fileIO import csv_to_protein

from functions.base import *
from functions.formatting import *
from functions.logging import write_log

write_log(__file__)


output = 'save'
output_dir = current_dir / 'figures' / 'Fig_2'
output_dir.mkdir(parents=True, exist_ok=True)
fname = 'Fig_2a_SecA_linear_bars'

figure_width = 90  # Width in mm
figure_width /= 25.4


states = [
    'SecA-monomer',
    'SecA mono ADP',
    'SecA1-901 wt apo',
    'SecA wt ADP',
    'SecA1-834 apo',
    'SecA 1-834 ADP'
    ]

labels = [
    'SecA$_{1}$',
    'Monomer ADP',
    'SecA$_{2}$',
    'Dimer ADP',
    'SecA-ΔC$_{2}$',
    'ΔCtail ADP',
]

# Itervals are inclusive, inclusive
structural_motifs = {
    'Q': (80, 87),
    'I': (104, 109),
    'Ia': (130, 136),
    'Ib': (159, 160),
    'Ic': (175, 179),
    'II': (207, 214),
    'III': (391, 393),
    'IV': (453, 458),
    'IVa': (484, 487),
    '': (491, 491),
    'V': (500, 504),
    'Va': (508, 511),
    'Vb': (556, 558),
    'VI': (570, 577),
}


fit_kwargs = settings_dict['seca_all_states']


protein = csv_to_protein(current_dir / 'batch_fits' / 'SecA_states' / f'fit_result.csv')

res = 20
p = 0.90
#%%

array = np.zeros((len(states), res))
for i, row in enumerate(array):
    row[:] = i + 1
    row[-1] = 0

fig, axes = pplt.subplots(array, aspect=20, width=figure_width,
                          hspace=(0, None, 0, None, 0), top=0.25, left='2em')

extent = [protein.index.min()-0.5, protein.index.max()+0.5, 0, 1]
for i, (ax, state, label) in enumerate(zip(axes, states, labels)):
    if i % 2 == 0:
        yvals = protein[state]['deltaG']*1e-3
        img = np.expand_dims(yvals, 0)

        Axes.imshow(ax, rgb_norm(img), aspect='auto', cmap=rgb_cmap, vmin=0, vmax=1, interpolation='None', extent=extent)
        ax.format(yticks=[])
        Axes.text(ax, 1.02, 0.5, label, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
        Axes.text(ax, -0.01, 0.2, 'ΔG', horizontalalignment='right', verticalalignment='bottom',
                  transform=ax.transAxes, fontsize=6)

    else:
        diff_vals = (protein[state]['deltaG'] - protein[states[i - 1]]['deltaG'])*1e-3
        img = np.expand_dims(diff_vals, 0)

        Axes.imshow(ax, diff_norm(img), aspect='auto', cmap=diff_cmap, vmin=0, vmax=1, interpolation='None', extent=extent)
        ax.format(yticks=[])
        Axes.text(ax, 1.02, 0.5, '+ ADP', horizontalalignment='left',
          verticalalignment='center', transform=ax.transAxes)
        Axes.text(ax, -0.01, 0.2, 'ΔΔG', horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes,
                  fontsize=6)


axes.format(xlabel=r_xlabel, xlim=(0, 950))

top_ax = axes[0]
split_annotations = ['I', 'IV', 'V']
for i, (k, v) in enumerate(structural_motifs.items()):
    x1, x2 = v

    x1 -= 0.5
    x2 += 0.5
    y = 1.2
    for ax in axes[::2]:
        line = mpl.lines.Line2D([x1, x2], [y, y], color='k', lw=1.5, solid_capstyle='butt')
        line.set_clip_on(False)
        ax.add_line(line)

    xmid = np.mean([x1, x2-2])
    ytext = 1.4

    motif = k.strip('abc')
    if motif in split_annotations:
        label = k[-1] if k.endswith(('a', 'b', 'c')) else ''
    else:
        label = motif
#        ytext += 0.7

    trans = top_ax.get_xaxis_transform()
    if label == 'VI':
        xmid += 8  # Nudge last label
    top_ax.text(xmid, ytext, label, transform=trans, horizontalalignment='center', fontsize=7)


for s in split_annotations:
    elements = [(a, b) for k, (a, b) in structural_motifs.items() if k.strip('abc') == s]
    rmin = elements[0][0]
    rmax = elements[-1][1]
    rmin -= 0.5
    rmax += 0.5
    y = 1.2

    #top_ax.arrow(rmin, y, rmax-rmin, 0, head_width=0)
    y = 2
    line = mpl.lines.Line2D([rmin, rmax], [y, y], color='k', lw=0.5, solid_capstyle='butt')
    line.set_clip_on(False)
    top_ax.add_line(line)

    cap_size = 0.2
    cap1 = mpl.lines.Line2D([rmin, rmin], [y-cap_size, y+cap_size], color='k', lw=0.3)
    cap1.set_clip_on(False)
    top_ax.add_line(cap1)

    cap2 = mpl.lines.Line2D([rmax, rmax], [y-cap_size, y+cap_size], color='k', lw=0.3)
    cap2.set_clip_on(False)
    top_ax.add_line(cap2)

    xmid = (rmin + rmax)/2
    trans = top_ax.get_xaxis_transform()
    top_ax.text(xmid, 2.1, s, horizontalalignment='center', fontsize=7)

if output == 'show':
    plt.show()
elif output == 'save':
    plt.savefig(output_dir / f'{fname}.png')
    plt.savefig(output_dir / f'{fname}.pdf')

