from functions.logging import write_log
from functions.base import *
from functions.rainbows import label_axes, kdeplot, boxplot, stripplot
from functions.formatting import *
import proplot as pplt
import matplotlib.pyplot as plt
import numpy as np
from pyhdx.fileIO import csv_to_protein

# Rainclouds:
# Allen M, Poggiali D, Whitaker K et al. Raincloud plots: a multi-platform tool for robust data visualization [version 1; peer review: 2 approved]. Wellcome Open Res 2019, 4:63. DOI: 10.12688/wellcomeopenres.15191.1

write_log(__file__)

output = 'save'
output_dir = current_dir / 'figures' / 'Fig_2'
output_dir.mkdir(parents=True, exist_ok=True)
fname = 'Fig_2fg_rainbow_proteins'


reg = 1
plt.rcParams["image.composite_image"] = False

f_states = [
    'TF_monomer_4C',
    'SecA_monomer',
    'MBP_wt',
    'PPiA_WT',
    'PPiB_WT',
    'ecSecB',
    'hPREP_default',
    'bcl2',
    'EscV'
]

f_labels = [
    '$\it{ec}$TF',
    '$\it{ec}$SecA',
    '$\it{ec}$MBP',
    '$\it{ec}$PpiA',
    '$\it{ec}$PpiB',
    '$\it{ec}$SecB',
    '$\it{h}$PREP',
    '$\it{h}$Bcl-2',
    '$\it{ec}$SctV'
]

f_labels = np.array(f_labels)

# Load data
fit_dir = 'single_fits'
reg_dir = f"r1_{reg}"
files = [fitresults_dir / fit_dir / state / reg_dir / 'fit_result.csv' for state in f_states]

# List of deltaG values as numpy arrays
deltaGs = np.array([csv_to_protein(f)['deltaG'].dropna().to_numpy()*1e-3 for f in files], dtype=object)

# Sort by mean deltaG
idx = np.argsort([np.mean(arr) for arr in deltaGs])
f_data = deltaGs[idx]
f_labels = f_labels[idx]

g_states = ['TF_dimer_4C', 'SecA_1-901_wt_apo']
g_labels = np.array(['$\it{ec}$TF$_{2}$', '$\it{ec}$SecA$_{2}$'])

files = [fitresults_dir / fit_dir / state / reg_dir / 'fit_result.csv' for state in g_states]
deltaGs = np.array([csv_to_protein(f)['deltaG'].dropna().to_numpy()*1e-3 for f in files], dtype=object)

idx = np.argsort([np.mean(arr) for arr in deltaGs])
g_data = deltaGs[idx]
g_labels = g_labels[idx]

#%%

boxplot_width = 0.1
orientation = 'vertical'

strip_kwargs = dict(offset=0.0, orientation=orientation, s=2, colors='k', jitter=0.2, alpha=0.25)
kde_kwargs = dict(linecolor='k', offset=0.15, orientation=orientation, fillcolor=False, fill_cmap=rgb_cmap,
        fill_norm=rgb_norm, y_scale=5, y_norm=None, linewidth=1)
boxplot_kwargs = dict(offset=0.2, sym='', linewidth=1., linecolor='k', orientation=orientation, widths=boxplot_width)

fig, axes = pplt.subplots(ncols=2, width=page_width/25.4, aspect=5, ref=1, wratios=[4, 1], sharey=False)
ax = axes[0]
stripplot(f_data, ax=ax, **strip_kwargs)
kdeplot(f_data, ax=ax, **kde_kwargs)
boxplot(f_data, ax=ax, **boxplot_kwargs)
label_axes(f_labels, ax=ax, rotation=0)
ax.format(xlim=(-0.75, len(f_data) - 0.5), ylabel=dG_ylabel, yticklabelloc='left', ytickloc='left',
          ylim=ax.get_ylim()[::-1])


ax = axes[1]
stripplot(g_data, ax=ax, **strip_kwargs)
kdeplot(g_data, ax=ax, **kde_kwargs)
boxplot(g_data, ax=ax, **boxplot_kwargs)
#
label_axes(g_labels, ax=ax, rotation=0)
ax.format(xlim=(-0.5, len(g_data) - 0.5), yticklabelloc='right', ytickloc='right', ylabel='',
          ylim=axes[0].get_ylim())
#
tick_labels = [0, 20, 40]
add_colorbar(fig, ax, rgb_cmap, rgb_norm, tick_labels=tick_labels)

if output == 'show':
    plt.show()
elif output == 'save':
    plt.savefig(output_dir / f'{fname}.png')
    plt.savefig(output_dir / f'{fname}.pdf')

