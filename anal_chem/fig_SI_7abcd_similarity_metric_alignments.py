import matplotlib.pyplot as plt
from pyhdx.fileIO import csv_to_protein
from pyhdx.alignment import align_dataframes
from matplotlib.colors import to_hex

from functions.align import alignments
from functions.base import *
from functions.formatting import *
from functions.logging import write_log
from functions.plotting import plot_aligned

import pandas as pd

write_log(__file__)


output = 'save'
output_dir = current_dir / 'figures' / 'Supplement'
output_dir.mkdir(parents=True, exist_ok=True)
fname = 'Fig_SI_7_alignments_metric'

fit_kwargs = settings_dict['ecsecb_mtsecb']

output_folder = 'ecSecB_mtSecB'

states = ['ecSecB', 'mtSecB']
labels = ['$\it{ec}$SecB', '$\it{mt}$SecB']

# List of fit results: secondary structure
"""
List of input dataframes by fitting method
Batch fit (secondary structure alignment)
Batch fit (clustal)
Single fit 
Batch fit (alignment by residue number)

"""

fit_result_dict = {}

df = csv_to_protein(fitresults_dir / 'batch_fits' / output_folder / 'secondary_structure' / f'fit_result.csv')
df.rename(columns={'SecB WT apo': 'ecSecB', 'SecB_apo': 'mtSecB'}, inplace=True)
fit_result_dict['Secondary Structure'] = df

df = csv_to_protein(fitresults_dir / 'batch_fits' / output_folder / 'clustal' / f'fit_result.csv')
df.rename(columns={'SecB WT apo': 'ecSecB', 'SecB_apo': 'mtSecB'}, inplace=True)
fit_result_dict['Clustal Alignment'] = df

proteins = {state: csv_to_protein(fitresults_dir / 'ecSecB_mtSecB_single' / state / 'fit_result.csv').df for state in states}
df = pd.concat(proteins.values(), keys=proteins.keys(), axis=1)
fit_result_dict['Single fits'] = df

df = csv_to_protein(fitresults_dir / 'batch_fits' / output_folder / 'unaligned' / f'fit_result.csv')
df.rename(columns={'SecB WT apo': 'ecSecB', 'SecB_apo': 'mtSecB'}, inplace=True)
fit_result_dict['Unaligned'] = df

names = ['Single fits', 'Secondary Structure', 'Clustal Alignment', 'Unaligned']

fig, axes = pplt.subplots(ncols=1, nrows=4, aspect=7, axwidth=160/25.4, sharey=1)
ax_iter = iter(axes)

colors = ['pdark_blue', 'pdark_green']

format_kwargs = [
    {'color': color_dict['pdark_blue'], 'label': '$\it{ec}$SecB'},
    {'edgecolors': color_dict['pdark_green'], 'facecolors': 'none', 'label': '$\it{mt}$SecB'}
]

#indices = [0, 1, 3, 5]
#for idx, name in zip(indices, names):
for name in names:
    df = fit_result_dict[name]
    df_dict = {state: df[state] for state in states}

    aligned_dataframe = align_dataframes(df_dict, alignments['secondary_structure'])
    diffs = (aligned_dataframe[states[0]]['deltaG'] - aligned_dataframe[states[1]]['deltaG']).abs()
    title = f'{name}, Mean absolute difference: {diffs.mean()*1e-3:.2f} kJ/mol'
    ax = next(ax_iter)
    #ax = axes[idx]
    for kw, state in zip(format_kwargs, states):
        protein = aligned_dataframe[state]
        ax.axvspan(101, 110, color='#e8d261', alpha=0.5, lw=0, zorder=-10)
        ax.axvspan(173, 187, color='#e8d261', alpha=0.5, lw=0, zorder=-10)

        yvals = protein['deltaG']*1e-3
        ax.errorbar(protein.index, yvals, yerr=protein['covariance'] * 1e-3, **errorbar_kwargs, zorder=-1)
        ax.scatter(protein.index, yvals, **scatter_kwargs, **kw)

    ax.format(title=title)
axes.format(ylim=(45, 0), ylabel=dG_ylabel, xlabel='Alignment index')


hs = axes[0].get_legend_handles_labels()
axes[0].legend(*hs, loc='ul')


if output == 'show':
    plt.show()
elif output == 'save':
    plt.savefig(output_dir / f'{fname}.png')
    plt.savefig(output_dir / f'{fname}.pdf')
