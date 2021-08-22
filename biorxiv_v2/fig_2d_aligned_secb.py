from functions.base import *
from functions.formatting import *
from functions.plotting import plot_aligned
from functions.align import alignments
import matplotlib.pyplot as plt
from pyhdx.fileIO import csv_to_protein
from pyhdx.alignment import align_dataframes

output = 'save'
output_dir = current_dir / 'figures' / 'Fig_2'
output_dir.mkdir(parents=True, exist_ok=True)
fname = 'Fig_2d_ecSecB_mtSecB_aligned'

fit_kwargs = settings_dict['ecsecb_mtsecb']

output_folder = 'ecSecB_mtSecB'
alignment_folder = 'secondary_structure'

df = csv_to_protein(current_dir / 'batch_fits' / output_folder / alignment_folder / f'fit_result.csv')
df.rename(columns={'SecB WT apo': 'ecSecB', 'SecB_apo': 'mtSecB'}, inplace=True)

states = ['ecSecB', 'mtSecB']
labels = ['$\it{ec}$SecB', '$\it{mt}$SecB']

dfs = {state: df[state] for state in states}
aligned_df = align_dataframes(dfs, alignments[alignment_folder])

fig, axes = plt.subplots(2, figsize=(18.3/2.54, 3/2.54))
plt.subplots_adjust(left=0.1, right=0.99, hspace=0.9)

plot_aligned(axes, states, labels, aligned_df, alignments[alignment_folder])

if output == 'show':
    plt.show()
elif output == 'save':
    plt.savefig(output_dir / f'{fname}.png')
    plt.savefig(output_dir / f'{fname}.pdf')
