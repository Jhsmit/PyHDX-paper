import matplotlib.pyplot as plt
from pyhdx.fileIO import csv_to_protein

from functions.align import alignments
from functions.base import *
from functions.formatting import *
from functions.logging import write_log
from functions.plotting import plot_aligned

write_log(__file__)


output = 'save'
output_dir = current_dir / 'figures' / 'Supplement'
output_dir.mkdir(parents=True, exist_ok=True)
fnames = ['Fig_SI_5_' + n for n in ['individual_fits', 'clustal', 'no_alignment']]

fit_kwargs = settings_dict['ecsecb_mtsecb']

output_folder = 'ecSecB_mtSecB'


states = ['ecSecB', 'mtSecB']
labels = ['$\it{ec}$SecB', '$\it{mt}$SecB']

alignments_names = list(alignments.keys())

# List of fit results: no alignment, single fits; clustal alignment; no alignment (batch fit)
fit_result_list = []

proteins = {state: csv_to_protein(current_dir / 'fits' / state / f"r1_{fit_kwargs['r1']}" / 'fit_result.csv').df for state in states}
fit_result_list.append(proteins)

for alignment_folder in alignments_names[1:]:
    df = csv_to_protein(current_dir / 'batch_fits' / output_folder / alignment_folder / f'fit_result.csv')
    df.rename(columns={'SecB WT apo': 'ecSecB', 'SecB_apo': 'mtSecB'}, inplace=True)
    fit_result_list.append(df)

    # These lines were missing, resulting in incorrect alignment in the biorxiv figure:
    # align_dataframes is a function from pyhdx

    # dfs = {state: df[state] for state in states}
    # aligned_df = align_dataframes(dfs, alignments[alignment_folder])
    #
    # fit_result_list.append(aligned_df)

for aligned_df, alignment, fname in zip(fit_result_list, alignments.values(), fnames):

    fig, axes = plt.subplots(2, figsize=(18.3/2.54, 3/2.54))
    plt.subplots_adjust(left=0.1, right=0.99, hspace=0.9)

    plot_aligned(axes, states, labels, aligned_df, alignment)

    if output == 'show':
        plt.show()
    elif output == 'save':
        plt.savefig(output_dir / f'{fname}.png')
        plt.savefig(output_dir / f'{fname}.pdf')
