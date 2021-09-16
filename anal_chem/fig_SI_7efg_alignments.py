import matplotlib.pyplot as plt
from pyhdx.fileIO import csv_to_protein
from pyhdx.alignment import align_dataframes

from functions.align import alignments
from functions.base import *
from functions.formatting import *
from functions.logging import write_log
from functions.plotting import plot_aligned_nocolor

write_log(__file__)


output = 'save'
output_dir = current_dir / 'figures' / 'Supplement'
output_dir.mkdir(parents=True, exist_ok=True)
names = ['secondary_structure', 'clustal', 'unaligned']
fnames = ['Fig_SI_7_' + n for n in names]

titles = {
    'secondary_structure': 'Secondary structure alignment',
    'clustal': 'Clustal alignment',
    'unaligned': 'No alignment (aligned by residue number)'
}

fit_kwargs = settings_dict['ecsecb_mtsecb']

output_folder = 'ecSecB_mtSecB'

states = ['ecSecB', 'mtSecB']
labels = ['$\it{ec}$SecB', '$\it{mt}$SecB']

sizes = iter([95, 95, 95])
alignment_names = iter(names)
for fname in fnames:
    fig, axes = plt.subplots(2, figsize=(page_width/25.4, 3/2.54))
    plt.subplots_adjust(left=0.1, right=0.99, hspace=0.2)

    name = next(alignment_names)
    alignment = alignments[name]
    plot_aligned_nocolor(axes, states, labels, alignment, size=next(sizes))
    plt.suptitle(titles[name])

    if output == 'show':
        plt.show()
    elif output == 'save':
        plt.savefig(output_dir / f'{fname}.png')
        plt.savefig(output_dir / f'{fname}.pdf')

    print(fname)

