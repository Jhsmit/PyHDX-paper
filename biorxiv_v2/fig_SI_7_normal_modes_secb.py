import matplotlib.pyplot as plt
import numpy as np
import proplot as pplt
from pyhdx.alignment import align_dataframes

from functions.align import *
from functions.base import *
from functions.fileIO import nma_result_to_df
from functions.formatting import *
from functions.logging import write_log


write_log(__file__)


output = 'save'
output_dir = current_dir / 'figures' / 'Supplement' / 'SecB_structures'
output_dir.mkdir(parents=True, exist_ok=True)
fname = 'Fig_SI_7_aligned_normal_modes'

width = 183 / 25.4 / 2

states = ['ecSecB', 'mtSecB']
color_dict = {'pdark_blue': [34, 34, 85], 'pdark_green': [34, 85, 34]}

fit_kwargs = settings_dict['ecsecb_mtsecb']

alignment_folder = 'secondary_structure'

input_folder = 'ecSecB_mtSecB'

df_nma = {name: nma_result_to_df(current_dir / 'nma' / name.lower()) for name in states}
aligned = align_dataframes(df_nma, alignments['secondary_structure'])
#
fix, axes = pplt.subplots(nrows=3, ncols=2, width=width*2, sharey=True, aspect=5)
for mode, ax in zip(range(7, 13), axes):
    col = f'displacement_mode_{mode}'

    for name, color in zip(states, color_dict.values()):
        vals = aligned[name][col].copy()
        line, = ax.plot(aligned.index, vals, label=name, c=np.array(color)/255)
        ax.format(title=f'Mode {mode}')
axes.format(ylabel='Displacement', xlabel='Aligned residue index')
axes[0].legend()

if output == 'show':
    plt.show()
elif output == 'save':
    plt.savefig(output_dir / f'{fname}.png')
    plt.savefig(output_dir / f'{fname}.pdf')

