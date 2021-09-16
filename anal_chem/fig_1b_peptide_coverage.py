import matplotlib as mpl
import matplotlib.pyplot as plt
import proplot as pplt
from pyhdx.batch_processing import yaml_to_hdxm
from pyhdx.plot import plot_peptides

from functions.base import *
from functions.formatting import *
from functions.logging import write_log

write_log(__file__)

output = 'save'
output_dir = current_dir / 'figures' / 'Fig_1'
output_dir.mkdir(parents=True, exist_ok=True)
fname = 'Fig_1b_peptide_coverage'


figure_width = 60 # Width in mm
figure_width /= 25.4

figure_height = 30
figure_height /= 25.4


hdxm = yaml_to_hdxm(data_dict['ecSecB'], data_dir=input_data_dir)

fig, ax = pplt.subplots(ncols=1, nrows=1, height=figure_height, sharex=False, sharey=False, width=figure_width,
                          wspace=0.1)


plot_peptides(hdxm[2], ax, intervals='corrected')
ax.format(xlabel=r_xlabel, xlim=(0, 160))

cmap = mpl.cm.get_cmap('jet')
cbar = fig.colorbar(cmap, width=cbar_width, ticks=[0, 1], space=0)
cbar.ax.set_yticklabels([0, 1])
cbar.set_label('RFU', labelpad=-0)


if output == 'show':
    plt.show()
elif output == 'save':
    plt.savefig(output_dir / f'{fname}.pdf')
    plt.savefig(output_dir / f'{fname}.png')

