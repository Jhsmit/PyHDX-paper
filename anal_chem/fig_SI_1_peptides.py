
from pyhdx.batch_processing import yaml_to_hdxm
from pyhdx.plot import plot_peptides

from functions.base import *
from functions.formatting import *
from functions.logging import write_log

write_log(__file__)

output = 'save'
output_dir = current_dir / 'figures' / 'Supplement'
output_dir.mkdir(parents=True, exist_ok=True)
fname_time = 'Fig_SI_1a_Peptides_SecB_time'
fname_coverage = 'Fig_SI_1b_Peptides_SecB_coverage'

figure_width = page_width/25.4
hdxm = yaml_to_hdxm(data_dict['ecSecB'], input_data_dir)

values = [0, 10, 30, 1, 5, 10, 100]
units = 3*['s'] + 4*['min']
labels = [f't = {v} {u}' for v, u in zip(values, units)]

fig, axes = pplt.subplots(ncols=2, nrows=4, sharex=True, width=figure_width, aspect=4)
axes_list = list(axes[:, 0]) + list(axes[:, 1])

for label, ax, pm in zip(labels, axes_list, hdxm):
    plot_peptides(pm, ax, linewidth=0.5)
    ax.format(xlim=(0, 160), title=label)

axes.format(xlabel=r_xlabel)
axes_list[-1].set_axis_off()
#%%
cmap = mpl.cm.get_cmap('jet')
cbar = fig.colorbar(cmap, width=cbar_width, label='RFU', tickminor=True)
cbar.ax.set_yticklabels(np.array([0, 20, 40, 60, 80, 100])/100.)

if output == 'show':
    plt.show()
elif output == 'save':
    plt.savefig(output_dir / f'{fname_coverage}.png')
    plt.savefig(output_dir / f'{fname_coverage}.pdf')


#
# Np = hdxm.coverage.X.shape[0]
# Np // 16
#
# gen = count(0, 3)
# idx = list(islice(count(0, 5), 12))
# #%%
# fig, axes = pplt.subplots(ncols=4, nrows=3, width=figure_width, aspect=2)
#
# for peptide_index, ax in zip(idx, axes):
#     ax.set_xscale('log')
#     ax.plot(hdxm.timepoints, [p.data['uptake'][peptide_index] for p in hdxm],
#                  label=f'Peptide {peptide_index}',
#                  marker='o', markersize=4, color='k')
#     ax.format(title=f'Peptide {peptide_index}')
#
# axes.format(xlabel='Time (min)', ylabel='Peptide Deuteration (Da)')
#
# if output == 'show':
#     plt.show()
# elif output == 'save':
#     plt.savefig(output_dir / f'{fname_time}.png')
#     plt.savefig(output_dir / f'{fname_time}.pdf')
#
