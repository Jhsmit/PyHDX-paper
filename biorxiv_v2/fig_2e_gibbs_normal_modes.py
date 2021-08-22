import pandas as pd
import scipy
from pyhdx.fileIO import csv_to_protein

from functions.base import *
from functions.fileIO import nma_result_to_df
from functions.formatting import *
from functions.logging import write_log

write_log(__file__)


output = 'save'
output_dir = current_dir / 'figures' / 'Fig_2'
output_dir.mkdir(parents=True, exist_ok=True)
fname = 'Fig_2e_normal_modes_gibbs'

width = 183 / 25.4 / 2

fit_kwargs = settings_dict['ecsecb_mtsecb']

output_folder = 'ecSecB_mtSecB'
alignment_folder = 'secondary_structure'

df = csv_to_protein(current_dir / 'batch_fits' / output_folder / alignment_folder / f'fit_result.csv')
df.rename(columns={'SecB WT apo': 'ecSecB', 'SecB_apo': 'mtSecB'}, inplace=True)
names = ['ecSecB', 'mtSecB']

df_nma = {name: nma_result_to_df(current_dir / 'nma' / name.lower()) for name in names}
combined_dfs = {name: pd.concat([df_nma[name], df[name]], axis=1) for name in names}

fig, axes = pplt.subplots(nrows=1, ncols=2, aspect=4, width=2*width, sharey=False, sharex=False)
for name, ax in zip(names, axes):
    df = combined_dfs[name]

    na_removed = df.dropna(how='any')
    x = na_removed['displacement']
    y = na_removed['deltaG']
    a, b = np.polyfit(x, y, 1)
    rho, p = scipy.stats.pearsonr(x, y)
    print(name, rho, p)

    ax.plot(df.index, df['displacement'], color='k', label='Displacement')

    yvals = df['deltaG'] * 1e-3
    rgba_colors = rgb_cmap(rgb_norm(yvals), bytes=True)
    hex_colors = rgb_to_hex(rgba_colors)

    ax1 = ax.twinx()
    ax1.errorbar(df.index, yvals, yerr=df['covariance'] * 1e-3, **errorbar_kwargs, zorder=-1)
    ax1.scatter(df.index, yvals, c=hex_colors, **scatter_kwargs)
    ax1.format(ylabel='', ylim=(40, 0))

    if name == 'mtSecB':
        tick_labels = [0, 20, 40]
        add_colorbar(ax1, ax1, rgb_cmap, rgb_norm, tick_labels=tick_labels, label=dG_ylabel)

    if name == 'ecSecB':
        ax1.set_ylabel(dG_ylabel)
        prefix = '$\it{ec}$'
    else:
        prefix = '$\it{mt}$'

    title = prefix + f"SecB: r={rho:.2f}"
    ax.format(title=title, ylabel='Displacement', xlabel=r_xlabel)

axes[0].format(ylim=(0, 7))
axes[1].format(ylim=(0, 7))

if output == 'show':
    plt.show()
elif output == 'save':
    plt.savefig(output_dir / f'{fname}.png')
    plt.savefig(output_dir / f'{fname}.pdf')

