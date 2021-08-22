#%%

from functions.formatting import *
from functions.base import *
#%%

aspect = 17.5
width = 50
output_dir = current_dir / 'figures'
output_dir.mkdir(parents=True, exist_ok=True)

current_dir = Path().resolve()
fig, ax = pplt.subplots(width=width/25.4, aspect=aspect)

ticks = np.linspace(10, 40, endpoint=True, num=4).astype(int)
tick_pos = rgb_norm(ticks)
cbar = ax.colorbar(rgb_cmap, loc='fill', ticks=tick_pos)
cbar.set_label(dG_ylabel, labelpad=0)
cbar.ax.set_xticklabels(ticks)
ax.format(xlabel=dG_ylabel)

plt.savefig(output_dir / 'Colorbar_abs_gibbs.png')
plt.savefig(output_dir / 'Colorbar_abs_gibbs.pdf')
#
fig, ax = pplt.subplots(width=width/25.4, aspect=aspect)

ticks = np.linspace(-10, 10, endpoint=True, num=5).astype(int)
tick_pos = diff_norm(ticks)
cbar = ax.colorbar(diff_cmap, loc='fill', ticks=tick_pos)
cbar.set_label(ddG_ylabel, labelpad=0)
cbar.ax.set_xticklabels(ticks)

plt.savefig(output_dir / 'Colorbar_diff_gibbs.png')
plt.savefig(output_dir / 'Colorbar_diff_gibbs.pdf')