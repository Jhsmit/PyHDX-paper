import proplot as pplt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import to_hex
from pyhdx.support import rgb_to_hex

from anal_chem.functions.base import current_dir

font_pth = current_dir / 'fonts' / 'COURIER.TTF'
page_width = 177.8 # full anal chem two column width in mm

settings = {
    'font.family': 'Arial',
    'font.size': 7}

pplt.rc.update(settings)

node_pos = [10, 25, 40]  # in kJ/mol
linear_colors = ['#ff0000', '#00ff00', '#0000ff']  # red, green, blue
#linear_colors = ['#ee3377', '#009988', '#33bbee']  # magenta, teal, cyan
no_coverage = '#8c8c8c'

color_dict = {'pdark_blue': [34, 34, 85], 'pdark_green': [34, 85, 34]}
color_dict = {k: to_hex(np.array(v)/255.) for k, v in color_dict.items()}

rgb_norm = plt.Normalize(node_pos[0], node_pos[-1], clip=True)
rgb_cmap = mpl.colors.LinearSegmentedColormap.from_list("rgb_cmap", list(zip(rgb_norm(node_pos), linear_colors)))
rgb_cmap.set_bad(color=no_coverage)

diff_colors = ['#54278e', '#ffffff', '#006d2c'][::-1]
diff_node_pos = [-10, 0, 10]
diff_norm = plt.Normalize(diff_node_pos[0], diff_node_pos[-1], clip=True)
diff_cmap = mpl.colors.LinearSegmentedColormap.from_list("diff_cmap", list(zip(diff_norm(diff_node_pos), diff_colors)))
diff_cmap.set_bad(color=no_coverage)

cbar_width = 0.075

dG_ylabel = 'ΔG (kJ/mol)'
ddG_ylabel = 'ΔΔG (kJ/mol)'

r_xlabel = 'Residue Number'


errorbar_kwargs = {
    'fmt': 'o',
    'ecolor': 'k',
    'elinewidth': 0.2,
    'markersize': 0,
    'alpha': 0.5
}

scatter_kwargs = {
    's': 7
}


def add_colorbar(fig, ax, cmap, norm, tick_labels, label=None, num=100):
    ymin, ymax = ax.get_ylim()
    values = np.linspace(ymin, ymax, endpoint=True, num=num)
    colors = cmap(norm(values))
    if ymin > ymax:  # Reversed y axis
        colors = colors[::-1]
    cbar = fig.colorbar(colors, values=values, ticks=tick_labels, space=0, width=cbar_width, label=label)
    ax.format(yticklabelloc='None', ytickloc='None')

    return cbar
