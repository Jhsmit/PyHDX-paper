import matplotlib.pyplot as plt
import proplot as pplt
from scipy.stats import kde
import matplotlib as mpl
from matplotlib.axes import Axes

import numpy as np


def stripplot(data, ax=None, jitter=0.25, colors=None, offset=0., orientation='vertical', **scatter_kwargs):
    ax = ax or plt.gca()
    color_list = _prepare_colors(colors, len(data))

    for i, (d, color) in enumerate(zip(data, color_list)):
        jitter_offsets = (np.random.rand(d.size) - 0.5) * jitter
        cat_var = i * np.ones_like(d) + jitter_offsets + offset  # categorical axis variable
        if orientation == 'vertical':
            ax.scatter(cat_var, d, color=color, **scatter_kwargs)
        elif orientation == 'horizontal':
            ax.scatter(d, len(data) - cat_var, color=color, **scatter_kwargs)


def _prepare_colors(colors, N):
    if not isinstance(colors, list):
        return [colors]*N
    else:
        return colors


# From joyplot
def _x_range(data, extra=0.2):
    """ Compute the x_range, i.e., the values for which the
        density will be computed. It should be slightly larger than
        the max and min so that the plot actually reaches 0, and
        also has a bit of a tail on both sides.
    """
    try:
        sample_range = np.nanmax(data) - np.nanmin(data)
    except ValueError:
        return []
    if sample_range < 1e-6:
        return [np.nanmin(data), np.nanmax(data)]
    return np.linspace(np.nanmin(data) - extra*sample_range,
                       np.nanmax(data) + extra*sample_range, 1000)


def kdeplot(data, ax=None, offset=0., orientation='vertical',
            linecolor=None, linewidth=None, zero_line=True, x_extend=1e-3, y_scale=None, y_norm=None, fillcolor=False, fill_cmap=None,
            fill_norm=None):
    assert not (y_scale and y_norm), "Cannot set both 'y_scale' and 'y_norm'"
    y_scale = 1. if y_scale is None else y_scale

    color_list = _prepare_colors(linecolor, len(data))

    for i, (d, color) in enumerate(zip(data, color_list)):
        #todo remove NaNs?

        # Perhaps also borrow this part from joyplot
        kde_func = kde.gaussian_kde(d)
        kde_x = _x_range(d, extra=0.4)
        kde_y = kde_func(kde_x)*y_scale
        if y_norm:
            kde_y = y_norm*kde_y / kde_y.max()
        bools = kde_y > x_extend * kde_y.max()
        kde_x = kde_x[bools]
        kde_y = kde_y[bools]

        cat_var = len(data) - i + kde_y + offset # x in horizontal
        cat_var_zero = (len(data) - i)*np.ones_like(kde_y) + offset

        # x = i * np.ones_like(d) + jitter_offsets + offset  # 'x' like, could be y axis
        if orientation == 'horizontal':
            plot_x = kde_x
            plot_y = cat_var
            img_data = kde_x.reshape(1, -1)
        elif orientation == 'vertical':
            plot_x = len(data) - cat_var
            plot_y = kde_x
            img_data = kde_x[::-1].reshape(-1, 1)
        else:
            raise ValueError(f"Invalid value '{orientation}' for 'orientation'")

        line, = ax.plot(plot_x, plot_y, color=color, linewidth=linewidth)
        if zero_line:
            ax.plot([plot_x[0], plot_x[-1]], [plot_y[0], plot_y[-1]], color=line.get_color(), linewidth=linewidth)

        if fillcolor:
            #todo refactor to one if/else orientation
            color = line.get_color() if fillcolor is True else fillcolor
            if orientation == 'horizontal':
                ax.fill_between(kde_x, plot_y, np.linspace(plot_y[0], plot_y[-1], num=plot_y.size, endpoint=True),
                                color=color)
            elif orientation == 'vertical':
                ax.fill_betweenx(kde_x, len(data) - cat_var, len(data) - cat_var_zero, color=color)

        if fill_cmap:
            fill_norm = fill_norm or (lambda x: x)
            color_img = fill_norm(img_data)

            xmin, xmax = np.min(plot_x), np.max(plot_x)
            ymin, ymax = np.min(plot_y), np.max(plot_y)
            extent = [xmin-offset, xmax-offset, ymin, ymax] if orientation == 'horizontal' else [xmin, xmax, ymin-offset, ymax-offset]
            im = Axes.imshow(ax, color_img, aspect='auto', cmap=fill_cmap, extent=extent)  # left, right, bottom, top
            fill_line, = ax.fill(plot_x, plot_y, facecolor='none')
            im.set_clip_path(fill_line)


def boxplot(data, ax, offset=0., orientation='vertical', widths=0.25, linewidth=None, linecolor=None, **kwargs):
    if orientation == 'vertical':
        vert = True
        positions = np.arange(len(data)) + offset
    elif orientation == 'horizontal':
        vert = False
        positions = len(data) - np.arange(len(data)) - offset
    else:
        raise ValueError(f"Invalid value '{orientation}' for 'orientation', options are 'horizontal' or 'vertical'")

    #todo for loop
    boxprops = kwargs.pop('boxprops', {})
    whiskerprops = kwargs.pop('whiskerprops', {})
    medianprops = kwargs.pop('whiskerprops', {})

    boxprops['linewidth'] = linewidth
    whiskerprops['linewidth'] = linewidth
    medianprops['linewidth'] = linewidth

    boxprops['color'] = linecolor
    whiskerprops['color'] = linecolor
    medianprops['color'] = linecolor

    Axes.boxplot(ax, data, vert=vert, positions=positions, widths=widths, boxprops=boxprops, whiskerprops=whiskerprops,
               medianprops=medianprops, **kwargs)


def label_axes(labels, ax, offset=0., orientation='vertical', **kwargs):
    #todo check offset sign
    if orientation == 'vertical':
        ax.set_xticks(np.arange(len(labels)) + offset)
        ax.set_xticklabels(labels, **kwargs)
    elif orientation == 'horizontal':
        ax.set_yticks(len(labels) - np.arange(len(labels)) + offset)
        ax.set_yticklabels(labels, **kwargs)


if __name__ == '__main__':
    np.random.seed(43)

    num_samples = 5
    data = [np.random.normal(np.random.randint(3, 9), scale=np.random.rand() + 0.5, size=np.random.randint(20, 500)) for
            f in range(num_samples)]

    node_pos = [1, 5, 9]
    norm = plt.Normalize(node_pos[0], node_pos[-1], clip=True)
    nodes = norm(node_pos)
    linear_colors = np.array(['#ff0000', '#00ff00', '#0000ff'])
    custom_cmap = mpl.colors.LinearSegmentedColormap.from_list("custom_cmap", list(zip(nodes, linear_colors)))

    labels = 'ABCDE'

    fig, ax = pplt.subplots()
    orientation = 'horizontal'
    stripplot(data, ax=ax, offset=0.0, orientation=orientation, s=5)
    kdeplot(data, ax=ax[0], linecolor='k', offset=0.2, orientation=orientation, fillcolor=False, fill_cmap=custom_cmap,
            fill_norm=norm, y_scale=None, y_norm=0.3)
    boxplot(data, ax=ax, offset=0.3, sym='', linewidth=1.5, linecolor='k', orientation=orientation)
    label_axes(labels, ax=ax, orientation=orientation)

    ax.format(xlim=(-2, 11))

    cbar = fig.colorbar(custom_cmap, width=0.1, ticks=[0, 0.5, 1])
    cbar.ax.set_yticklabels([1, 2, 3])
    ylim = ax.get_ylim()
    ax.set_ylim(None, ylim[1] + 0.25)
    plt.show()