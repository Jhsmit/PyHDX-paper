from Bio.SubsMat.MatrixInfo import pam250
from matplotlib.patches import Rectangle, Patch
from pyhdx.fileIO import csv_to_dataframe, load_fitresult
from pyhdx.plot import plot_peptides

from anal_chem.functions.base import settings_dict, fitresults_dir
from anal_chem.functions.formatting import *


try:
    from pyhdx.output import Output, Report
except ModuleNotFoundError:
    print('no pylatex ')

def alignment_score(s1, s2):
    """Calculate PAM250 similarty score between two FASTA strings and return the resutl as string indicators"""
    alignment = ''
    for s in zip(s1, s2):

        if '-' in s:
            alignment += ' '
        elif s[0] == s[1]:
            alignment += '*'
        else:
            try:
                score = pam250[s]
            except KeyError:
                score = pam250[s[::-1]]
            if score == 0.:
                alignment += ' '
            elif score <= 0.5:
                alignment += '.'
            elif score > 0.5:
                alignment += ':'
            else:
                raise ValueError('oei')
    return alignment


def plot_aligned(axes, names, labels, aligned_dataframe, alignments_dict, size=95):
    aa_font_size = 5

    #size = size  # dont chagne this unless manually changeing y labels
    start = -0.5
    ranges = [(start, start + size), (start + size, start + 2 * size)]
    for ax, r in zip(axes, ranges):
        for i, (name, _label) in enumerate(zip(names, labels)):

            vals = aligned_dataframe[name]['deltaG']*1e-3
            #sequence = aligned_dataframe[name]['sequence']

            colors_rgba = rgb_cmap(rgb_norm(vals), bytes=True)
            color_arr = rgb_to_hex(colors_rgba)
            #color_arr[np.isnan(vals)] = '#8c8c8c'

            for j, (c, s) in enumerate(zip(color_arr, alignments_dict[name])):

                #Place colored rectangle
                color = '#ffffff' if s == '-' else c
                rect = Rectangle((j - 0.5, 2 - 2 * i - 0.5), width=1, height=1, color=color)
                ax.add_patch(rect)
                text_color = 'black'  # if np.mean(hex_to_rgb(color)) < 200 else 'white'

                # Place residue one letter codes
                if j > r[0] and j < r[1]:
                    ax.annotate(s.upper(), (j, 2 - i * 2), ha='center', va='center', color=text_color, size=aa_font_size,
                                font=font_pth)

            # Place labels
            ax.text(r[0] - 10, 2 - 2 * i, _label, size=7, verticalalignment='center')

        # add alignment text
        score = alignment_score(*[alignments_dict[name].upper() for name in names])


        for j, s in enumerate(score):
            if j > r[0] and j < r[1]:
                ax.annotate(s.upper(), (j, 1), ha='center', va='center', color=text_color, size=aa_font_size,
                            font=font_pth)
                #ax.text(j, 1, s, horizontalalignment='center', verticalalignment='center', color='k', size=aa_font_size)

        ax.set_xlim(r)
        ax.set_ylim(-0.5, 2.9)
        ax.set_yticks([])
        ax.set_xticks([])
        ax.tick_params(axis=u'both', which=u'both', length=0)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)

    top_part = alignments_dict['ecSecB'][:size]
    top_num = sum([c.isalpha() for c in top_part])

    bot_part = alignments_dict['mtSecB'][:size]
    bot_num = sum([c.isalpha() for c in bot_part])

    n = 0
    axes[n].text(ranges[n][0] - 3, 2, '1', size=7, verticalalignment='center')
    axes[n].text(ranges[n][0] - 3, 0, '1', size=7, verticalalignment='center')

    n = 1
    axes[n].text(ranges[n][0] - 3, 2, top_num + 1, size=7, verticalalignment='center')
    axes[n].text(ranges[n][0] - 3, 0, bot_num + 1, size=7, verticalalignment='center')


def plot_aligned_nocolor(axes, names, labels, alignments_dict, size=95):
    aa_font_size = 5

    #size = size  # dont chagne this unless manually changeing y labels
    start = -0.5
    ranges = [(start, start + size), (start + size, start + 2 * size)]
    for ax, r in zip(axes, ranges):
        # for i, (name, _label) in enumerate(zip(names, labels)):
        #
        #     vals = aligned_dataframe[name]['deltaG']*1e-3
        #     #sequence = aligned_dataframe[name]['sequence']
        #
        #     colors_rgba = rgb_cmap(rgb_norm(vals), bytes=True)
        #     color_arr = rgb_to_hex(colors_rgba)
        #     #color_arr[np.isnan(vals)] = '#8c8c8c'
        #
        #     for j, (c, s) in enumerate(zip(color_arr, alignments_dict[name])):
        #
        #         #Place colored rectangle
        #         color = '#ffffff' if s == '-' else c
        #         rect = Rectangle((j - 0.5, 2 - 2 * i - 0.5), width=1, height=1, color=color)
        #         ax.add_patch(rect)
        #         text_color = 'black'  # if np.mean(hex_to_rgb(color)) < 200 else 'white'
        #
        #         # Place residue one letter codes
        #         if j > r[0] and j < r[1]:
        #             ax.annotate(s.upper(), (j, 2 - i * 2), ha='center', va='center', color=text_color, size=aa_font_size,
        #                         font=font_pth)
        #
        #     # Place labels
        #     ax.text(r[0] - 10, 2 - 2 * i, _label, size=7, verticalalignment='center')

        # add alignment text
        score = alignment_score(*[alignments_dict[name].upper() for name in names])

        for j, s in enumerate(score):
            if j > r[0] and j < r[1]:
                text_color = 'black'
                ax.annotate(s.upper(), (j, 1), ha='center', va='center', color=text_color, size=aa_font_size,
                            font=font_pth)
                for i, name in enumerate(names):
                    aa = alignments_dict[name].upper()[j]
                    ax.annotate(aa.upper(), (j, 2 - i * 2), ha='center', va='center', color=text_color, size=aa_font_size,
                                    font=font_pth)
                if j % 10 == 0:
                    ax.annotate(str(j), (j, 3), ha='center', va='center', color=text_color, size=aa_font_size,
                                font=font_pth)

        # Place labels
        for i, _label in enumerate(labels):
            ax.text(r[0] - 8, 2 - 2 * i, _label, size=7, verticalalignment='center')
                #ax.text(j, 1, s, horizontalalignment='center', verticalalignment='center', color='k', size=aa_font_size)

        ax.set_xlim(r)
        ax.set_ylim(-0.5, 3.9)
        ax.set_yticks([])
        ax.set_xticks([])
        ax.tick_params(axis=u'both', which=u'both', length=0)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)

    top_part = alignments_dict['ecSecB'][:size]
    top_num = sum([c.isalpha() for c in top_part])

    bot_part = alignments_dict['mtSecB'][:size]
    bot_num = sum([c.isalpha() for c in bot_part])

    # n = 0
    # axes[n].text(ranges[n][0] - 3, 2, '1', size=7, verticalalignment='center')
    # axes[n].text(ranges[n][0] - 3, 0, '1', size=7, verticalalignment='center')
    #
    # n = 1
    # axes[n].text(ranges[n][0] - 3, 2, top_num + 1, size=7, verticalalignment='center')
    # axes[n].text(ranges[n][0] - 3, 0, bot_num + 1, size=7, verticalalignment='center')
