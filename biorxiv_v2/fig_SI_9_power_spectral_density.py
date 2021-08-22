import pandas as pd
from pycorrelate import ucorrelate
from pyhdx.fileIO import csv_to_protein

from functions.formatting import *
from functions.logging import write_log
import pyhdx
print(pyhdx.__version__)
write_log(__file__)


output = 'save'
output_dir = current_dir / 'figures' / 'Supplement'
output_dir.mkdir(parents=True, exist_ok=True)
fname = 'Fig_SI_9_PSD'


r1 = 0.05
N = 300


states = [
          'MBP_wt',
          'PPiA_WT',
          'PPiB_WT',
          'random'
    ]

labels = [
    'MBP',
    'PpiA',
    'PpiB',
    'Random Uniform'
]

fit_dir = 'fits'

proteins = {}
for state in states:
    if state == 'random':
        y = np.random.random_sample(N)
        protein = pd.DataFrame({'_deltaG': y})

    else:
        f = current_dir / fit_dir / state / f'r1_{r1}' / 'fit_result.csv'
        protein = csv_to_protein(f)

    proteins[state] = protein

keys = proteins.keys()
for i, key in enumerate(keys):
    protein = proteins[key]
    y = protein['_deltaG'].to_numpy()
    zero = y - np.mean(y)
    norm = zero / zero.std()
    corr = ucorrelate(norm, norm)
    corr = corr / corr[0]

    ft = np.fft.fft(corr)
    freq = np.fft.fftfreq(corr.shape[-1])
    fig, axes = pplt.subplots(nrows=2, aspect=3, width=183/25.4/2, sharex=False, sharey=False)

    # We look for peak in the 10-40 region
    masked_corr = corr.copy()
    masked_corr[:10] = 0
    masked_corr[40:] = 0
    idx = np.argmax(masked_corr)
    ax1 = axes[0]
    ax1.plot(corr)

    ax1.axhline(0, color='k')

    if key != 'random':
        ax1.scatter(idx, np.max(masked_corr), color='r')
        title = f'Autocorrelation: Period {idx}'
    else:
        title = f'Autocorrelation'

    ax1.format(ylabel='ACF Amplitude (norm)', xlabel='Lag (Residues)', title=title)

    max = np.nan
    ax2 = axes[1]
    ax2.plot(1 / freq, ft.real)
    idx = np.argmax(ft.real)

    pos_freq = freq[:freq.size // 2]
    cutoff = np.argmin(np.abs(pos_freq - 1 / 100))
    selected_freq = pos_freq[cutoff:]
    selected_amp = ft.real[:freq.size // 2][cutoff:]
    idx = np.argmax(selected_amp)
    max_period = 1 / selected_freq[idx]

    max_freq = freq[idx]
    print(key, max_period)
    ax2.scatter(max_period, np.max(selected_amp), color='r')
    ax2.format(xlim=(0, 100), xlabel='1/Frequency (Residues)',
               title=f'Power Spectral Density: Period {max_period:.2f}', ylabel='PSD (a.u.)')
    axes.format(suptitle=labels[i])

    if output == 'show':
        plt.show()
    elif output == 'save':
        plt.savefig(output_dir / f'{fname}_{labels[i]}.png')
        plt.savefig(output_dir / f'{fname}_{labels[i]}.pdf')
