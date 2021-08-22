import scipy
from pyhdx.fileIO import csv_to_protein

from functions.formatting import *
from functions.logging import write_log

write_log(__file__)


output = 'save'
output_dir = current_dir / 'figures' / 'supplement'
output_dir.mkdir(parents=True, exist_ok=True)
fname = 'Fig_SI_8_flexibility_NMA'

nma_data_dir = current_dir / 'nma'

r1 = 0.5
fit_dir = 'fits'

label_lut = {
    'TF Dimer': ('tf_dimer', 'TF_dimer_4C'),
    'TF': ('tf_monomer_5owi', 'TF_monomer_4C'),
    'SecA': ('seca_monomer', 'SecA_monomer'),
    'SecA Dimer': ('seca_dimer_1nl3', 'SecA_1-901_wt_apo'),
    'MBP': ('mbp', 'MBP_wt'),
    'PPiA': ('ppia', 'PPiA_WT'),
    'PPiB': ('ppib', 'PPiB_WT'),
    'ecSecB': ('ecsecb', 'ecSecB'),
    'mtSecB': ('mtsecb', 'mtSecB'),
    'hPREP': ('hprep', 'hPREP_default'),
    'Bcl-2': ('bcl2_1gjh', 'Bcl2'),
    'EscV': ('escv', 'EscV')
}

flex_dict = {}
for name, v in label_lut.items():
    nma_label, hdx_label = v
    f_pth = nma_data_dir / nma_label / 'modes.txt'
    eigenvalues = np.array(f_pth.read_text().split('\n')[4].split(' ')).astype(float)[6:]
    N = len(f_pth.read_text().split('\n')[5].split(' '))
    flexibility = np.sqrt(np.sum(1/eigenvalues)/N)
    print(name, N, flexibility, eigenvalues.min(), eigenvalues.max())
    flex_dict[name] = flexibility

#%%

gibbs_dfs = {}
for k, v in label_lut.items():
    print(k)
    nma_label, hdx_label = v
    gibbs_pth = current_dir / fit_dir / hdx_label / f'r1_{r1}' / 'fit_result.csv'

    protein = csv_to_protein(gibbs_pth)
    gibbs_dfs[k] = protein

#%%
excluded = []
selected_gibbs = {k:v for k, v in gibbs_dfs.items() if k not in excluded}
selected_flex = {k:v for k, v in flex_dict.items() if k not in excluded}

medians = [np.nanmedian(d['deltaG'])*1e-3 for d in selected_gibbs.values()]
means = [np.mean(d['deltaG'])*1e-3 for d in selected_gibbs.values()]
nma_values = list(selected_flex.values())

x = medians
y = nma_values

rho, p = scipy.stats.pearsonr(x, y)

print(rho, p)

fig, ax = pplt.subplots(width=80/25.4, aspect=1.2)
ax.scatter(x, y, color='k')

a, b = np.polyfit(x, y, 1)
xvec = np.array([5, 40])
yvec = a*xvec + b
offsets = {
    'TF Dimer': (-2, 0.03),
    'Bcl-2': (-5.8, -0.02),
    'ecSecB': (-7.5, -0.0),
    'mtSecB': (0, -0.025),
    'PPiB': (-3.2, -0.05),
    'PPiA': (0, 0.02),
    'MBP': (0.5, -0.01)
}

xlim = ax.get_xlim()
ylim = ax.get_ylim()
ax.plot(xvec, yvec, color='r', zorder=-5)
ax.set_xlim(5, 40)
ax.set_ylim(0.1, 0.8)

for i, txt in enumerate(selected_flex.keys()):
    x_offset, y_offset = offsets.get(txt, (0, 0))
    ax.annotate(txt, (x[i]+0.5+x_offset, y[i]+y_offset), fontsize=10, rotation=0, horizontalalignment='left')

ax.set_ylabel('Normal mode flexibility')
ax.set_xlabel('Mean ΔG (kJ/mol)')
ax.set_title(f'Protein flexibility, pearson r={rho:.2f} (p={p:.1e})')

if output == 'show':
    plt.show()
elif output == 'save':
    plt.savefig(output_dir / f'{fname}.png')
    plt.savefig(output_dir / f'{fname}.pdf')
