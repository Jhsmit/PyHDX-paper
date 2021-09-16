import numpy as np
from pyhdx.fileIO import csv_to_protein
from pyhdx.support import color_pymol, apply_cmap
from pymol import cmd

from functions.base import *
from functions.formatting import *
from functions.pymol import *
from functions.logging import write_log


write_log(__file__)


output_dir = current_dir / 'figures' / 'Fig_2'
output_dir.mkdir(parents=True, exist_ok=True)
fname_diff = 'Fig_2c_seca_seca_adp_diff'

#%%

name = 'SecA_monomer'
views = structure_views[name]['views']
pdb_file = structure_views[name]['filename']
pdb_path = current_dir / 'structures' / pdb_file

fit_kwargs = settings_dict['seca_all_states']

names = ['SecA_monomer', 'SecA_monomer_ADP']
states = ['SecA-monomer', 'SecA mono ADP']

sclf = 1
px, py = 1000, 1000
px *= sclf
py *= sclf

protein = csv_to_protein(fitresults_dir / 'batch_fits' / 'SecA_states' / f'fit_result.csv')

cmd.reinitialize()
cmd.load(str(pdb_path))
cmd.set_view(views[0])
cmd.set('antialias', 2)
cmd.set('fog', 0)
cmd.color('magenta')  # Color full protein magenta to show any potential residues missed by ΔG coloring

# Individual states absolute ΔG render
for name, state in zip(names, states):
    # Reindex such that pandas series includes uncovered N-terminal part (padded with NaN's)
    dG = protein[state]['deltaG'].reindex(np.arange(1, protein.index.max() + 1)) * 1e-3
    # Convert deltaG values to colors
    colors = apply_cmap(dG, rgb_cmap, rgb_norm)
    color_pymol(colors, cmd)

    cmd.ray(px, py, renderer=0, antialias=2)
    output_file = output_dir / f'Fig_2b_{name}_pymol_render.png'
    cmd.png(output_file)


# Differences render
ddG = protein[states[1]]['deltaG'] - protein[states[0]]['deltaG']
ddG = ddG.reindex(np.arange(1, protein.index.max() + 1)) * 1e-3

cmd.color('magenta')
colors = apply_cmap(ddG, diff_cmap, diff_norm)
color_pymol(colors, cmd)

cmd.ray(px, py, renderer=0, antialias=2)
output_file = output_dir / f'{fname_diff}.png'
cmd.png(output_file)



