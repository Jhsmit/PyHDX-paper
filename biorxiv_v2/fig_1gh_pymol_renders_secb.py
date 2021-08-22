from pyhdx.fileIO import csv_to_protein
from pyhdx.support import color_pymol, apply_cmap
from pymol import cmd

from functions.base import *
from functions.formatting import *
from functions.logging import write_log
from functions.pymol import structure_views

write_log(__file__)

output = 'save'
output_dir = current_dir / 'figures' / 'Fig_1'
output_dir.mkdir(parents=True, exist_ok=True)
fname_tetramer = 'Fig_1g_ecsecb_tetramer'
fname_tetramer_dimer = 'Fig_1h_ecsecb_tetramer_dimer'

name = 'ecSecB'
views = structure_views[name]['views']
pdb_file = structure_views[name]['filename']
pdb_path = current_dir / 'structures' / pdb_file

fit_kwargs = settings_dict['ecsecb_tetramer_dimer']

sclf = 2
px, py = 640, 480
px *= sclf
py *= sclf


def do_render(colors, views, fname):
    for i, view in enumerate(views):
        cmd.set_view(view)
        cmd.set('antialias', 2)
        cmd.set('fog', 0)

        color_pymol(colors, cmd)
        cmd.ray(px, py, renderer=0, antialias=2)
        output_file = output_dir / f'{fname}_pymol_render_v{i}.png'
        cmd.png(output_file)


#-------------------------------------------------------------------------#
# ecSecB tetramer Single
#-------------------------------------------------------------------------#

state = 'ecSecB'

# Load fitted deltaG's
protein = csv_to_protein(current_dir / 'fits' / state / f"r1_{fit_kwargs['r1']}" / 'fit_result.csv')
# Reindex such that pandas series includes uncovered N-terminal part (padded with NaN's)
dG = protein['deltaG'].reindex(np.arange(1, protein.index.max() + 1)) * 1e-3
# Convert deltaG values to colors
colors = apply_cmap(dG, rgb_cmap, rgb_norm)

cmd.load(str(pdb_path))
do_render(colors, views, fname_tetramer)

#-------------------------------------------------------------------------#
# ecSecB tetramer/dimer delta
#-------------------------------------------------------------------------#

states = ['SecB WT apo', 'SecB his dimer apo']
fit_output_file = current_dir / 'batch_fits' / 'SecB_tetramer_dimer' / f"r2_{fit_kwargs['r2']}" / f'fit_result.csv'
protein = csv_to_protein(fit_output_file)

p1, p2 = protein[states[0]], protein[states[1]]
dG = (p1['deltaG'] - p2['deltaG'])*1e-3
dG = dG.reindex(np.arange(1, dG.index.max() + 1))

colors = apply_cmap(dG, diff_cmap, diff_norm)

cmd.reinitialize()
cmd.load(str(pdb_path))
do_render(colors, views, fname_tetramer_dimer)

