from pyhdx.fileIO import csv_to_protein
from pyhdx.support import color_pymol, apply_cmap
from pymol import cmd

from functions.base import *
from functions.formatting import *
from functions.pymol import *
from functions.logging import write_log


write_log(__file__)

output = 'save'
output_dir = current_dir / 'figures' / 'Supplement' / 'SecB_structures'
output_dir.mkdir(parents=True, exist_ok=True)
fname = 'Fig_SI_6_aligned_secb'

states = ['ecSecB', 'mtSecB']
color_dict = {'pdark_blue': [34, 34, 85], 'pdark_green': [34, 85, 34]}

fit_kwargs = settings_dict['ecsecb_mtsecb']

sclf = 0.5
print("SET SCAL FACETOR!!")
px, py = 640, 480
px *= sclf
py *= sclf

alignment_folder = 'secondary_structure'
df = csv_to_protein(current_dir / 'batch_fits' / 'ecSecB_mtSecB' / alignment_folder / f'fit_result.csv')
df.rename(columns={'SecB WT apo': 'ecSecB', 'SecB_apo': 'mtSecB'}, inplace=True)

pse_file = 'ecSecB_mtSecB_aligned.pse'
pdb_path = current_dir / 'structures' / pse_file
cmd.load(str(pdb_path))

for name, c in color_dict.items():
    cmd.set_color(name, c)
cmd.set('antialias', 2)

for orientation, view in zip(['top', 'side'], aligned_views):
    cmd.set_view(view)

    for name, color in zip(states, color_dict.keys()):
        cmd.color(color, selection=f'model {name.lower()}')

    cmd.ray(px, py, renderer=0, antialias=2)
    cmd.png(output_dir / f'{fname}_{orientation}_no_gibbs.png')

    for name in states:
        dG = df[name]['deltaG'] * 1e-3
        colors = apply_cmap(dG, rgb_cmap, rgb_norm)
        color_pymol(colors, cmd, model=name.lower())

    cmd.ray(px, py, renderer=0, antialias=2)
    cmd.png(output_dir / f'{fname}_{orientation}_gibbs.png')

    cmd.disable('mtSecB')
    cmd.ray(px, py, renderer=0, antialias=2)
    cmd.png(output_dir / f'{fname}_{orientation}_gibbs_ecsecb.png')

    cmd.disable('ecSecB')
    cmd.enable('mtSecB')
    cmd.ray(px, py, renderer=0, antialias=2)
    cmd.png(output_dir / f'{fname}_{orientation}_gibbs_mtsecb.png')

    cmd.enable('ecSecB')


