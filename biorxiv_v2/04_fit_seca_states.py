from functions.logging import write_log
from functions.base import settings_dict, data_dict, input_data_dir, current_dir
from pyhdx.batch_processing import load_from_yaml
from pyhdx.models import HDXMeasurementSet
from pyhdx.fileIO import csv_to_protein, save_fitresult
from pyhdx.fitting import fit_gibbs_global_batch

import time

write_log(__file__)

output_dir = current_dir / 'batch_fits' / 'SecA_states'
output_dir.mkdir(parents=True, exist_ok=True)

fit_kwargs = settings_dict['seca_all_states']

states = ['SecA_monomer',
          'SecA_monomer_ADP',
          'SecA_1-834_ADP',
          'SecA_wt_ADP',
          'SecA_1-834_apo',
          'SecA_1-901_wt_apo',
          ]

hdxm_list = [load_from_yaml(data_dict[state], data_dir=input_data_dir) for state in states]
hdx_set = HDXMeasurementSet(hdxm_list)

guesses = [csv_to_protein(current_dir / 'guesses' / f'{state}_initial_guess.csv')['rate'] for state in states]
gibbs_guess = hdx_set.guess_deltaG(guesses)

t0 = time.time()
fr = fit_gibbs_global_batch(hdx_set, gibbs_guess, **fit_kwargs)
t1 = time.time()

log_lines = [f"Time elapsed: {(t1 - t0):.2f} s"]
save_fitresult(output_dir, fr, log_lines=log_lines)