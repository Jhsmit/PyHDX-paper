from functions.logging import write_log
from functions.base import settings_dict, data_dict, input_data_dir, current_dir
from pyhdx.batch_processing import load_from_yaml
from pyhdx.models import HDXMeasurementSet
from pyhdx.fileIO import csv_to_protein, save_fitresult
from pyhdx.fitting import fit_gibbs_global_batch

import time

write_log(__file__)

output_dir = current_dir / 'batch_fits' / 'SecB_tetramer_dimer'
output_dir.mkdir(parents=True, exist_ok=True)

fit_kwargs = settings_dict['ecsecb_tetramer_dimer']

states = ['ecSecB', 'ecSecB_dimer']

hdxm_list = [load_from_yaml(data_dict[state], data_dir=input_data_dir) for state in states]
hdx_set = HDXMeasurementSet(hdxm_list)

guesses = [csv_to_protein(current_dir / 'guesses' / f'{state}_initial_guess.csv')['rate'] for state in states]
gibbs_guess = hdx_set.guess_deltaG(guesses)

for r2 in [0, 0.05, 0.1, 0.5, 1, 2, 5, 10]:
    fit_kwargs['r2'] = r2

    t0 = time.time()
    fr = fit_gibbs_global_batch(hdx_set, gibbs_guess, **fit_kwargs)
    t1 = time.time()

    fit_output_dir = output_dir / f"r2_{fit_kwargs['r2']}"

    log_lines = [f"Time elapsed: {(t1-t0):.2f} s"]
    save_fitresult(fit_output_dir, fr, log_lines=log_lines)
