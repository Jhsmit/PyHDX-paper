from functions.logging import write_log
from functions.base import settings_dict, data_dict, states, input_data_dir, current_dir
from pyhdx.batch_processing import load_from_yaml
from pyhdx.fileIO import csv_to_protein, save_fitresult
from pyhdx.fitting import fit_gibbs_global
import time

write_log(__file__)

output_dir = current_dir / 'fits'
output_dir.mkdir(exist_ok=True)


def load_and_fit(state, fit_kwargs):
    hdxm = load_from_yaml(data_dict[state], data_dir=input_data_dir)
    guesses = csv_to_protein(current_dir / 'guesses' / f'{state}_initial_guess.csv')['rate']
    gibbs_guess = hdxm.guess_deltaG(guesses)

    t0 = time.time()
    fr = fit_gibbs_global(hdxm, gibbs_guess, **fit_kwargs)
    t1 = time.time()

    fit_output_dir = output_dir / state / f"r1_{fit_kwargs['r1']}"

    log_lines = [f"Time elapsed: {(t1-t0):.2f} s"]
    save_fitresult(fit_output_dir, fr, log_lines=log_lines)


rainbow_fit_kwargs = {
    'r1': 0.5,
    **settings_dict['single_fit']
}


# Individual fits for rainbowclouds plot
for state in states:
    print(state)
    load_and_fit(state, rainbow_fit_kwargs)

# Single ecSecB fit for figure 1e, supplementary figure 3
state = 'ecSecB'
ecsecb_fit_kwargs = settings_dict['ecsecb_tetramer_dimer']
ecsecb_fit_kwargs.pop('r2')
for r1 in [0, 0.05, 0.1, 0.5, 1, 2, 5, 10]:
    fit_kwargs = {'r1': r1, **settings_dict['single_fit']}

    load_and_fit(state, fit_kwargs)

# Individual mtSecB for SI fig 5 (top) and SI fig 9
ecsecb_fit_kwargs = settings_dict['ecsecb_tetramer_dimer']
ecsecb_fit_kwargs.pop('r2', None)
r1 = 0.05
fit_kwargs = {'r1': r1, **settings_dict['single_fit']}
for state in ['mtSecB', 'MBP_wt', 'PPiA_WT', 'PPiB_WT']:
    print(state)
    load_and_fit(state, fit_kwargs)