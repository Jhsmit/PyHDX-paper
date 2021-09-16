from functions.logging import write_log
from functions.base import settings_dict, data_dict, states, input_data_dir, current_dir, fitresults_dir
from pyhdx.batch_processing import yaml_to_hdxm
from pyhdx.fileIO import csv_to_protein, save_fitresult, dataframe_to_file
from pyhdx.fitting import fit_gibbs_global
from pyhdx.fitting_torch import CheckPoint
import time

write_log(__file__)

output_dir = fitresults_dir / 'single_fits'
output_dir.mkdir(exist_ok=True)


def load_and_fit(fit_output_dir, fit_kwargs):
    hdxm = yaml_to_hdxm(data_dict[state], data_dir=input_data_dir)
    guesses = csv_to_protein(current_dir / 'guesses' / f'{state}_initial_guess.csv')['rate']
    gibbs_guess = hdxm.guess_deltaG(guesses)

    t0 = time.time()
    checkpoint = CheckPoint(epoch_step=1000)
    fr = fit_gibbs_global(hdxm, gibbs_guess, callbacks=[checkpoint], **fit_kwargs)
    t1 = time.time()

    final_output_dir = fit_output_dir / f"r1_{fit_kwargs['r1']}"
    print(final_output_dir)
    log_lines = [f"Time elapsed: {(t1-t0):.2f} s"]
    save_fitresult(final_output_dir, fr, log_lines=log_lines)
    history_df = checkpoint.to_dataframe()
    dataframe_to_file(final_output_dir / 'model_history.csv', history_df)
    dataframe_to_file(final_output_dir / 'model_history.txt', history_df, fmt='pprint')


rainbow_fit_kwargs = {
    'r1': 1,
    **settings_dict['single_fit']
}


# # Individual fits for rainbowclouds plot
for state in states:
    print(state)
    #rainbow_fit_kwargs['epochs'] = 20
    load_and_fit(output_dir / state, rainbow_fit_kwargs)


# Single ecSecB fit for figure 1e, supplementary figure 3, supplementary figure 5
state = 'ecSecB'
output_dir = fitresults_dir / 'ecSecB_r1'
output_dir.mkdir(exist_ok=True)
ecsecb_fit_kwargs = settings_dict['ecsecb_tetramer_dimer']
ecsecb_fit_kwargs.pop('r2')
for r1 in [0, 0.05, 0.2, 0.5, 2, 10, 50, 100]:
    fit_kwargs = ecsecb_fit_kwargs
    fit_kwargs['r1'] = r1

    load_and_fit(output_dir, fit_kwargs)
