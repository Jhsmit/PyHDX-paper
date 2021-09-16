import time

from pyhdx.batch_processing import yaml_to_hdxm
from pyhdx.fileIO import csv_to_protein, save_fitresult, dataframe_to_file
from pyhdx.fitting import fit_gibbs_global_batch_aligned, fit_gibbs_global
from pyhdx.models import HDXMeasurementSet
from pyhdx.support import pprint_df_to_file
from pyhdx.fitting_torch import CheckPoint


from functions.align import alignments
from functions.base import *

from functions.logging import write_log


write_log(__file__)


fit_kwargs = settings_dict['ecsecb_mtsecb']

output_base_dir = fitresults_dir / 'batch_fits' / 'ecSecB_mtSecB'
output_base_dir.mkdir(parents=True, exist_ok=True)

output_base_dir_single = fitresults_dir / 'ecSecB_mtSecB_single'
output_base_dir_single.mkdir(parents=True, exist_ok=True)


# Load the HDX-MS data and initial guesses
states = ['ecSecB', 'mtSecB']
hdxm_list = [yaml_to_hdxm(data_dict[state], data_dir=input_data_dir) for state in states]
hdx_set = HDXMeasurementSet(hdxm_list)

guesses = [csv_to_protein(f'guesses/{state}_initial_guess.csv')['rate'] for state in states]
gibbs_guess = hdx_set.guess_deltaG(guesses)

for output_folder, alignment in alignments.items():
    output_dir = output_base_dir / output_folder
    output_dir.mkdir(parents=True, exist_ok=True)
    hdx_set.add_alignment(list(alignment.values()))

    sequence_alignment = hdx_set.aligned_dataframes['sequence']
    pprint_df_to_file(sequence_alignment, output_dir / 'sequence_alignment.txt')
    t0 = time.time()
    checkpoint = CheckPoint(epoch_step=1000)
    fr = fit_gibbs_global_batch_aligned(hdx_set, gibbs_guess, callbacks=[checkpoint], **fit_kwargs)
    t1 = time.time()

    log_lines = [f"Time elapsed: {(t1 - t0):.2f} s"]
    save_fitresult(output_dir, fr, log_lines=log_lines)

    history_df = checkpoint.to_dataframe(names=hdx_set.names)
    dataframe_to_file(output_dir / 'model_history.csv', history_df)
    dataframe_to_file(output_dir / 'model_history.txt', history_df, fmt='pprint')

# individual fits
fit_kwargs.pop('r2')
for state, hdxm, guess in zip(states, hdxm_list, guesses):
    gibbs_guess = hdxm.guess_deltaG(guess)

    t0 = time.time()
    checkpoint = CheckPoint(epoch_step=1000)
    fr = fit_gibbs_global(hdxm, gibbs_guess, callbacks=[checkpoint], **fit_kwargs)
    t1 = time.time()

    log_lines = [f"Time elapsed: {(t1 - t0):.2f} s"]
    output_dir = output_base_dir_single / state
    save_fitresult(output_dir, fr, log_lines=log_lines)

    history_df = checkpoint.to_dataframe(names=hdx_set.names)
    dataframe_to_file(output_dir / 'model_history.csv', history_df)
    dataframe_to_file(output_dir / 'model_history.txt', history_df, fmt='pprint')