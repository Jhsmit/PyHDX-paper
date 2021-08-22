import time

from pyhdx.batch_processing import load_from_yaml
from pyhdx.fileIO import csv_to_protein, save_fitresult
from pyhdx.fitting import fit_gibbs_global_batch_aligned
from pyhdx.models import HDXMeasurementSet
from pyhdx.support import pprint_df_to_file

from functions.align import alignments
from functions.base import *

fit_kwargs = settings_dict['ecsecb_mtsecb']

output_base_dir = current_dir / 'batch_fits' / 'ecSecB_mtSecB'
output_base_dir.mkdir(parents=True, exist_ok=True)

# Load the HDX-MS data and initial guesses
states = ['ecSecB', 'mtSecB']
hdxm_list = [load_from_yaml(data_dict[state], data_dir=input_data_dir) for state in states]
hdx_set = HDXMeasurementSet(hdxm_list)

guesses = [csv_to_protein(f'guesses/{state}_initial_guess.csv')['rate'] for state in states]
gibbs_guess = hdx_set.guess_deltaG(guesses)


r_numbers = iter([[1, 1], [1, -2], [1, -2]])
for output_folder, alignment in alignments.items():
    output_dir = output_base_dir / output_folder
    output_dir.mkdir(parents=True, exist_ok=True)
    hdx_set.add_alignment(list(alignment.values()), first_r_numbers=next(r_numbers))

    sequence_alignment = hdx_set.aligned_dataframes['sequence']
    pprint_df_to_file(sequence_alignment, output_dir / 'sequence_alignment.txt')
    t0 = time.time()
    fr = fit_gibbs_global_batch_aligned(hdx_set, gibbs_guess, **fit_kwargs)
    t1 = time.time()

    log_lines = [f"Time elapsed: {(t1 - t0):.2f} s"]
    save_fitresult(output_dir, fr, log_lines=log_lines)
