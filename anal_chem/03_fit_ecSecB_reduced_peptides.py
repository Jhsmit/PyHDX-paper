import numpy as np
from dask.distributed import Client
from pyhdx.fileIO import csv_to_protein, save_fitresult, read_dynamx
from pyhdx.fitting import fit_gibbs_global, fit_rates_weighted_average
from pyhdx.local_cluster import default_cluster
from pyhdx.models import PeptideMasterTable, HDXMeasurement

from functions.base import settings_dict, input_data_dir, fitresults_dir

rng = np.random.default_rng(43)

do_guesses = False
if __name__ == '__main__':
    if do_guesses:
        cluster = default_cluster()
        client = Client(cluster)


    fit_kwargs = settings_dict['ecsecb_tetramer_dimer']
    fit_kwargs.pop('r2')

    output_dir = fitresults_dir / 'fits_ecsecb_reduced'
    output_dir.mkdir(exist_ok=True)

    state = 'ecSecB'
    fpath = input_data_dir / 'SecB apo full series.csv'
    data = read_dynamx(fpath)
    control = ('Full deuteration control', 0.167)
    temperature, pH = 273.15 + 30, 8.
    d_percentage = 90
    pmt = PeptideMasterTable(data, drop_first=1, ignore_prolines=True, remove_nan=False, d_percentage=d_percentage)

    pmt.set_control(control)
    data = pmt.get_state('SecB WT apo')
    ids = data['start'] + 1000 * data['end']

    # ids = np.array([f'{s}, e) for s, e in zip(data['start'], data['end']))])
    # sequences = data['sequence']
    Np = len(np.unique(ids))
    fracs = [0.75, 0.5, 0.25]
    for frac in fracs:
        for i in range(3):
            print(frac, i)
            fit_output_dir = output_dir / f"frac_{frac}_{i}"
            fit_output_dir.mkdir(exist_ok=True)
            n = int(np.round(frac*Np))
            selected = rng.choice(np.unique(ids), n, replace=False)
            bools = np.isin(ids, selected)
            selected_data = data.copy()[bools]
            hdxm = HDXMeasurement(selected_data, pH=pH, temperature=temperature)

            guess_name = "guesses"

            if do_guesses:
                guess_result = fit_rates_weighted_average(hdxm, client=client)
                guess_output = guess_result.output
                guess_output.to_file(fit_output_dir / (guess_name + '.csv'))
                guess_output.to_file(fit_output_dir / (guess_name + '.txt'), fmt='pprint')
            else:
                guess_output = csv_to_protein(fit_output_dir / (guess_name + '.csv'))

            gibbs_guess = hdxm.guess_deltaG(guess_output['rate'])
            fr_global = fit_gibbs_global(hdxm, gibbs_guess, **fit_kwargs)

            save_fitresult(fit_output_dir, fr_global)

            hdxm.to_file(fit_output_dir / 'hdxm.csv')

