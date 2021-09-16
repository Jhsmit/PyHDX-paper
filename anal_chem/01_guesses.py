from functions.logging import write_log

from functions.base import current_dir, input_data_dir, data_dict, states
from pyhdx.local_cluster import default_cluster, default_client
from pyhdx.fitting import fit_rates_weighted_average
from pyhdx.batch_processing import yaml_to_hdxm
from dask.distributed import Client

write_log(__file__)
output_dir = current_dir / 'guesses'
output_dir.mkdir(exist_ok=True)

if __name__ == '__main__':
    client = Client()
    print(client)

    info_dir = current_dir / 'info'
    info_dir.mkdir(exist_ok=True)

    for state in states:
        print(state)
        hdxm = yaml_to_hdxm(data_dict[state], data_dir=input_data_dir)
        guesses = fit_rates_weighted_average(hdxm, client=client)
        guesses.output.to_file(current_dir / 'guesses' / f'{state}_initial_guess.csv')
        guesses.output.to_file(current_dir / 'guesses' / f'{state}_initial_guess.txt', fmt='pprint')
        hdxm.coverage.protein.to_file(info_dir / f'{state}_info.csv')
        hdxm.coverage.protein.to_file(info_dir / f'{state}_info.txt', fmt='pprint')

