import numpy as np
import pandas as pd


def nma_result_to_df(nma_path, chains=None):
    modes_pdb = nma_path / 'modes_CA.pdb'

    text = modes_pdb.read_text()
    lines = text.split('\n')

    residue_numbers = np.array([int([entry for entry in line.split(' ') if entry][5]) for line in lines if line.startswith('ATOM')])
    repeated_idx = np.where(np.diff(residue_numbers) == 0)[0]

    if chains is None:
        chains = np.array([[entry for entry in line.split(' ') if entry][4] for line in lines if line.startswith('ATOM')])
        chains = np.unique(chains)

    if repeated_idx.size > 0:
        residue_numbers = np.delete(residue_numbers, repeated_idx)

    all_indices = [0] + list(np.where(np.diff(residue_numbers) < 0)[0] + 1) + [len(residue_numbers)]
    split_residues = [residue_numbers[i0: i1] for i0, i1 in zip(all_indices[:-1], all_indices[1:])]

    series_list = []
    for i in range(6):
        displacement = np.genfromtxt(nma_path / 'displacement' / 'displacements.txt', skip_header=1)[:, i + 1]
        displacement = np.sqrt(displacement) # Take square root as displacements are squared displacements
        mode = i + 7

        split_displacement = [displacement[i0: i1] for i0, i1 in zip(all_indices[:-1], all_indices[1:])]
        for r, d, chain in zip(split_residues, split_displacement, chains):
            index = pd.Index(r, name='r_number')
            series = pd.Series(d, name=f"displacement_mode_{mode}_chain_{chain}", index=index)
            series_list.append(series)

    df = pd.concat(series_list, axis=1)

    # Sum all chains
    cols = [col for col in df.columns if f'chain_{chains[0]}' in col]
    series_list = []  # List of summed displacements per mode
    for entry in cols:
        base = entry[:-1]
        selection = [base + letter for letter in chains]

        name = base[:-7]
        sub_df = df[selection]
        series = sub_df.mean(axis=1) # Average of all chains
        series.name = name
        series_list.append(series)
    df_chains = pd.concat(series_list, axis=1)

    # Sum all modes
    overall_sum = df_chains.sum(axis=1)
    overall_sum.name = 'displacement'

    total = pd.concat([df, df_chains, overall_sum], axis=1)

    return total