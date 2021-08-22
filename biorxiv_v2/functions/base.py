import numpy as np
import torch
from pathlib import Path
import yaml

np.random.seed(43)
torch.manual_seed(43)

root_dir = Path(__file__).parent.parent.parent
input_data_dir = root_dir / 'data'
current_dir = root_dir / 'biorxiv_v2'

yaml_stream = (root_dir / 'datasets.yaml').read_text()
data_dict = yaml.safe_load(yaml_stream)

settings_stream = yaml_stream = (current_dir / 'settings.yaml').read_text()
settings_dict = yaml.safe_load(settings_stream)


cluster = '127.0.0.1:52123'

states = ['ecSecB',
          'ecSecB_dimer',
          'SecA_monomer',
          'SecA_monomer_ADP',
          'SecA_1-901_wt_apo',
          'SecA_wt_ADP',
          'SecA_1-834_apo',
          'SecA_1-834_ADP',
          'mtSecB',
          'TF_monomer_4C',
          'bcl2',
          'hPREP_default',
          'EscV',
          'PPiB_WT',
          'PPiA_WT',
          'MBP_wt',
          'TF_dimer_4C',
          ]
