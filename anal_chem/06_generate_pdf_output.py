
from pyhdx.fileIO import load_fitresult
from pyhdx.output import Output, Report

from functions.base import settings_dict, fitresults_dir
from functions.formatting import *
from functions.logging import write_log

write_log(__file__)

output_dir = current_dir / 'figures' / 'Supplement'
output_dir.mkdir(parents=True, exist_ok=True)
fname = 'Fig_SI_4_peptides_fit'

fit_kwargs = settings_dict['ecsecb_tetramer_dimer']


f = 'ecSecB_r1'
fit_dir = fitresults_dir / f / f"r1_{fit_kwargs['r1']}"

fit_result = load_fitresult(fit_dir)

output = Output(fit_result)
report = Report(output, title='Supporting Information ecSecB')
report.add_peptide_figures()
report.generate_pdf(output_dir / 'SecB_fit_report')

# Adapted code from add_peptide_figures to export each subplot figure separately
# (Figure S3)
Np = report.output.fit_result.data_obj.Np
indices = range(Np)
ncols = 4
nrows = 5
n = ncols * nrows
chunks = [indices[i:i + n] for i in range(0, len(indices), n)]
for i, chunk in enumerate(chunks):
    print(i)

    fig = report.output._make_peptide_subplots(chunk, ncols=ncols, nrows=nrows)

    plt.savefig(output_dir / f'{fname}_{i}.pdf')
    plt.savefig(output_dir / f'{fname}_{i}.png')

    plt.close(fig)