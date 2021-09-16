# from pyhdx.support import pprint_df_to_file
# from pathlib import Path
# import pandas as pd
# import yaml
# from datetime import datetime
# from pyhdx import VERSION_STRING
#
#
# def save_fitresult(fr, fit_kwargs, output_dir, time_elapsed=None):
#     output_dir.mkdir(parents=True, exist_ok=True)
#     output_dir = Path(output_dir)
#
#     yaml_file_out = output_dir / 'settings.yaml'
#     yaml_file_out.write_text(yaml.dump(fit_kwargs))
#
#     fit_file_out = output_dir / 'deltaG.csv'
#     fr.output.to_file(fit_file_out)
#
#     fit_file_out_pprint = output_dir / 'deltaG.txt'
#     fr.output.to_file(fit_file_out_pprint, fmt='pprint')
#
#     fr.losses.to_csv(output_dir / 'losses.csv')
#
#     pprint_df_to_file(fr.losses, output_dir / 'losses.txt')
#
#     loss = f'Total_loss {fr.total_loss:.2f}, mse_loss {fr.mse_loss:.2f}, reg_loss {fr.reg_loss:.2f}' \
#            f'({fr.regularization_percentage:.2f}%)'
#     epochs = f"Number of epochs: {len(fr.metadata['total_loss'])}"
#     version = VERSION_STRING
#     now = datetime.now()
#     date = f'# {now.strftime("%Y/%m/%d %H:%M:%S")} ({int(now.timestamp())})'
#
#     lines = [date, version, loss, epochs]
#     if time_elapsed is not None:
#         lines.append(f"Time elapsed: {(time_elapsed):.2f} s")
#     log_file_out = output_dir / 'log.txt'
#     log_file_out.write_text('\n'.join(lines))