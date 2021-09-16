from pathlib import Path
from pyhdx import VERSION_STRING
from importlib.metadata import version
from datetime import datetime
import sys
import os
import platform


def write_log(script_path):
    log_name = Path(script_path).stem
    out_path = Path(script_path).parent / 'logs'
    out_path.mkdir(exist_ok=True)

    lines = [f"Log file for python script: {log_name}.py"]

    now = datetime.now()
    date = f'{now.strftime("%Y/%m/%d %H:%M:%S")} ({int(now.timestamp())})'
    lines.append(f"Executed at: {date}")

    lines.append(f"Python version: {sys.version}")
    lines.append(f"OS: {platform.system()}, release {platform.release()}")
    try:
        lines.append(f"Conda env: {os.environ['CONDA_DEFAULT_ENV']}")
    except KeyError:
        pass

    lines.append("")
    lines.append(VERSION_STRING)

    lines.append("")
    lines.append("Dependencies versions:")
    packages = ['numpy', 'torch', 'pandas', 'hdxrate', 'scipy', 'symfit', 'scikit-image', 'dask', 'distributed']
    for package in packages:
        ver = version(package)
        line = f"{package}: {ver}"
        lines.append(line)

    s = '\n'.join(lines)

    (out_path / (log_name + '.log')).write_text(s)




















if __name__ == '__main__':
    write_log(__file__)