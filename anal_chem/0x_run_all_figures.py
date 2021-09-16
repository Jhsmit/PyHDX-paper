import importlib
from functions.base import current_dir

print(current_dir)

scripts = [f.stem for f in current_dir.iterdir() if f.stem.startswith('fig') and f.suffix == '.py']
print(scripts)

for s in scripts:
    print(s)
    importlib.import_module(s)