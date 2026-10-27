#!/usr/bin/env bash
# Run in Linux/WSL from any directory. Models/data are provisioned separately.
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
python_bin=${PYTHON_BIN:-python3}
venv_dir=${VENV_DIR:-.venv}
"$python_bin" -c 'import sys; assert sys.version_info[:2] in ((3,9), (3,10)), "Use Python 3.9 or 3.10"'
"$python_bin" -m venv "$venv_dir"
if [[ ${EXPITOPE_CPU_ONLY:-0} == 1 ]]; then
    "$venv_dir/bin/python" -m pip install 'torch==1.13.1+cpu' --index-url https://download.pytorch.org/whl/cpu
fi
"$venv_dir/bin/python" -m pip install -r requirements.txt
printf '\nDependencies installed. Supply the assets listed in README.md before starting the worker.\n'
