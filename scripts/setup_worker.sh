#!/usr/bin/env bash
# Run in Linux/WSL from any directory. Models/data are provisioned separately.
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
python_bin=${PYTHON_BIN:-python3.9}
"$python_bin" -c 'import sys; assert sys.version_info[:2] == (3,9), "Historical worker requires Python 3.9"'
"$python_bin" -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
printf '\nDependencies installed. Supply the assets listed in README.md before starting the worker.\n'
