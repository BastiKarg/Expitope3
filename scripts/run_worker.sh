#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
[[ "$(uname -s)" == Linux ]] || { echo 'Use Linux or WSL for the worker.' >&2; exit 2; }
venv_dir=${VENV_DIR:-.venv}
[[ -x "$venv_dir/bin/python" ]] || { echo 'Run bash scripts/setup_worker.sh first.' >&2; exit 2; }
exec "$venv_dir/bin/python" -u main_final.py "$@"
