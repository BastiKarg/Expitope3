#!/usr/bin/env bash
# Build a user-supplied, licensed Rosetta source tree on Linux/WSL.
set -euo pipefail
source_dir=${1:?Usage: bash scripts/setup_rosetta.sh /path/to/rosetta/main/source [build_jobs]}
jobs=${2:-2}
[[ "$jobs" =~ ^[1-9][0-9]*$ ]] || { echo 'build_jobs must be positive' >&2; exit 2; }
[[ -f "$source_dir/scons.py" ]] || { echo 'Expected main/source/scons.py' >&2; exit 2; }
command -v python3 >/dev/null
command -v g++ >/dev/null
cd -- "$source_dir"
python3 scons.py -j"$jobs" mode=release relax docking_prepack_protocol docking_protocol
printf '\nBuilt core applications. Locate binaries in %s/bin and database in ../database.\n' "$PWD"
printf 'This does not install or validate the historical mutation/alanine-scan pipeline.\n'
