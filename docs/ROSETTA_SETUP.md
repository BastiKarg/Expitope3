# Rosetta setup for the Linux worker

Rosetta is installed separately; this repository does not include its source,
binaries or database. Obtain a release under the applicable terms through the
[official download page](https://rosettacommons.org/software/download/).
Follow the [official build guide](https://docs.rosettacommons.org/docs/latest/build_documentation/Build-Documentation)
for your selected release. Record its version and compiler in your environment.

On Ubuntu/WSL, install build tools in your Linux environment:

```bash
sudo apt-get update
sudo apt-get install build-essential python3 zlib1g-dev
bash scripts/setup_rosetta.sh /path/to/rosetta/main/source 2
```

The helper builds `relax`, `docking_prepack_protocol` and `docking_protocol`
from an already downloaded source tree. It does not accept licenses, download
Rosetta, modify the operating system or install AlphaFold. Adjust the build job
count to available RAM/CPU. Run each resulting executable with `-help` and set
`-database /path/to/rosetta/main/database` for a smoke test.

The copied scripts reference both the 2021.16.61629 bundle and release-371,
including MPI builds and a 126-process/10,000-structure workload. These are
historical assumptions, not recommended defaults for a workstation. Start with
one structure in an isolated directory; compare to a known original result.
Do not launch the legacy queue consumer against user jobs.

Mutation and alanine scans need additional restoration: identify the exact
mutation executable or protocol, chain naming, residue numbering, relaxation
flags and docking partners used for each workflow. Merely installing the three
core binaries does not make those workflows operational. MPI requires a
separate MPI toolchain/build as described in the official guide.

Status: setup helper prepared, but no Linux/WSL Rosetta build or scientific
end-to-end run has been executed on this Windows machine.

The copied `kanu/3d/rosetta_bin_linux_2021.16.61629_bundle/main/` currently
contains documentation/tests/tools but lacks `source/bin/` and `database/`.
Its directory name alone does not indicate a complete executable installation.
