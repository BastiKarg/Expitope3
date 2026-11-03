# Running the interface and computation worker

There are two independent processes. The Shiny interface writes complete job
files into `running/`. The persistent Python worker reads them, performs the
prediction, writes `results/<project>/` and updates `recent_finish/last_finish`.
Rosetta/TCRdock is a separate optional pipeline and is not needed for Expitope
predictions or viewing existing results.

## View existing results (Windows)

Install R separately, then run these commands from the repository root:

```powershell
Rscript scripts/install_shiny_dependencies.R
Rscript scripts/run_shiny.R
```

Open http://127.0.0.1:3838. Select a finished run and a result file. A Python
worker, CUDA and Rosetta are not required to inspect existing output.
Rscript must be on PATH, or use its full installed path. RStudio can also run
`shiny::runApp('.')`. RStudio Server is not required to serve the interface.

## Run new predictions (Linux/WSL)

Use the same checkout and shared data directories as the Windows UI, e.g.
`C:\Users\cakes\Documents\chatty` maps to
`/mnt/c/Users/cakes/Documents/chatty` in WSL. Use Python 3.9 or 3.10 in Linux, then:

```bash
bash scripts/setup_worker.sh
bash scripts/run_worker.sh
```

Supply the database, precomputed cleavage tables, checkpoints and model caches
listed in the root README before starting. The installer installs Python
packages; it does not provision research data or silently download all model
weights. Keep the worker terminal running while submitting from Shiny. Its
startup loads models and can take time before it consumes jobs.

For this first deployment, use the default `running/`, `results/` and
`recent_finish/` directories in the shared checkout. UI environment overrides
alone do not reconfigure the existing Python worker's hard-coded paths.
Do not start multiple workers on the same queue: job claiming is not atomic.

The Python requirements are reconstructed historical versions, not a verified
GPU environment lock. CUDA, TensorFlow/PyTorch coexistence and one full reference
prediction must be validated on the target Linux host before calling the full
installation reproducible. The launcher has been exercised under WSL; the
environment was installed directly, rather than through the setup wrapper.

### Local WSL CPU validation environment

Ubuntu 22.04 / WSL2 is installed locally. The isolated environment is
`/opt/expitope3-venv` (Python 3.10, historical dependencies, PyTorch 1.13.1 CPU).
Its model caches belong to the WSL root account, used for this initial test.
GPU passthrough is detected, but GPU model inference is not validated.

Update 2026-11-03: ProtT5 GPU inference has now been tested separately (see below).
The full worker remains in the validated CPU environment.

On 2026-10-27, a Shiny-generated queue file for the synthetic custom sequence
`MKTIIALSYIFCLVFADYKDDDDK` and peptide `ALSYIFCLV` completed with real cleavage,
MHCflurry, TAP, ProtT5 and classifier computations on CPU. Output is local under
`results/wsl_smoke_20261027/Mutation_epitope_search_scores_mlp.tsv`; the job moved
to `recent_finish/last_finish`. This validates the custom-sequence queue route,
not the database-wide search/expression workflow or scientific accuracy.

To process exactly one existing job without consuming the rest of the queue:

```bash
cd /mnt/c/Users/cakes/Documents/chatty
CUDA_VISIBLE_DEVICES=-1 VENV_DIR=/opt/expitope3-venv \
  bash scripts/run_worker.sh --job running/YOUR_JOB.job
```

Run this in the same WSL account/environment used to install the models.
`--once` instead processes one queued job and exits. Without either option,
the worker remains active. A successful job is moved to `recent_finish/last_finish`.

The current FASTA is `database5/Expitope3_final_clean_database_singlev3.fa_new_clean`.
Missing/empty NN cleavage variants 4 and 5 explicitly fall back to the default
table, as requested by the dataset owner. Tables are loaded only when selected.
The supplied precomputed cleavage tables currently end at protein ID 329469,
whereas the current FASTA has 315451 entries: database-wide jobs need matching
tables and matching protein order. An ID count check rejects this known mismatch;
equal counts alone would not prove correct ordering. Custom-sequence jobs compute
cleavage directly and do not use these tables.

### Isolated ProtT5 GPU check (2026-11-03)

The separate `/opt/expitope3-gpu-test` environment uses Python 3.10,
PyTorch 2.7.1+cu128, transformers 4.29.2, sentencepiece 0.1.99 and numpy 1.24.3.
It reuses the same root-account model cache. PyTorch's
[2.7 release](https://pytorch.org/blog/pytorch-2-7/) introduced Blackwell support
with CUDA 12.8. The historical CPU environment is unchanged.

Reproduce from the repository root under the same WSL account:

```bash
CUDA_VISIBLE_DEVICES=-1 /opt/expitope3-venv/bin/python scripts/check_embedding_device.py \
  --device cpu --output .cleanup-backup/embedding-cpu.npz
/opt/expitope3-gpu-test/bin/python scripts/check_embedding_device.py \
  --device cuda --reference .cleanup-backup/embedding-cpu.npz \
  --output .cleanup-backup/embedding-gpu.npz
```

The script calls the actual application embedding helper on six peptides,
checks all 1024 components per peptide, and reports the median of three warmed
inference runs. Model loading and device transfer are excluded. The observed
CPU time was 5.008 s (PyTorch 1.13.1, four CPU threads); the RTX 5070 Ti Laptop
GPU took 0.147 s (PyTorch 2.7.1). Peak allocated GPU memory was 4.84 GB.
Maximum absolute difference was 2.65e-7; all components passed comparison
at rtol=1e-3, atol=1e-4. TF32 was disabled. These numbers compare the two tested
environments, not GPU versus CPU under one identical software version.

This validates the ProtT5 stage only. The isolated GPU environment does not
contain the full worker dependencies. TensorFlow/MHCflurry coexistence, classifier
loading and a complete GPU queue run must be tested before changing the normal
worker environment. The FASTA/cleavage-table matching requirement still applies.

## GitHub distribution

Publish source, setup/launch scripts, tests and installation instructions.
Distribute large research assets separately with versions, checksums and
provenance. Keep private/user demo results out of Git unless explicitly selected
and suitable for redistribution. A small redistributable example dataset would
allow a documented demo mode without installing prediction models.

For an eventual server deployment, run Shiny and the Python worker as separate
managed services with a shared data volume and logs. Add containers or service
units after the Linux environment has passed the full prediction test; an
untested container would not resolve the current dependency uncertainty.
