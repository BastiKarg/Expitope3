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
