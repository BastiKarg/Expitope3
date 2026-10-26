# Shiny restoration (work in progress)

The Windows UI is based on `app_local.R`. All three supplied variants parse,
but this variant contains the most useful Windows filesystem adaptations.
This is a functional selection, not a claim about which file was last edited.
Originals are retained locally in `.cleanup-backup/shiny-originals/`.

## Launch

From the repository root, run `Rscript scripts/install_shiny_dependencies.R`,
then `Rscript scripts/run_shiny.R`. Open http://127.0.0.1:3838.
The server binds only to localhost. R 4.5.1 and existing user-library packages
were used for the restoration checks; some packages were built for R 4.5.2.

`app.R` loads `R/bootstrap.R`. UI, server, legacy plotting helpers, session
handlers and filesystem helpers are separated under `R/`.
Nineteen unmounted output handlers were removed. The historically disabled
tumor/upload tabs are not restored features.

## Local assets and current validation

Update after the completed copy: all three UI/help images, tissue metadata,
abbreviations and the historical reference database folder are present.
`tcga_read_counts/` contains 33 TPM tables (about 3.3 GB). The app discovers these
without copying them and derives 60,660 gene mappings from their annotated first
column when the standalone gene metadata files are absent. Raw-count tables in
the same folder are not mixed into the TPM explorer. Existing sample column
labels are preserved; labels such as `TCGA-ACC_1` are not verified patient IDs.

The ACC test confirms 79 sample values and the TSPAN6 mapping. All 18 combinations
of boxplot/violin/beeswarm, raw/log2 and ordering options render in a Shiny test
session with that gene and tissue. This is a representative functional test,
not a content-integrity audit of every copied table. The plot code shares a
single transformation: `log2(expression + 1)` for the log display.

Still missing for a complete restoration: sample-level GTEx tissue tables and
AlphaFold model weights. The GTEx files under the reference database folder
contain tissue-level columns and cannot substitute for individual samples.
The copied `kanu/alphafold_finetune/` contains outputs, not a complete fine-tuning
installation. Rosetta is not installed; see [Rosetta setup](ROSETTA_SETUP.md).
The old script points to `/data/project/kanu/alphafold_download_dir` as the
AlphaFold data directory; search the original storage for its `params/` child.

Expected legacy filenames (standalone gene metadata are now optional):

- `ensg_list`: ENSG identifiers.
- `ensg_genename.tsv`: columns `ENSG`, `gene`.
- `tcga_gtex_tissues.tsv`: tissue/body mapping used by the expression views.
- `tcga_abbreviations`: the original tab-separated abbreviation table.
- `tcga/`: sample-level expression tables with `id` and sample columns.
- `Datenbanken_Ordner_shiny/`: original visualization distributions.
- `www/expitope_logo_v3.png`, `www/exp_scoring_new.png`,
  `www/model_score_table.png`: original UI/help images.

Missing data are reported; empty metadata are not synthetic expression data.
Full plot/export validation needs the original assets and representative
completed Expitope runs. Keep these large datasets outside Git.

Environment variables configure paths: `EXPITOPE_DATA_DIR`, `EXPITOPE_EXPRESSION_DIR`,
`EXPITOPE_STRUCTURE_DIR`, `EXPITOPE_QUEUE_DIR`, `EXPITOPE_RESULTS_DIR`.
Defaults are the repository root, `3d/testing_2`, `running`, and `results`.
The expression directory defaults to `tcga/` when present, otherwise
`tcga_read_counts/`. It reads `*_tpm.tsv` files preferentially.
Windows and the Linux worker must see the same job/result directories, with
their respective Windows and Linux path mappings.

## Structure worker restoration

The supplied pipeline is **TCRdock plus Rosetta**, not interchangeable with a
generic AlphaFold or ColabFold invocation. The local TCRdock checkout identifies
commit `758ffb7e1d784a05b488f2720dc20ff534f64365`.
The [upstream TCRdock instructions](https://github.com/phbradley/TCRdock)
describe a later AlphaFold 2.3.2 update. Select and pin a tested TCRdock revision
together with its compatible AlphaFold code and weights; do not mix revisions
without a reference prediction comparison.

Before enabling new structure jobs:

1. Locate the original AlphaFold `params/`, any fine-tuned weights and actual
   fine-tuning source, a working Rosetta installation, and the old environment
   specification if available. The supplied fine-tuning directory primarily
   contains outputs; these do not establish a working installation.
2. Verify WSL/Linux GPU access and compatible CUDA/JAX dependencies using that
   chosen release's installation instructions. Large weights and Rosetta must
   be installed separately, subject to their own terms; do not commit them here.
3. Replace `/data/project/kanu` and `/data2/project/kanu` paths in historical
   scripts with configuration. Give each job its own working directory.
4. Replace the old queue consumer: it executes shell command files and deletes
   jobs even after failures. Record explicit queued/running/succeeded/failed
   states, preserve logs, validate input and propagate process exit codes.
5. Test one known TCRdock target, one mutation, one Rosetta relaxation/docking
   and one alanine scan against expected outputs before enabling submissions.

TCRdock target columns are `organism`, `mhc_class`, `mhc`, `peptide`, `va`, `ja`,
`cdr3a`, `vb`, `jb`, `cdr3b`. The supplied `testing_2/test.tsv` has placeholder
columns and is not a valid target template.

Structure submission is currently explicitly unavailable. Existing PDB viewing
and downloads can be tested independently. Chain B is treated as the peptide,
matching the historical pipeline; arbitrary uploaded PDBs may use another chain.
Hard-coded example biotin contacts/residue numbers were removed from the viewer.

## Verification

Run `Rscript tests/test_shiny_runtime.R`. It checks exact gene matching and
sample/value alignment, isolated prediction queue writes, invalid inputs,
mutated PDB resolution, gapped peptide residue selections, and missing metadata.
These checks do not validate model accuracy or a real GPU execution.

The four supplied Rosetta score files were also read successfully (20, 1,
170 and 3 rows). Trailing empty header fields are removed before parsing.
The local `wsl --list --quiet` check reports that WSL is not installed.

The interface is published with these limitations; full worker end-to-end
validation remains pending. No AlphaFold/Rosetta execution has been certified
by the UI smoke test.

The restored `distPlot_server_new` histogram is loaded by the normal `app.R`
entry point. It reads only the chosen reference tissue, groups hover data once,
and sends bin summaries to Plotly. The historical per-row scan is removed.
Counts retain reference/match join multiplicity. One consistent 99-bin scheme
replaces the old inconsistent tooltip/colour/histogram bins; maximum and
constant values are included. `tests/test_distribution.R` checks these cases.
