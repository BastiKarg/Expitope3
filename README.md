# Expitope3

Python research pipeline for searching protein sequences for potential epitope
cross-reactivity and combining cleavage, MHC presentation, TAP transport,
ProtT5 embeddings, neural-network classification, and tissue expression scores.

This repository contains the Python backend reachable from `main_final.py`.
The R application, historical benchmarks, duplicate entry points, databases,
trained checkpoints, and generated results are outside this source release.

## Pipeline

`main_final.py` runs a persistent file-queue worker. It loads the reference
proteome, cleavage tables, MHCflurry predictor, ProtT5 encoder, and MLP checkpoints
once, then processes jobs placed in `running/`.

Supported branches include whole-epitope searches, fixed amino-acid positions,
indel/multigap searches, k-mer searches, and matching against a custom protein
sequence. The code retains the original scoring formulas and output schemas.

## Environment

The historical deployment used **Linux and Python 3.9**. The multiprocessing
workers rely on inherited globals and the `fork` start method. Native Windows
and macOS default `spawn` execution are not supported by this implementation.
Use a Linux environment, for example WSL2 on Windows.

```bash
python3.9 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

`requirements.txt` reconstructs the direct runtime dependencies from the supplied
environment exports; it is not a tested lockfile. Those exports disagree on
some versions, including PyTorch Lightning. GPU/CUDA compatibility must be
validated in the deployment environment. No dependency upgrades were attempted
as part of the source cleanup.

## Required assets

Run from the repository root. Supply these assets separately, keeping the
reference-protein order consistent with the precomputed cleavage tables:

- `database5/database_october.fa`
- `database5/pepsickle_default.tsv`
- `database5/pepsickle_gbm_C.tsv` and `pepsickle_gbm_I.tsv`
- `database5/pepsickle_nn_C.tsv` and `pepsickle_nn_I.tsv`
- `best_models/best_val_f1.2_9547_113.ckpt`
- `best_models/atlas_best_val_f1.2_08846_47.ckpt`
- Tissue expression tables in `Datenbanken_Ordner_v5/`
- MHCflurry presentation models in its configured model cache
- The Hugging Face model `Rostlab/prot_t5_xl_half_uniref50-enc`

The small tissue-weight lookup `tissue_weights_exp4.unix_lowercase_unique` is
included. Large datasets and model weights are not uploaded. Their provenance,
redistribution terms, download locations, and checksums still need a separate
data release. ProtT5 downloads its weights on first use unless already cached.

## Running jobs

```bash
python main_final.py
```

The worker creates `running/` and `recent_finish/` if needed. In another shell,
write one whitespace-delimited line to a file, then move the complete file into
the queue:

```bash
printf '%s\n' 'SLLMWITQC HLA-A*02:01 0 example 1 1 normal_seq' > example.job
mv example.job running/
```

Fields, in order:

1. Peptide sequence. Mixed case marks uppercase positions as fixed (one-based).
   Entirely uppercase or lowercase peptides have no fixed positions.
2. HLA allele, such as `HLA-A*02:01`.
3. Maximum mismatches.
4. Output directory name beneath `results/`.
5. Cleavage model: `1` default, `2` GBM constitutive, `3` GBM immunoproteasome,
   `4` neural constitutive, or `5` neural immunoproteasome.
6. MLP model: `1` default, `2` atlas, or `3` (currently the same checkpoint as `1`).
7. Search mode: `normal_seq`, `indel`, `multigap`, or an integer k-mer length.
8. Optional custom protein sequence.

Queue files contain a single line without a separate header. The worker reads
the fields as pandas column names. It selects the newest queued file first and
moves each completed job to `recent_finish/last_finish`, replacing that marker.
Use trusted local job files: input validation and multi-user queue isolation are
not implemented.

Outputs are written below `results/<output-name>/`, including epitope score
tables, expression tables, log sums, and tissue/ICR summaries where applicable.
The custom-sequence path writes its own score and benchmark-format exports.

## Verification and limitations

```bash
python -m unittest discover -s tests -v
```

The lightweight tests run without downloading models. They exercise fixed
positions, sequence matching helpers, TAP/Q-score arithmetic, syntax, and
unresolved global-name detection. Full protein-search/ML inference and numerical
equivalence of complete outputs have **not** been validated in this cleanup
environment.

Known inherited limitations remain: hard-coded database offsets in the normal
search path, global multiprocessing state, queue waits without propagated worker
exceptions, inconsistent no-hit handling in the custom-sequence path, and positional
lists for score storage. These require separate behavioral changes and validation
against real reference outputs. See [cleanup notes](docs/cleanup.md).

## Licensing

No project license has been selected. Public visibility alone does not grant an
open-source license. Third-party data, models, and dependencies retain their own
terms.
