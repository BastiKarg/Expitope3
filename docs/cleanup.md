# Source cleanup

## Scope and entry point

The retained root is `main_final.py` → `run_waiting_loaded_multicore()`.
`main_final_tac.py` was identical to `main_final.py`. The worker loops forever;
the positional-argument dispatch following that call was unreachable during
normal execution, and the later benchmark code additionally followed `exit(0)`.
Those blocks and obsolete entry points were removed.

The traversal follows local function references, imported aliases, the
`mlp_predict.MLP` re-export, and multiprocessing initializer/callback references.
It retains all job branches: normal search, indel/multigap, k-mers, and custom
protein sequences. Lightning model methods remain intact because the framework
and checkpoint loading invoke methods indirectly. Unreferenced ordinary Epitope
setters were removed.

## Retained flow

- The worker loads models and dispatches jobs.
- Search callbacks call `Search_Uniprot` and construct `Epitope` objects.
- `Run_Pepsickle`, `MHCflurry`, `TAP_calculation`, and `Epitope` append scores.
- `embedding` computes ProtT5 features; `mlp_predict` classifies and exports them.
- `Get_Expression` joins tissue tables; `Logsums_and_ICR_calculation` computes
  summaries; `Output_File` writes final tables.

## Removed material

The cleanup reduces the supplied 30 Python files (16,674 lines) to 12 backend
modules (approximately 1,700 lines), before adding lightweight tests. Removed
material includes duplicate/backup main files, the obsolete command-line script,
old epitope/embedding implementations, unreferenced model and dataset variants,
unused NetChop and alternative Pepsickle wrappers, experimental scripts,
unreachable functions, commented experiments, and inert code strings.

Four excluded experiments were already syntactically invalid: root `__init__.py`,
`all_hla_code.py`, `mlp_prediction.py`, and `tmp_testing.py`.

Original Python sources are preserved in the local, Git-ignored
`.cleanup-backup/` directory. They are intentionally absent from public Git
history. The original environment exports, R application, data, and checkpoints
remain local and are ignored by Git.

## Small operational corrections

- Set `TF_USE_LEGACY_KERAS` before imports that can load TensorFlow.
- Create the queue and completion-marker directories before polling.
- Restore the missing `pandas` binding used by the retained v4 expression branch.
- Remove unused model initialization and overwritten search-index assignments.
- Remove the normal-search alternatives guarded by constant `parallel=True` and
  `new=False`, including their unused serial-search helper.

The effective hard-coded normal-search offsets are unchanged. Changing those
offsets, the multiprocessing design, or score layouts would require a separate
behavioral review with reference results.

## Validation

Six dependency-light unit tests cover fixed positions, fixed-position mismatch
rejection, k-mers, TAP arithmetic, Q-score zero-denominator handling, compilation,
and unresolved global names. GitHub Actions runs them with Python 3.9.

A local AST comparison against the original backup checked all 36 retained
free functions: 33 have identical executable syntax after excluding imports and
inert strings. The three differences are the worker dispatcher, normal-search
orchestration, and k-mer orchestration described above. This is structural
evidence, not proof of end-to-end numerical equivalence.

The full ML stack was not installed and complete inference was not executed.
The requirements file reconstructs historical versions and still needs an
installation/inference check on Linux with the project assets.

## Inherited issues requiring follow-up

- Normal-search database offsets are hard-coded; k-mer offsets are calculated
  from chunks. Database ordering must remain consistent with cleavage IDs.
- Workers depend on fork-inherited global state. Queue polling can hang after
  worker exceptions; retries do not guarantee recovery.
- Custom-sequence search returns a string on no hits, which downstream routines
  do not consistently handle.
- Input parsing, peptide length constraints, malformed queue files, and output
  directory names need validation before multi-user deployment.
- Cleavage-model fallback catches broad exceptions, and some prediction paths
  emit extensive debug output.
- Data provenance, model distribution, reproducible dependency locking, and a
  project license remain to be supplied by the project maintainers.
