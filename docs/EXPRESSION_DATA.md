# TCGA / GTEx explorer input

TCGA TPM tables and metadata are available locally. The missing part is
**individual-sample GTEx expression**, grouped into one file per tissue.
The existing reference database tables contain tissue summaries and cannot
reconstruct sample distributions for boxplots, violin plots or beeswarms.

Place files alongside the existing TPM files in `tcga_read_counts/`, or choose
a directory with `EXPITOPE_EXPRESSION_DIR` before starting Shiny.
Example filename: `GTEX_liver.tsv` (match the `tissue` column in
`tcga_gtex_tissues.tsv` exactly for anatomical grouping).

Tab-separated, header included:

```text
id	GTEX_sample_1	GTEX_sample_2	GTEX_sample_3
ENSG00000000003	12.4	9.8	0
ENSG00000000005	0	0.2	1.1
```

- One unique, unversioned ENSG ID per row, in a column named `id`.
- One actual sample per remaining column; retain unique sample identifiers.
- Numeric, finite, non-negative TPM values, before log transformation.
- Do not mix TPM with raw counts or tissue averages. The UI applies
  `log2(expression + 1)` only when selected.
- Provide all genes needed for comparisons; a missing gene produces an error
  rather than fabricated zero expression.
- `tcga_gtex_tissues.tsv` supplies `tissue`, `project`, `body` mappings.
  The existing TCGA-derived gene names are reused when standalone
  `ensg_genename.tsv` / `ensg_list` are absent.

The importer includes `GTEX*.tsv` alongside `*_tpm.tsv`; raw TCGA count files in
that directory are excluded. Files ending `_tpm.tsv` use the existing annotated
TCGA first-column format, so use plain `.tsv` for the GTEx format above.
If you have an original GTEx GCT matrix plus a sample-to-tissue annotation file,
those can be converted into this layout without inventing sample values.
