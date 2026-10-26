# TCGA / GTEx explorer input

TCGA TPM tables and metadata are available locally. The added 54 GTEx v11
`gene_reads_adult_gtex_v11_*.gct` files contain raw counts, not TPM.
The matching 54 individual-sample GTEx v11 TPM GCT files have also been supplied locally.
Use the importer below to create the explorer tables; raw counts remain excluded.
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

## GTEx v11 import and anatomical grouping

Download **Gene TPMs by tissue**, `gene_tpm_adult_gtex_v11_*.gct.gz`, from
https://www.gtexportal.org/home/downloads/adult-gtex/overview (not gene reads,
median TPMs, transcript TPMs or LCM). Compressed or extracted files are supported:

```sh
python scripts/prepare_gtex.py tcga_read_counts tcga_read_counts
```

The converter streams rows, strips ENSG versions, preserves numeric TPM and
original sample IDs, and combines sample columns according to
`docs/gtex_tissue_groups.tsv`. It checks gene alignment, sample uniqueness,
dimensions and numeric values. Missing tissues are reported; existing output
tables are not overwritten. Restart Shiny after import. Each output has a
`.samples.tsv` provenance file, excluded from expression-file discovery.

Grouping pools samples; it does not average values or perform statistical
clustering. Brain regions, skin sites, kidney regions, colon sites and other
sub-tissues share organ groups. Sub-tissue identities remain in the provenance
files. Samples from one donor across sub-tissues are not independent donors.
Ileum stays separate from colon; salivary gland stays separate from esophagus;
spleen, whole blood and the two cultured-cell types stay separate.

The historical converter is `get_gtex_tcga_counts.R` (local original retained).
Its `grouped start` section contains the old pooling rules; the original
`tcga_gtex_tissues.tsv` supplies TCGA anatomical facets. The bootstrap now uses
the explicit GTEx mapping above for GTEx facet labels. Legacy TCGA mappings
remain, including coarse assignments such as MESO to Lung and UVM to Skin;
these are display groups, not validated matched normal controls. No appropriate
GTEx counterpart is assumed for every cancer type.

Do not run the old converter for normalization: it estimates lengths from
Ensembl gene start/end coordinates, rather than the release's quantification
annotation, and its v10 filename handling does not match the supplied v11 data.
Official TPM avoids this approximation. TPM units alone do not remove TCGA–GTEx
processing/batch differences; these plots are descriptive comparisons.
