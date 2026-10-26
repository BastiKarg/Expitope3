"""Stream tissue-specific GTEx v11 TPM GCTs into explorer tables (stdlib only)."""
import argparse
import csv
import gzip
import math
from contextlib import ExitStack
from itertools import zip_longest
from pathlib import Path
import tempfile


def convert(paths, target, tissues):
    """Concatenate samples, never average tissues or normalize raw counts."""
    with ExitStack() as stack:
        readers, samples, provenance, dimensions = [], [], [], []
        for path, tissue in zip(paths, tissues):
            if not path.name.startswith('gene_tpm_adult_gtex_v11_'):
                raise ValueError('Only official gene_tpm_adult_gtex_v11 files are accepted')
            opener = gzip.open if path.suffix == '.gz' else open
            handle = stack.enter_context(opener(path, 'rt', newline=''))
            if handle.readline().strip() != '#1.2':
                raise ValueError(f'Unsupported GCT version: {path}')
            nrows, ncols = map(int, handle.readline().split())
            reader = csv.reader(handle, delimiter='\t')
            header = next(reader)
            if header[:2] != ['Name', 'Description'] or len(header) != ncols + 2:
                raise ValueError(f'Invalid GCT header: {path}')
            samples.extend(header[2:])
            provenance.extend((sample, tissue, path.name) for sample in header[2:])
            readers.append(reader)
            dimensions.append((nrows, ncols))
        if len(samples) != len(set(samples)):
            raise ValueError('Duplicate sample IDs within pooled tissue')
        target.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=target.parent) as temporary:
            staged = Path(temporary) / target.name
            seen = set()
            with staged.open('w', newline='') as out:
                writer = csv.writer(out, delimiter='\t', lineterminator='\n')
                writer.writerow(['id', *samples])
                for rows in zip_longest(*readers):
                    if any(row is None for row in rows):
                        raise ValueError('GCT gene counts differ')
                    ids = [row[0].split('.')[0] for row in rows]
                    if len(set(ids)) != 1 or ids[0] in seen:
                        raise ValueError('Gene order differs or normalized gene IDs are duplicated')
                    seen.add(ids[0])
                    values = []
                    for row, (_, ncols) in zip(rows, dimensions):
                        if len(row) != ncols + 2:
                            raise ValueError('Invalid GCT row width')
                        if any(not math.isfinite(float(v)) or float(v) < 0 for v in row[2:]):
                            raise ValueError('TPM must be finite and nonnegative')
                        values.extend(row[2:])
                    writer.writerow([ids[0], *values])
            if any(nrows != len(seen) for nrows, _ in dimensions):
                raise ValueError('GCT declared gene count differs from actual rows')
            staged.replace(target)
        with target.with_suffix('.samples.tsv').open('w', newline='') as out:
            writer = csv.writer(out, delimiter='\t', lineterminator='\n')
            writer.writerow(['sample_id', 'source_tissue', 'source_file'])
            writer.writerows(provenance)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--mapping', type=Path, default=Path(__file__).resolve().parents[1] / 'docs/gtex_tissue_groups.tsv')
    args = parser.parse_args()
    groups = {}
    with args.mapping.open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            stem = 'gene_tpm_adult_gtex_v11_' + row['source_tissue'] + '.gct'
            candidates = [args.input / stem, args.input / (stem + '.gz')]
            found = [p for p in candidates if p.exists()]
            if found:
                groups.setdefault(row['group'], []).append((found[0], row['source_tissue']))
            else:
                print('Missing TPM:', row['source_tissue'])
    if not groups:
        parser.error('No GTEx v11 tissue TPM GCT files found; raw gene_reads files cannot be substituted')
    for group, inputs in groups.items():
        target = args.output / (group + '.tsv')
        if target.exists():
            raise FileExistsError(f'Refusing to overwrite {target}')
        convert([p for p, _ in inputs], target, [t for _, t in inputs])
        print('Created', target)


if __name__ == '__main__':
    main()
