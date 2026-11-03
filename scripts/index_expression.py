"""Build byte-offset indexes for the Shiny expression explorer without changing TPM files."""
import argparse
import os
from pathlib import Path
import re
import tempfile


def build_index(path):
    path = Path(path)
    before = path.stat()
    target = Path(str(path) + '.gidx')
    temporary = None
    try:
        with path.open('rb') as source, tempfile.NamedTemporaryFile(
                mode='w', encoding='utf-8', newline='', dir=path.parent, delete=False) as out:
            temporary = Path(out.name)
            out.write(f'{before.st_size}\t{before.st_mtime:.6f}\n')
            out.write('id\toffset\tlength\n')
            header = source.readline()
            if b'\t' not in header:
                raise ValueError(f'Not a sample expression table: {path}')
            offset = len(header)
            for line in source:
                if not line.strip():
                    raise ValueError(f'Blank row in {path}')
                gene = line.split(b'\t', 1)[0].decode('utf-8')
                if path.name.endswith('_tpm.tsv'):
                    gene = re.split(r'[._]', gene)[0]
                out.write(f'{gene}\t{offset}\t{len(line)}\n')
                offset += len(line)
        after = path.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise RuntimeError(f'Source changed while indexing: {path}')
        os.replace(temporary, target)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()
    return target


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder', type=Path)
    args = parser.parse_args()
    paths = sorted(p for p in args.folder.glob('*.tsv') if
                   p.name.endswith('_tpm.tsv') or
                   (p.name.upper().startswith('GTEX') and not p.name.endswith('.samples.tsv')))
    if not paths:
        parser.error('No expression tables found')
    for path in paths:
        print(build_index(path), flush=True)
