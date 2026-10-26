"""Small conversion regression: pooling, provenance, and invalid input rejection."""
import importlib.util
from pathlib import Path
import tempfile

spec = importlib.util.spec_from_file_location('prepare', Path(__file__).resolve().parents[1] / 'scripts/prepare_gtex.py')
prepare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare)
with tempfile.TemporaryDirectory() as folder:
    root = Path(folder)
    paths = [root / f'gene_tpm_adult_gtex_v11_skin_{i}.gct' for i in range(2)]
    for i, path in enumerate(paths):
        path.write_text(f'#1.2\n2 1\nName\tDescription\tsample{i}\nENSG1.2\tA\t{i}.5\nENSG2.1\tB\t0\n')
    out = root / 'GTEX_skin.tsv'
    prepare.convert(paths, out, ['skin_a', 'skin_b'])
    assert out.read_text().splitlines() == ['id\tsample0\tsample1', 'ENSG1\t0.5\t1.5', 'ENSG2\t0\t0']
    assert 'sample1\tskin_b' in out.with_suffix('.samples.tsv').read_text()
    original = out.read_bytes()
    paths[1].write_text(paths[1].read_text().replace('ENSG2.1', 'ENSG3.1'))
    try:
        prepare.convert(paths, out, ['a', 'b'])
    except ValueError:
        pass
    else:
        raise AssertionError('Mismatched genes accepted')
    assert out.read_bytes() == original
    try:
        prepare.convert([root / 'gene_reads_adult_gtex_v11_liver.gct'], out, ['liver'])
    except ValueError:
        pass
    else:
        raise AssertionError('Raw counts accepted')
print('GTEx conversion tests passed')
