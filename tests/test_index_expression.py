import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('index_expression', Path(__file__).resolve().parents[1] / 'scripts/index_expression.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class IndexTests(unittest.TestCase):
    def test_offsets_preserve_rows_and_duplicate_ids(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'TCGA-test_tpm.tsv'
            for newline in (b'\n', b'\r\n'):
                rows = [b'ENSG1.1_A_coding\t1.5', b'ENSG1.2_B_coding\t2']
                path.write_bytes(newline.join([b'id\tsample', *rows]) + newline)
                index = module.build_index(path).read_text().splitlines()[2:]
                with path.open('rb') as handle:
                    for entry, row in zip(index, rows):
                        gene, offset, size = entry.split('\t')
                        self.assertEqual(gene, 'ENSG1')
                        handle.seek(int(offset))
                        self.assertEqual(handle.read(int(size)), row + newline)


if __name__ == '__main__':
    unittest.main()
