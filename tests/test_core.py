"""Lightweight checks for the retained backend; no model downloads required."""

import ast
import builtins
import contextlib
import io
import math
from pathlib import Path
import symtable
import sys
from types import SimpleNamespace
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import Epitope
import TAP_calculation


def load_function(filename, name):
    """Isolate pure helpers from modules that eagerly load ML dependencies."""
    tree = ast.parse((ROOT / filename).read_text(encoding='utf-8'))
    fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name)
    scope = {}
    exec(compile(ast.Module(body=[fn], type_ignores=[]), filename, 'exec'), scope)
    return scope[name]


class CoreTests(unittest.TestCase):
    def test_cleavage_cache_and_fasta_count_guard(self):
        tree = ast.parse((ROOT / 'main_final.py').read_text())
        fn = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == 'get_pepsickle')
        reads = []
        def read_table(path):
            reads.append(path)
            return {'protein_id': SimpleNamespace(max=lambda: 2 if path == 'default' else 3)}
        scope = dict(pepsickle_filese=['default', 'default', 'mismatch'],
                     loaded_pepsickle={}, proteins=[None, None],
                     pandas=SimpleNamespace(read_table=read_table))
        exec(compile(ast.Module(body=[fn], type_ignores=[]), 'main_final.py', 'exec'), scope)
        load = scope['get_pepsickle']
        self.assertIs(load(0), load(1))
        self.assertEqual(reads, ['default'])
        with self.assertRaisesRegex(ValueError, 'do not match'):
            load(2)
        self.assertNotIn('mismatch', scope['loaded_pepsickle'])

    def test_fixed_positions(self):
        fixed = load_function('main_final.py', 'get_fixed_pos')
        self.assertEqual(fixed('SLLMWITQC'), ('SLLMWITQC', []))
        self.assertEqual(fixed('sllmwitqc'), ('SLLMWITQC', []))
        self.assertEqual(fixed('sLlMwitqc'), ('SLLMWITQC', [2, 4]))

    def test_fixed_position_mismatches(self):
        mismatch = load_function('Search_Uniprot.py', 'missmatches_fix')
        self.assertEqual(mismatch('ACDEFGHIK', 'ACDEYGHIK', [1]), [4])
        self.assertEqual(mismatch('ACDEFGHIK', 'ACDEFGHIK', [5]), [])
        self.assertGreater(len(mismatch('ACDEFGHIK', 'ACDEYGHIK', [5])), 4)

    def test_kmers(self):
        kmers = load_function('MHCflurry.py', 'get_all_kmers')
        self.assertEqual(kmers('ACDEFG', 3), ['ACD', 'CDE', 'DEF', 'EFG'])
        self.assertEqual(kmers('AC', 3), [])

    def test_tap_scores(self):
        # Independent expected values for a nine-alanine peptide.
        hit = ['AAAAAAAAA']
        with contextlib.redirect_stdout(io.StringIO()):
            TAP_calculation.calculate([SimpleNamespace(protein_hits=[hit])], 0.2)
        self.assertAlmostEqual(hit[1], 0.55 + 0.2 * (-1.56 - 0.25 - 0.10))
        self.assertAlmostEqual(hit[2], 10 ** hit[1])
        self.assertAlmostEqual(hit[3], 0.5 * .53 + .25 * .58 + .2 * .45 + .25 * .46 + .43)

    def test_q_scores_and_zero_denominators(self):
        for mhc, tap10 in [(100.0, 2.0), (0.0, 0.0)]:
            hit = [0.0] * 18
            hit[6], hit[11], hit[15], hit[16], hit[17] = .5, mhc, 1.0, tap10, .75
            with contextlib.redirect_stdout(io.StringIO()):
                Epitope.calculate_q_score_new2([SimpleNamespace(protein_hits=[hit])])
            self.assertEqual(len(hit), 21)
            self.assertAlmostEqual(hit[18], .5 / ((tap10 or 1e-6) * (mhc or 1e-6)))
            self.assertAlmostEqual(hit[19], .5 / (math.log10(9.75) * (mhc or 1e-6)))
            self.assertAlmostEqual(hit[20], .5 * .75 / (mhc or 1e-6))

    def test_sources_compile_and_have_resolved_globals(self):
        paths = list(ROOT.glob('*.py')) + list((ROOT / 'model').glob('*.py'))
        for path in paths:
            if path.name.startswith('_'):
                continue
            with self.subTest(path=path.name):
                source = path.read_text(encoding='utf-8')
                compile(source, str(path), 'exec')
                table = symtable.symtable(source, str(path), 'exec')
                known = set(table.get_identifiers()) | set(dir(builtins))
                pending = list(table.get_children())
                while pending:
                    child = pending.pop()
                    pending.extend(child.get_children())
                    missing = {s.get_name() for s in child.get_symbols()
                               if s.is_global() and s.is_referenced() and s.get_name() not in known}
                    self.assertFalse(missing, f'{path.name}: {missing}')


if __name__ == '__main__':
    unittest.main()
