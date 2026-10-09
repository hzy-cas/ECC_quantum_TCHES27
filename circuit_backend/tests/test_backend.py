"""Behavioral checks for the shared Python and C++ NCT backends."""
import importlib.util
import itertools
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from circuit_backend.python.qasm import FlatGateManager
from circuit_backend.python.stats_utils import get_exact_resources_optimized


def build(sequence, chunk_size=4):
    gm = FlatGateManager()
    gm.CHUNK_LIMIT = chunk_size
    for op, *args in sequence:
        getattr(gm, {'X': 'add_X', 'C': 'add_CNOT', 'T': 'add_Toffoli'}[op])(*args)
    return gm


def metrics(gm, width=6):
    counts, fd, cd, td = get_exact_resources_optimized(gm, width)
    return [fd, cd, td, counts['Toffoli_count'], counts['CNOT_count']]


def truth(sequence, state):
    state = list(state)
    for op, *a in sequence:
        if op == 'X': state[a[0]] ^= 1
        elif op == 'C': state[a[1]] ^= state[a[0]]
        else: state[a[2]] ^= state[a[0]] & state[a[1]]
    return state


class PythonBackendTests(unittest.TestCase):
    def test_swapped_controls_cancel_in_imported_stream(self):
        gm = build([])
        # Bypass add_Toffoli's normalization to exercise imported old gate streams.
        for a, b in [(0, 1), (1, 0)]:
            gm._append((gm.OP_TOFFOLI << gm.SHIFT_OP) | (a << gm.SHIFT_ARG) | b)
            gm._append(2)
        self.assertEqual(metrics(gm), [0, 0, 0, 0, 0])

    def test_both_cross_chunk_cancellation_cases(self):
        for prefix in [1, 3]:
            with self.subTest(prefix=prefix):
                gm = build([('X', 3)] * prefix + [('T', 0, 1, 2), ('T', 1, 0, 2)])
                self.assertEqual(metrics(gm), [1, 1, 0, 0, 0])
                self.assertEqual(gm.get_stats(), (0, 0, 1))

    def test_intervening_gate_blocks_cancellation(self):
        sequence = [('T', 0, 1, 2), ('X', 0), ('T', 1, 0, 2)]
        gm = build(sequence)
        self.assertEqual(metrics(gm)[3], 2)
        for state in itertools.product(range(2), repeat=3):
            actual = list(state)
            gm.simulate(actual)
            self.assertEqual(actual, truth(sequence, state))

    def test_reverse_across_chunks_restores_every_basis_state(self):
        sequence = [('X', 3), ('T', 1, 0, 2), ('C', 2, 3), ('T', 3, 2, 1)]
        gm = build(sequence)
        end = gm.current_pointer()
        gm.replay_reverse(0, end)
        for state in itertools.product(range(2), repeat=4):
            actual = list(state)
            gm.simulate(actual)
            self.assertEqual(actual, list(state))
        self.assertEqual(metrics(gm), [0, 0, 0, 0, 0])

    def test_optimization_preserves_truth_tables(self):
        rng = random.Random(163)
        for _ in range(30):
            sequence = []
            for _ in range(30):
                op = rng.choice('XCT')
                args = rng.sample(range(4), {'X': 1, 'C': 2, 'T': 3}[op])
                sequence.append(tuple([op] + args))
            gm = build(sequence)
            metrics(gm, 4)
            for state in itertools.product(range(2), repeat=4):
                actual = list(state)
                gm.simulate(actual)
                self.assertEqual(actual, truth(sequence, state))

    def test_old_entry_points_share_the_same_objects(self):
        for i, relative in enumerate([
            'jsb25/init/qasm.py',
            'mul_line/algorithm1_inplace_point_add/backend/qasm.py',
        ]):
            spec = importlib.util.spec_from_file_location(f'backend_alias_{i}', ROOT / relative)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            self.assertIs(module.FlatGateManager, FlatGateManager)


class CrossLanguageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler = shutil.which('c++') or shutil.which('g++')
        if not compiler:
            raise unittest.SkipTest('a C++17 compiler is required')
        cls.directory = tempfile.TemporaryDirectory(prefix='nct-backend-test-')
        cls.binary = Path(cls.directory.name) / 'driver'
        subprocess.run([compiler, '-std=c++17', '-O1',
                        str(ROOT / 'circuit_backend/tests/driver.cpp'),
                        str(ROOT / 'circuit_backend/cpp/GateManager.cpp'),
                        '-o', str(cls.binary)], check=True)

    @classmethod
    def tearDownClass(cls):
        cls.directory.cleanup()

    def test_identical_metrics_with_different_chunk_boundaries(self):
        rng = random.Random(20261005)
        cases = [
            [('T', 0, 1, 2), ('T', 1, 0, 2)],
            [('X', 3)] * 3 + [('T', 0, 1, 2)] * 2,
            [('X', 3)] + [('T', 0, 1, 2)] * 2,
        ]
        for _ in range(300):
            seq = []
            for _ in range(rng.randint(5, 50)):
                if seq and rng.random() < .2:
                    gate = seq[-1]
                    if gate[0] == 'T': gate = ('T', gate[2], gate[1], gate[3])
                    seq.append(gate)
                else:
                    op = rng.choice('XCT')
                    seq.append(tuple([op] + rng.sample(range(6), {'X':1,'C':2,'T':3}[op])))
            cases.append(seq)
        lines = []
        for seq in cases:
            lines.append(f'6 {len(seq)}')
            for op, *args in seq:
                lines.append(op + ' ' + ' '.join(map(str, args + [0] * (3-len(args)))))
        expected = [metrics(build(seq, 4)) for seq in cases]
        for size in [4, 5, 10000000]:
            with self.subTest(cpp_chunk_size=size):
                result = subprocess.run([str(self.binary), str(size)],
                                        input='\n'.join(lines)+'\n', text=True,
                                        capture_output=True, check=True)
                actual = [list(map(int, line.split())) for line in result.stdout.splitlines()]
                self.assertEqual(actual, expected)


if __name__ == '__main__':
    unittest.main()
