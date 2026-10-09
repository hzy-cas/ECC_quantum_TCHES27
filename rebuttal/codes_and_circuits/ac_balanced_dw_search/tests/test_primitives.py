import json
import os
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "build" / "ac_balanced_estimator"


def discover_data_root() -> Path:
    override = os.environ.get("ECC_SHOR_REPOSITORY_ROOT")
    if override:
        data_root = (
            Path(override).expanduser().resolve() / "mul_line" / "C++" / "data"
        )
        if data_root.is_dir():
            return data_root
        raise RuntimeError(
            "ECC_SHOR_REPOSITORY_ROOT does not contain mul_line/C++/data"
        )
    for candidate in (ROOT, *ROOT.parents):
        data_root = candidate / "mul_line" / "C++" / "data"
        if data_root.is_dir():
            return data_root
    raise RuntimeError(
        "could not locate mul_line/C++/data; set ECC_SHOR_REPOSITORY_ROOT"
    )


DATA = discover_data_root()


def run_json(*arguments: str) -> dict[str, int]:
    completed = subprocess.run(
        [str(BINARY), *arguments, "--data-root", str(DATA)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return json.loads(completed.stdout)


class PrimitiveRegressionTests(unittest.TestCase):
    def test_ac_multiplication_matches_table3_gate_stream(self):
        expected = {
            163: (2718, 906, 180537, 711),
            233: (4023, 1341, 377945, 935),
            283: (5004, 1668, 568948, 1119),
            571: (10707, 3569, 2490565, 2205),
        }
        for n, (qubits, toffoli, cnot, depth) in expected.items():
            with self.subTest(n=n):
                stats = run_json("primitive", "--n", str(n), "--kind", "multiplication")
                self.assertEqual(stats["qubits"], qubits)
                self.assertEqual(stats["toffoli"], toffoli)
                self.assertEqual(stats["cnot"], cnot)
                self.assertEqual(stats["full_depth"], depth)
                self.assertEqual(stats["current_depth"], depth)
                self.assertEqual(stats["toffoli_depth"], 1)

    def test_ac_inversion_chain_regression(self):
        stats = run_json("primitive", "--n", "163", "--kind", "inversion")
        self.assertEqual(stats["toffoli"], 8154)
        self.assertEqual(stats["cnot"], 1678371)
        self.assertEqual(stats["full_depth"], 6058)
        self.assertEqual(stats["toffoli_depth"], 8)

    def test_clean_point_add_k1(self):
        controlled = run_json(
            "layer", "--n", "163", "--lanes", "1", "--mode", "accumulation"
        )
        uncontrolled = run_json(
            "layer", "--n", "163", "--lanes", "1", "--mode", "reduction"
        )
        self.assertEqual(controlled["toffoli"], 20258)
        self.assertEqual(uncontrolled["toffoli"], 19932)
        self.assertEqual(controlled["toffoli"] - uncontrolled["toffoli"], 326)
        self.assertEqual(controlled["toffoli_depth"], 21)
        self.assertEqual(uncontrolled["toffoli_depth"], 20)


if __name__ == "__main__":
    unittest.main()
