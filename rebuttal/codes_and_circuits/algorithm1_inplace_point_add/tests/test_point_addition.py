"""Interface, functional, and resource regressions for the streamlined circuit."""

from __future__ import annotations

import inspect
from pathlib import Path
import sys
import unittest


HERE = Path(__file__).resolve().parent
PACKAGE_ROOT = HERE.parent
PACKAGE_PARENT = PACKAGE_ROOT.parent
if not __package__ and str(PACKAGE_PARENT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_PARENT))

if __package__:
    from ..backend.qasm import FlatGateManager, MatrixCache
    from ..circuit import circuit_library
    from ..circuit import optimal_inversion
    from ..circuit.config import (
        BINARY_ECC,
        CONFIGS,
        INVERSION_MODES,
        OPTIMAL_DEPTH,
        SUPPORTED_SIZES,
    )
    from ..circuit.inversion import Inverison_Itoh_Tsujii_based
    from ..circuit.point_addition import Point_addition
    from ..resources.point_addition_resources import (
        estimate_one,
        raw_toffoli_closed_form,
    )
    from ..verification.verify_point_addition import (
        algorithm1_classical,
        deterministic_fixture,
        simulate_packed,
        verify_one,
    )
else:
    from algorithm1_inplace_point_add.backend.qasm import (
        FlatGateManager,
        MatrixCache,
    )
    from algorithm1_inplace_point_add.circuit import (
        circuit_library,
        optimal_inversion,
    )
    from algorithm1_inplace_point_add.circuit.config import (
        BINARY_ECC,
        CONFIGS,
        INVERSION_MODES,
        OPTIMAL_DEPTH,
        SUPPORTED_SIZES,
    )
    from algorithm1_inplace_point_add.circuit.inversion import (
        Inverison_Itoh_Tsujii_based,
    )
    from algorithm1_inplace_point_add.circuit.point_addition import Point_addition
    from algorithm1_inplace_point_add.resources.point_addition_resources import (
        estimate_one,
        raw_toffoli_closed_form,
    )
    from algorithm1_inplace_point_add.verification.verify_point_addition import (
        algorithm1_classical,
        deterministic_fixture,
        simulate_packed,
        verify_one,
    )


class LocalComponentTests(unittest.TestCase):

    def test_local_backend_supports_stats_and_reverse_replay(self):
        gm = FlatGateManager()
        gm.add_CNOT(0, 1)
        forward_end = gm.current_pointer()
        gm.add_Toffoli(0, 1, 2)
        gm.replay_reverse(0, forward_end)
        self.assertEqual(gm.get_stats(), (1, 2))
        self.assertEqual(MatrixCache().get_matrix_data("/missing", 1), ([], []))

    def test_arithmetic_functions_keep_mul_line_api(self):
        expected_parameters = {
            "CNOT_n": ("gm", "a", "b", "n"),
            "Toffoli_gate": ("gm", "a", "b", "c"),
            "mul_matrix": ("gm", "matrix_loader", "x", "n", "file_path"),
            "mul": ("gm", "matrix_loader", "a0", "b0", "c0", "n", "path_config"),
            "Square": ("gm", "matrix_loader", "x", "n", "power", "base_path"),
        }
        for name, parameters in expected_parameters.items():
            with self.subTest(function=name):
                function = getattr(circuit_library, name)
                self.assertTrue(callable(function))
                self.assertEqual(tuple(inspect.signature(function).parameters), parameters)

    def test_only_two_new_public_circuit_functions(self):
        self.assertEqual(Point_addition.__name__, "Point_addition")
        self.assertEqual(
            Inverison_Itoh_Tsujii_based.__name__,
            "Inverison_Itoh_Tsujii_based",
        )

    def test_optimal_inversions_keep_mul_line_api(self):
        for name in (
            "Inversion_163",
            "Inversion_233",
            "Inversion_283",
            "Inversion_571",
        ):
            with self.subTest(function=name):
                function = getattr(optimal_inversion, name)
                self.assertTrue(callable(function))
                self.assertEqual(function.__name__, name)


class AlgebraAndLayoutTests(unittest.TestCase):

    def test_algorithm1_q0_q1_for_all_fields(self):
        for n in SUPPORTED_SIZES:
            with self.subTest(n=n):
                fixture = deterministic_fixture(n)
                self.assertEqual(
                    algorithm1_classical(
                        0,
                        fixture.x1,
                        fixture.y1,
                        fixture.x2,
                        fixture.y2,
                        fixture.curve_a,
                        n,
                    ),
                    (fixture.x1, fixture.y1, 0),
                )
                self.assertEqual(
                    algorithm1_classical(
                        1,
                        fixture.x1,
                        fixture.y1,
                        fixture.x2,
                        fixture.y2,
                        fixture.curve_a,
                        n,
                    ),
                    (fixture.x3, fixture.y3, 0),
                )

    def test_width_and_toffoli_closed_forms(self):
        expected_widths = {
            BINARY_ECC: {163: 11362, 233: 18133, 283: 24202, 571: 58818},
            OPTIMAL_DEPTH: {163: 13174, 233: 20815, 283: 27538, 571: 65956},
        }
        expected_toffoli = {163: 40190, 233: 64834, 283: 87302, 571: 215282}
        for n in SUPPORTED_SIZES:
            for inversion_mode in INVERSION_MODES:
                with self.subTest(n=n, inversion_mode=inversion_mode):
                    cfg = CONFIGS[n]
                    extra_workspace = (
                        2 * n + 2 * cfg["block_size"]
                        if inversion_mode == OPTIMAL_DEPTH
                        else 0
                    )
                    width = (
                        1
                        + 5 * n
                        + 2 * cfg["block_size"]
                        + (cfg["inversion_multiplications"] + 1)
                        * cfg["multiplication_targets"]
                        + extra_workspace
                    )
                    self.assertEqual(
                        width, expected_widths[inversion_mode][n]
                    )
                    self.assertEqual(
                        raw_toffoli_closed_form(n), expected_toffoli[n]
                    )


class OptimizedFrobeniusTests(unittest.TestCase):

    def test_saved_non_inplace_cnot_streams_match_their_matrices(self):
        required = {
            (163, 34),
            (233, 40),
            (233, 104),
            (283, 10),
            (283, 26),
            (571, 10),
            (571, 26),
        }
        files = sorted(
            (PACKAGE_ROOT / "square").glob(
                "RES_n*/result_square_Matrix_2_*.txt"
            )
        )
        self.assertTrue(files)
        found = {
            (
                int(path.parent.name.removeprefix("RES_n")),
                int(path.stem.rsplit("_", 1)[1]),
            )
            for path in files
        }
        self.assertTrue(required.issubset(found))
        for path in files:
            n = int(path.parent.name.removeprefix("RES_n"))
            dimension = 2 * n
            with self.subTest(n=n, filename=path.name):
                # The gate manager and matrix cache use the project's local backend;
                # the squaring emitter itself calls circuit_library.mul_matrix directly.
                gm = FlatGateManager()
                matrix_loader = MatrixCache()
                logical_outputs = circuit_library.mul_matrix(
                    gm,
                    matrix_loader,
                    list(range(dimension)),
                    dimension,
                    str(path),
                )
                _, matrix_rows = matrix_loader.get_matrix_data(
                    str(path), dimension
                )
                self.assertEqual(len(matrix_rows), dimension)
                self.assertEqual(len(logical_outputs), dimension)
                self.assertEqual(len(set(logical_outputs)), dimension)
                for index in range(n):
                    self.assertEqual(matrix_rows[index], [index])
                    self.assertEqual(
                        [bit for bit in matrix_rows[n + index] if bit >= n],
                        [n + index],
                    )

                case_count = 4
                initial = [
                    ((wire + 1) * 0x9E3779B1) & ((1 << case_count) - 1)
                    for wire in range(dimension)
                ]
                final = simulate_packed(gm, initial, case_count)
                for output_index, row in enumerate(matrix_rows):
                    expected = 0
                    for input_index in row:
                        expected ^= initial[input_index]
                    self.assertEqual(
                        final[logical_outputs[output_index]], expected
                    )
                self.assertEqual(gm.get_stats()[0], 0)

    def test_binary_ecc_inversion_has_no_matrix_power_fallback(self):
        source = inspect.getsource(Inverison_Itoh_Tsujii_based)
        self.assertNotIn("square_Matrix_2_1", source)
        self.assertNotIn("_FROBENIUS_ROWS", source)

        required_powers = {
            163: (1, 2, 4, 8, 16, 32, 34, 64),
            233: (1, 2, 4, 8, 16, 32, 40, 64, 104),
            283: (1, 2, 4, 8, 10, 16, 26, 32, 64, 128),
            571: (1, 2, 4, 8, 10, 16, 26, 32, 58, 64, 128, 256),
        }
        for n, powers in required_powers.items():
            missing = [
                power
                for power in powers
                if not (
                    PACKAGE_ROOT
                    / "square"
                    / f"RES_n{n}"
                    / f"result_square_Matrix_2_{power}.txt"
                ).is_file()
            ]
            if n == 571:
                self.assertIn(missing, ([], [58]))
            else:
                self.assertEqual(missing, [])


class FullGateStreamTests(unittest.TestCase):

    def test_q0_q1_and_all_ancillas_clean_for_n163(self):
        expected = {
            BINARY_ECC: (11362, 8517318),
            OPTIMAL_DEPTH: (13174, 8520802),
        }
        for inversion_mode in INVERSION_MODES:
            with self.subTest(inversion_mode=inversion_mode):
                row = verify_one(163, inversion_mode)
                self.assertTrue(row["algebra_passed"])
                self.assertTrue(row["gate_stream_passed"])
                self.assertEqual(row["dirty_ancilla_count"], 0)
                self.assertEqual(row["width"], expected[inversion_mode][0])
                self.assertEqual(row["raw_toffoli"], 40190)
                self.assertEqual(row["raw_cnot"], expected[inversion_mode][1])

    def test_resource_regression_for_n163(self):
        expected = {
            BINARY_ECC: (11362, 7375426, 32656, 32656, 46),
            OPTIMAL_DEPTH: (13174, 7065450, 29872, 30052, 42),
        }
        for inversion_mode in INVERSION_MODES:
            with self.subTest(inversion_mode=inversion_mode):
                row = estimate_one(163, inversion_mode)
                values = expected[inversion_mode]
                self.assertEqual(row["qubits"], values[0])
                self.assertEqual(row["toffoli"], 40190)
                self.assertEqual(row["cnot"], values[1])
                self.assertEqual(row["full_depth"], values[2])
                self.assertEqual(row["current_depth"], values[3])
                self.assertEqual(row["toffoli_depth"], values[4])


if __name__ == "__main__":
    unittest.main(verbosity=2)
