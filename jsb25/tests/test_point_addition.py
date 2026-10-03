"""Complete n=163 regression plus four-parameter interface checks."""

from __future__ import annotations

import unittest

from init.config import FIELD_CONFIGS, SUPPORTED_FIELDS, get_config
from point_addition import IN_PLACE, OUT_OF_PLACE, build_point_addition
from shor_table_resources import estimate_shor_resources
from verify_point_addition import run_verification


class ConfigurationTests(unittest.TestCase):
    def test_supported_fields_are_exactly_the_requested_four(self) -> None:
        self.assertEqual(SUPPORTED_FIELDS, (163, 233, 283, 571))
        self.assertEqual(tuple(FIELD_CONFIGS), SUPPORTED_FIELDS)

    def test_removed_small_parameter_interfaces_are_rejected(self) -> None:
        for n in (8, 16, 127):
            with self.subTest(n=n), self.assertRaises(ValueError):
                get_config(n)
            with self.subTest(n=n), self.assertRaises(ValueError):
                build_point_addition(n, "inplace")

    def test_width_formulas_for_all_supported_fields(self) -> None:
        for n in SUPPORTED_FIELDS:
            cfg = get_config(n)
            in_width = (
                1
                + 5 * n
                + cfg.karatsuba_ancillas
                + (cfg.inversion_multiplications + 1)
                * cfg.multiplication_targets
            )
            out_width = (
                1
                + 5 * n
                + cfg.karatsuba_ancillas
                + (cfg.inversion_multiplications + 2)
                * cfg.multiplication_targets
            )
            self.assertGreater(out_width, in_width)


class CompletePointAdditionTests(unittest.TestCase):
    def test_n163_both_modes_and_both_control_values(self) -> None:
        # n=163 is the default CI-sized full-stream test.  The same verifier
        # accepts 233/283/571 from its CLI, one field at a time.
        run_verification(163, "both", check_optimized=True)

    def test_n163_shor_resource_interfaces(self) -> None:
        for mode in (IN_PLACE, OUT_OF_PLACE):
            with self.subTest(mode=mode):
                estimate = estimate_shor_resources(163, mode)
                self.assertEqual(estimate.field_degree, 163)
                self.assertEqual(estimate.mode, mode)
                self.assertGreater(estimate.toffoli_count, 0)
                self.assertGreater(estimate.cnot_count, 0)
                self.assertGreater(estimate.full_depth, 0)
                self.assertGreater(estimate.toffoli_depth, 0)
                self.assertEqual(
                    estimate.dw_cost,
                    estimate.qubits * estimate.full_depth,
                )


if __name__ == "__main__":
    unittest.main()
