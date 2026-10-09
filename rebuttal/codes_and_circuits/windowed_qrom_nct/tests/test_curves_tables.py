from __future__ import annotations

import tempfile
import unittest
import json
from pathlib import Path

from windowed_qrom_nct.basis import l_to_poly, poly_to_l
from windowed_qrom_nct.curves import CURVES, affine_add, gf_inverse, gf_multiply, on_curve
from windowed_qrom_nct.tables import (
    PUBLIC_SEED,
    curve_manifest,
    descriptor,
    descriptors,
    iter_table,
    verify_manifest,
    write_manifest,
    write_manifest_index,
)
from windowed_qrom_nct.openssl_ec import PointGenerator


class CurveAndTableTests(unittest.TestCase):
    def test_standard_generators_are_on_curve(self) -> None:
        for curve in CURVES.values():
            self.assertTrue(on_curve(curve.gx, curve.gy, curve), curve.openssl_name)

    def test_manifest_fixes_both_p_and_q_coordinates(self) -> None:
        for curve in CURVES.values():
            record = curve_manifest(curve)
            qx = int(record["Q"]["x"], 16)
            qy = int(record["Q"]["y"], 16)
            self.assertEqual(record["Q_scalar_relative_to_P"], 7)
            self.assertTrue(on_curve(qx, qy, curve))

    def test_basis_round_trip(self) -> None:
        for n, curve in CURVES.items():
            for value in (0, 1, curve.a, curve.gx, curve.gy):
                self.assertEqual(l_to_poly(poly_to_l(value, n), n), value)

    def test_public_seed_is_frozen(self) -> None:
        self.assertEqual(
            PUBLIC_SEED.hex(),
            "de0faddf37eaf09c24fe74b8dc130a5f75e73234f62f1ca55dc66b7607a67340",
        )

    def test_descriptor_manifest_fixes_total_translation(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "manifest.json"
            payload = write_manifest(path, n=163, s=9, include_table_stats=False)
            expected = sum(row["delta_scalar"] for row in payload["tables"]) % CURVES[163].order
            self.assertEqual(payload["delta_total"]["scalar_relative_to_P"], expected)
            self.assertEqual(verify_manifest(path)["tables"], 38)

    def test_short_top_window(self) -> None:
        # m=164, s=9 gives eighteen full windows and a two-bit top window.
        rows = descriptors(163, 9)
        self.assertEqual(len(rows), 38)
        self.assertEqual(rows[18].window_bits, 2)
        self.assertEqual(rows[18].entries, 4)
        self.assertEqual(rows[-1].window_bits, 2)

    def test_small_shifted_table_has_no_infinity_and_starts_at_delta(self) -> None:
        desc = descriptor(CURVES[163], 2, "P", 0)
        table = list(iter_table(desc))
        self.assertEqual(len(table), 4)
        self.assertTrue(all(on_curve(x, y, CURVES[163]) for x, y in table))

    def test_public_q_coordinates_need_no_known_discrete_logarithm(self) -> None:
        curve = CURVES[163]
        with PointGenerator(curve) as generator:
            [(qx, qy)] = list(generator.progression(11, 1, 1))
        q_point = (qx, qy)
        q_descriptor = descriptor(curve, 2, "Q", 0, q_point)
        self.assertIsNone(q_descriptor.base_scalar)
        self.assertEqual(q_descriptor.scalar_reference, "Q")
        table = list(iter_table(q_descriptor))
        self.assertEqual(len(table), 4)
        self.assertTrue(all(on_curve(x, y, curve) for x, y in table))

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "custom-q.json"
            payload = write_manifest(
                path,
                n=163,
                s=18,
                include_table_stats=False,
                q_point=q_point,
            )
            self.assertIsNone(payload["curve"]["Q_scalar_relative_to_P"])
            self.assertEqual(int(payload["curve"]["Q"]["x"], 16), qx)
            self.assertEqual(int(payload["curve"]["Q"]["y"], 16), qy)
            self.assertEqual(verify_manifest(path)["tables"], 20)

    def test_manifest_index_can_freeze_the_active_study_scope(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_manifest(root / "n163_s2.json", n=163, s=2, include_table_stats=False)
            write_manifest(root / "n571_s2.json", n=571, s=2, include_table_stats=False)
            index_path = root / "study-index.json"
            result = write_manifest_index(root, index_path, field_sizes={163, 233, 283})
            self.assertEqual(result["field_sizes"], [163])
            self.assertEqual({row["n"] for row in result["manifests"]}, {163})
            self.assertEqual(json.loads(index_path.read_text())["field_sizes"], [163])

    def test_early_unlookup_formula_uses_only_r_z_v(self) -> None:
        curve = CURVES[163]
        table = list(iter_table(descriptor(curve, 2, "P", 0)))
        r = (curve.gx, curve.gy)
        addend = table[2]
        z = r[0] ^ addend[0]
        v = r[1] ^ addend[1]
        slope = gf_multiply(v, gf_inverse(z, curve), curve)
        x3 = gf_multiply(slope, slope, curve) ^ slope ^ z ^ curve.a
        y3 = gf_multiply(slope, r[0] ^ x3, curve) ^ x3 ^ r[1]
        self.assertEqual((x3, y3), affine_add(r, addend, curve))
        self.assertTrue(on_curve(x3, y3, curve))


if __name__ == "__main__":
    unittest.main()
