from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from windowed_qrom_nct.analytic import DEFAULT_RESULTS, FIELD_SIZES, run


EXPECTED = {
    (163, "minimum_toffoli"): (10, 34, 859352, 172273452, 126944, 4416, 283806),
    (163, "minimum_tdw"): (6, 56, 1140472, 169570892, 204856, 608, 264086),
    (163, "minimum_dw_barrier"): (8, 27, 927068, 155263818, 114260, 3372, 289358),
    (233, "minimum_toffoli"): (10, 48, 1673988, 441495408, 261705, 4440, 386846),
    (233, "minimum_tdw"): (5, 94, 2768248, 568697540, 501605, 560, 427272),
    (233, "minimum_dw_barrier"): (9, 32, 1739348, 432758628, 220895, 6448, 416854),
    (283, "minimum_toffoli"): (11, 52, 2434896, 862161814, 352297, 8568, 550770),
    (283, "minimum_tdw"): (5, 114, 4156864, 1026739070, 753590, 596, 560740),
    (283, "minimum_dw_barrier"): (8, 15, 2921228, 814667288, 247916, 9560, 731498),
    (571, "minimum_toffoli"): (11, 104, 9202596, 5703468882, 1472990, 8688, 1279034),
    (571, "minimum_tdw"): (6, 192, 14822796, 7441985718, 2685300, 848, 1395718),
    (571, "minimum_dw_barrier"): (10, 30, 10157204, 6145473890, 935355, 29092, 1730120),
}


class AnalyticScanTests(unittest.TestCase):
    def test_upper_model_reproduces_window_tex_rows_and_coverage(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            summary = run(DEFAULT_RESULTS, output, ("upper",))
            self.assertEqual(summary["field_sizes"], list(FIELD_SIZES))
            fields = summary["models"]["upper"]["fields"]
            self.assertEqual(
                {n: fields[str(n)]["evaluated_points"] for n in FIELD_SIZES},
                {163: 832, 233: 1178, 283: 1428, 571: 1734},
            )
            self.assertEqual(fields["571"]["grid_points"], 2866)
            self.assertEqual(fields["571"]["unavailable_points"], 1132)

            for (n, objective), expected in EXPECTED.items():
                row = fields[str(n)]["objectives_over_available_points"][objective]
                actual = (
                    row["s"], row["p"], row["toffoli"], row["cnot"],
                    row["width"], row["toffoli_depth"], row["nct_depth_barrier"],
                )
                self.assertEqual(actual, expected, (n, objective))
                self.assertEqual(row["dw_barrier"], row["width"] * row["nct_depth_barrier"])
                self.assertEqual(row["tdw"], row["width"] * row["toffoli_depth"])

            minima = output / "windowed_theory_upper_minima.csv"
            with minima.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 12)
            loaded = json.loads((output / "summary.json").read_text(encoding="utf-8"))
            self.assertFalse(loaded["gate_generation"])
            self.assertEqual(loaded["models"]["upper"]["rows"], 6304)


if __name__ == "__main__":
    unittest.main()

