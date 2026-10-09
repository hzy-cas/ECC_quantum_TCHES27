from __future__ import annotations

import csv
import json
import tempfile
import unittest
from dataclasses import asdict
from pathlib import Path

from windowed_qrom_nct.search import grid_points, index_artifacts, write_coverage
from windowed_qrom_nct.synthesis import ExactDesign


class ExactSearchAuditTests(unittest.TestCase):
    def _design(self) -> ExactDesign:
        return ExactDesign(
            model_version="test-flow",
            arithmetic="ac",
            n=163,
            scalar_bits=164,
            s=18,
            p=20,
            addends=20,
            phase1_layers=0,
            tree_layers=5,
            toffoli=1,
            cnot=2,
            width=3,
            toffoli_depth=4,
            nct_depth=5,
            current_depth=6,
            dw=15,
            tdw=12,
            table_seed="seed",
            data_basis="AC L-basis",
            qrom_data_write="test",
            scheduler="test",
            layer_records=(),
        )

    def test_grid_enumerates_every_lane_for_each_window(self) -> None:
        points = list(grid_points(163, 17, 18))
        # ceil(164/17)*2 = 20 and ceil(164/18)*2 = 20.
        self.assertEqual(len(points), 40)
        self.assertEqual(points[0], (17, 1))
        self.assertEqual(points[-1], (18, 20))

    def test_index_and_coverage_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            artifact = root / "n163_s18_p20_ac.json"
            artifact.write_text(json.dumps(asdict(self._design())), encoding="utf-8")
            index = root / "search.csv"
            self.assertEqual(
                index_artifacts(
                    n=163, arithmetic="ac", artifact_dir=root, output=index
                ),
                1,
            )
            with index.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual((rows[0]["s"], rows[0]["p"]), ("18", "20"))

            coverage = write_coverage(
                root / "coverage.json",
                n=163,
                arithmetic="ac",
                search_csv=index,
                s_min=18,
                s_max=18,
            )
            self.assertFalse(coverage["complete"])
            self.assertEqual(coverage["expected_points"], 20)
            self.assertEqual(coverage["exact_points"], 1)
            self.assertEqual(len(coverage["missing_points"]), 19)


if __name__ == "__main__":
    unittest.main()
