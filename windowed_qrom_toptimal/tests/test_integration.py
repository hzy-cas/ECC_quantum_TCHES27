"""Real C++ gate streams run only for n=163."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import model
import run
from support.fixtures import make_fixture
from support.basis import poly_to_l


@unittest.skipUnless(run.EXECUTABLE.is_file(), "run python3 run.py build first")
class IntegrationTests(unittest.TestCase):
    def test_parallel_generation_verification_and_resume(self):
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)
            command=[sys.executable,str(run.ROOT/"run.py"),"run","--n","163","--scope","smoke",
                "--jobs","2","--memory-budget-gb","3","--verify","--output",str(output)]
            first=subprocess.run(command,capture_output=True,text=True,check=True)
            files=sorted((output/"layers").glob("*.json"))
            self.assertEqual(len(files),6)
            snapshots={path.name:path.read_bytes() for path in files}
            for path in files:
                row=json.loads(path.read_text())
                self.assertTrue(row["basis_verified"])
                model.validate_record(row,163,row["kind"],row["lanes"])
            second=subprocess.run(command,capture_output=True,text=True,check=True)
            self.assertIn("pending: 0",second.stdout)
            self.assertEqual(snapshots,{path.name:path.read_bytes() for path in files})
            self.assertNotIn("FAILED",first.stderr)
            path=files[0]
            corrupted=json.loads(path.read_text()); corrupted["fingerprint"]="foreign"
            path.write_text(json.dumps(corrupted))
            rejected=subprocess.run(command,capture_output=True,text=True)
            self.assertNotEqual(rejected.returncode,0)
            self.assertIn("stale/foreign cache",rejected.stderr)

    def test_bad_fixture_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"bad.txt"
            fields=make_fixture(163,1).split()
            fields[-1]=hex(int(fields[-1],16)^1)
            path.write_text(" ".join(fields))
            command=[str(run.EXECUTABLE),"--n","163","--lanes","1","--kind","reduction",
                "--data-root",str(run.ROOT/"data"),"--curve-a-hex",hex(poly_to_l(1,163)),
                "--verify-fixture",str(path)]
            result=subprocess.run(command,capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn("output mismatch",result.stderr)

    def test_legacy_layers_reused_but_network_summary_recomputed(self):
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)
            command=[sys.executable,str(run.ROOT/'run.py'),'run','--n','163','--scope','smoke',
                     '--jobs','2','--output',str(output)]
            subprocess.run(command,capture_output=True,text=True,check=True)
            files=sorted((output/'layers').glob('*.json'))
            legacy=next(iter(run.LEGACY_LAYER_IDENTITIES))
            for path in files:
                record=json.loads(path.read_text()); record['fingerprint']=legacy
                path.write_text(json.dumps(record))
            snapshots={p.name:p.read_bytes() for p in files}
            (output/'summary.json').write_text(json.dumps({'model_version':model.LAYER_VERSION,
                'fields':{'163':{'minimum_toffoli_depth_over_available':{'toffoli_depth':-1}}}}))
            second=subprocess.run(command,capture_output=True,text=True,check=True)
            self.assertIn('pending: 0',second.stdout)
            self.assertEqual(snapshots,{p.name:p.read_bytes() for p in files})
            summary=json.loads((output/'summary.json').read_text())
            self.assertEqual(summary['model_version'],model.MODEL_VERSION)
            self.assertEqual(summary['windowing_model'],'lane_first')
            self.assertEqual(summary['fields']['163']['grid_points'],17*328)
            self.assertGreater(summary['fields']['163']['minimum_toffoli_depth_over_available']['toffoli_depth'],0)


if __name__=="__main__": unittest.main()
