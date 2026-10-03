"""Identity separation, portable build orchestration, and all-field plans."""
import argparse
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import model
import run


class RunnerTests(unittest.TestCase):
    def test_legacy_layers_only_reused_for_exact_unchanged_arithmetic(self):
        legacy,current=next(iter(run.LEGACY_LAYER_IDENTITIES.items()))
        self.assertEqual(run.layer_fingerprint(),current)
        self.assertTrue(run.compatible_layer_identity(legacy,current))
        self.assertTrue(run.compatible_layer_identity(current,current))
        self.assertFalse(run.compatible_layer_identity('foreign',current))
        self.assertFalse(run.compatible_layer_identity(legacy,'changed-backend'))

    def test_network_change_does_not_change_arithmetic_identity(self):
        previous_layer=run.layer_fingerprint()
        previous_network=run.fingerprint()
        original=run.sha256
        with patch('run.sha256',side_effect=lambda p:'changed-model' if p.name=='model.py' else original(p)):
            self.assertEqual(run.layer_fingerprint(),previous_layer)
            self.assertNotEqual(run.fingerprint(),previous_network)

    def test_winners_plan_fits_requested_200_gib_budget(self):
        args=argparse.Namespace(n=list(model.FIELDS),s_min=2,s_max=18,scope='winners')
        tasks,candidates=run.task_plan(args)
        self.assertEqual(len(tasks),74)
        self.assertEqual(len(candidates),971)
        self.assertEqual(max(map(run.estimated_gb,tasks)),193)
        self.assertTrue(all(kind=='reduction' for _,kind,_ in tasks))
        self.assertTrue(all(row['w']==row['p'] for row in candidates))
        for n,depth in ((163,248),(233,280),(283,308),(571,380)):
            self.assertEqual({r['predicted_toffoli_depth'] for r in candidates if r['n']==n},{depth})

    def test_auto_build_before_workers(self):
        with patch('run.check_build',side_effect=[OSError(8,'Exec format error'),None]) as check:
            with patch('run.build_project') as build:
                run.ensure_build(4)
                build.assert_called_once_with(4)
                self.assertEqual(check.call_count,2)

    def test_build_uses_ctest_316_compatible_working_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            with patch('run.ROOT',root),patch('run.EXECUTABLE',root/'build'/'retained_estimator'):
                with patch('run.shutil.which',return_value='/usr/bin/tool'),patch('run.subprocess.run') as command:
                    with patch('run.sha256',return_value='exe-hash'),patch('run.layer_fingerprint',return_value='layer-hash'):
                        run.build_project(2)
                self.assertEqual(command.call_args_list[-1].args[0],['ctest','--output-on-failure'])
                self.assertEqual(command.call_args_list[-1].kwargs['cwd'],root/'build')
                stamp=json.loads((root/'build'/'build_stamp.json').read_text())
                self.assertEqual(stamp['fingerprint'],'layer-hash')
                self.assertEqual(stamp['executable_sha256'],'exe-hash')

    def test_incompatible_build_preserved_before_reconfiguration(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            old=root/'build'; old.mkdir()
            (old/'CMakeCache.txt').write_text('CMAKE_HOME_DIRECTORY:INTERNAL=/a/different/computer\n')
            with patch('run.ROOT',root),patch('run.EXECUTABLE',old/'retained_estimator'):
                with patch('run.shutil.which',return_value='/usr/bin/tool'),patch('run.subprocess.run'):
                    with patch('run.sha256',return_value='exe-hash'),patch('run.layer_fingerprint',return_value='layer-hash'):
                        run.build_project(2)
            backups=list(root.glob('build.backup.*/build/CMakeCache.txt'))
            self.assertEqual(len(backups),1)
            self.assertIn('/a/different/computer',backups[0].read_text())


if __name__=='__main__': unittest.main()
