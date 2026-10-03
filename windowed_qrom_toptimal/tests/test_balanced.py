"""Balanced TDW search; actual gate generation is restricted to n=163."""
import argparse
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
from support.curves import CURVES


def synthetic_records(n, objective):
    """Only for testing the allocator/composer, never exported as results."""
    result = {}
    for kind in ("tail", "reduction"):
        for k in range(1, 2*n+3):
            depth = model.arithmetic_depth(n,k,objective)
            result[n,kind,k] = dict(model_version=model.layer_version(objective),n=n,kind=kind,lanes=k,
                **model.liveness(n,k,kind,objective),toffoli=model.arithmetic_count(n,k,objective),
                toffoli_depth=depth,cnot=0,full_depth=10*depth,current_depth=10*depth)
    return result


class BalancedModelTests(unittest.TestCase):
    def test_liveness_cleans_targets_and_square_but_retains_copies(self):
        for n in model.FIELDS:
            m,r=model.PARAMETERS[n]
            for k in (1,2,3,7,8,17):
                for kind in ("tail","reduction"):
                    a=model.liveness(n,k,kind,'tdw')
                    self.assertEqual(a['clean_qubits'],4*(m-n)*k+m*(r+3*(k-1))+4*n+n*(k//2))
                    self.assertEqual(a['retained_new_qubits'],(3 if kind=='tail' else 4)*n*k)
                    self.assertEqual(a['qubits'],a['input_qubits']+a['clean_qubits']+a['retained_new_qubits'])
                    self.assertEqual(model.arithmetic_depth(n,k,'tdw'),2*model.arithmetic_depth(n,k))
                    self.assertEqual(model.arithmetic_count(n,k,'tdw'),2*model.arithmetic_count(n,k))

    def test_tdw_minimizes_product_not_depth(self):
        rows=[model.evaluate(163,s,p,synthetic,'tdw')
              for synthetic in [synthetic_records(163,'tdw')] for s,p in model.grid(163)]
        best=min(row['tdw'] for row in rows)
        candidates=model.formula_minimizers(163,objective='tdw')
        self.assertEqual({(r['s'],r['p']) for r in candidates},
                         {(r['s'],r['p']) for r in rows if r['tdw']==best})
        self.assertEqual(best,120399200)
        self.assertEqual(candidates[0],dict(n=163,s=6,p=55,w=55,
            predicted_toffoli_depth=608,predicted_width=198025,predicted_tdw=120399200))
        # Minimum TD in the SAME Balanced family differs from minimum TDW.
        self.assertLess(min(r['toffoli_depth'] for r in rows),608)

    def test_qrom_and_window_schedule_unchanged(self):
        n,s,p=163,3,109  # includes a partial final accumulation round
        rows=[model.evaluate(n,s,p,synthetic_records(n,o),o) for o in ('td','tdw')]
        for key in ('original_addends','addends','accumulation_layers','tree_layers','windowing_model'):
            self.assertEqual(rows[0][key],rows[1][key])
        stages=model.schedule(n,s,p)
        arithmetic=sum(model.arithmetic_count(n,len(v) if kind=='tail' else v)
                       for kind,v in stages if kind!='initialization')
        self.assertEqual(rows[1]['toffoli']-rows[0]['toffoli'],2*arithmetic)
        qrom_forward=sum(sum(model.qrom(n,b,kind=='tail')['toffoli'] for b in v)
                         for kind,v in stages if kind!='reduction')
        self.assertEqual(rows[0]['toffoli'],2*(arithmetic+qrom_forward))
        self.assertEqual(rows[1]['toffoli'],2*(2*arithmetic+qrom_forward))
        self.assertLess(rows[1]['width'],rows[0]['width'])

    def test_all_field_plans_formula_only(self):
        args=argparse.Namespace(n=list(model.FIELDS),s_min=2,s_max=18,scope='winners',objective='tdw')
        tasks,candidates=run.task_plan(args)
        self.assertEqual(len(tasks),25)
        self.assertEqual(len(candidates),53)
        self.assertEqual(max(run.estimated_gb(t,'tdw') for t in tasks),191)
        for n,s,w,td,width in [(163,6,55,608,198025),(233,5,94,560,501605),
                               (283,6,95,716,624661),(571,6,191,848,2658604)]:
            row=next(r for r in candidates if r['n']==n)
            self.assertEqual((row['s'],row['w'],row['predicted_toffoli_depth'],row['predicted_width']),
                             (s,w,td,width))
            self.assertTrue(all(r['predicted_tdw']==td*width for r in candidates if r['n']==n))

    def test_missing_records_do_not_become_complete_results(self):
        row=model.evaluate(163,6,55,{},'tdw')
        self.assertFalse(row['available'])
        self.assertEqual(row['predicted_tdw'],120399200)
        for key in ('toffoli','toffoli_depth','width','cnot','tdw','dw','nct_depth_barrier'):
            self.assertNotIn(key,row)

    def test_cross_architecture_records_rejected(self):
        old=synthetic_records(163,'td')[163,'reduction',1]
        new=synthetic_records(163,'tdw')[163,'reduction',1]
        with self.assertRaisesRegex(ValueError,'mismatched'):
            model.validate_record(old,163,'reduction',1,'tdw')
        with self.assertRaisesRegex(ValueError,'mismatched'):
            model.validate_record(new,163,'reduction',1,'td')


@unittest.skipUnless(run.EXECUTABLE.is_file(), 'run python3 run.py build first')
class BalancedIntegrationTests(unittest.TestCase):
    def test_winners_n163_complete_resources(self):
        with tempfile.TemporaryDirectory() as directory:
            command=[sys.executable,str(run.ROOT/'run.py'),'run','--n','163',
                '--scope','winners','--objective','tdw','--jobs','1','--verify',
                '--memory-budget-gb','6','--output',directory]
            subprocess.run(command,capture_output=True,text=True,check=True)
            summary=json.loads((Path(directory)/'summary.json').read_text())
            info=summary['fields']['163']
            self.assertTrue(info['formula_minimizers_all_available'])
            self.assertFalse(info['complete_full_resource_grid'])
            self.assertEqual(len(info['formula_minimizers']),13)
            best=info['selected_optimum']
            expected=dict(s=6,w=55,toffoli=1122328,cnot=167316138,width=198025,
                toffoli_depth=608,nct_depth_barrier=264098,dw=52298006450,tdw=120399200)
            self.assertEqual({key:best[key] for key in expected},expected)
            self.assertEqual(best,info['minimum_tdw_over_available'])
            layers=list((Path(directory)/'layers').glob('*.json'))
            self.assertEqual(len(layers),6)
            self.assertTrue(all(json.loads(path.read_text())['basis_verified'] for path in layers))

    def test_basis_and_clean_workspace_odd_even_batches_n163(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture=Path(directory)/'fixture.txt'
            for k in (1,2,3,4,5,7,8,9):
                fixture.write_text(make_fixture(163,k))
                for kind in ('tail','reduction'):
                    command=[str(run.EXECUTABLE),'--n','163','--lanes',str(k),'--kind',kind,
                        '--architecture','balanced','--data-root',str(run.ROOT/'data'),
                        '--curve-a-hex',hex(poly_to_l(CURVES[163].a,163)),
                        '--verify-fixture',str(fixture)]
                    row=json.loads(subprocess.run(command,capture_output=True,text=True,check=True).stdout)
                    model.validate_record(row,163,kind,k,'tdw')
                    self.assertTrue(row['basis_verified'])

    def test_balanced_resume_and_output_isolation(self):
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)
            command=[sys.executable,str(run.ROOT/'run.py'),'run','--n','163',
                '--scope','smoke','--objective','tdw','--jobs','2','--verify','--output',str(output)]
            first=subprocess.run(command,capture_output=True,text=True,check=True)
            self.assertIn('all TDW formula minimizers available: False',first.stdout)
            summary=json.loads((output/'summary.json').read_text())
            info=summary['fields']['163']
            self.assertEqual(summary['point_addition'],'balanced')
            self.assertEqual(info['selected_optimum'],info['minimum_tdw_over_available'])
            files=list((output/'layers').glob('*.json'))
            snapshots={p.name:p.read_bytes() for p in files}
            second=subprocess.run(command,capture_output=True,text=True,check=True)
            self.assertIn('pending: 0',second.stdout)
            self.assertEqual(snapshots,{p.name:p.read_bytes() for p in files})
            before=(output/'summary.json').read_bytes()
            command[command.index('tdw')]='td'
            wrong=subprocess.run(command,capture_output=True,text=True)
            self.assertNotEqual(wrong.returncode,0)
            self.assertIn('another objective',wrong.stderr)
            self.assertEqual(before,(output/'summary.json').read_bytes())
            self.assertEqual(snapshots,{p.name:p.read_bytes() for p in files})

    def test_balanced_wrong_expected_output_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture=Path(directory)/'bad.txt'
            values=make_fixture(163,1).split(); values[-1]=hex(int(values[-1],16)^1)
            fixture.write_text(' '.join(values))
            command=[str(run.EXECUTABLE),'--n','163','--lanes','1','--kind','tail',
                '--architecture','balanced','--data-root',str(run.ROOT/'data'),
                '--curve-a-hex',hex(poly_to_l(CURVES[163].a,163)), '--verify-fixture',str(fixture)]
            result=subprocess.run(command,capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('output mismatch',result.stderr)


if __name__=='__main__': unittest.main()
