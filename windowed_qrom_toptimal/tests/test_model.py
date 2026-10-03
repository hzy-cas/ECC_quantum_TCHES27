import unittest
from collections import Counter
import model


class ModelTests(unittest.TestCase):
    def test_n163_grid_and_winner(self):
        self.assertEqual(len(list(model.grid(163))), 17*328)
        candidates = model.formula_minimizers(163)
        self.assertEqual(len(candidates), 16*19)
        self.assertEqual({(r['s'],r['p']) for r in candidates},
                         {(s,p) for s in range(3,19) for p in range(110,129)})
        self.assertTrue(all(r['predicted_toffoli_depth']==248 and r['w']==r['p'] for r in candidates))
        self.assertEqual(model.required_layers(163, 3, 110),
                         {(163, "reduction", k) for k in (1, 2, 3, 7, 14, 27, 55)})

    def test_short_windows_and_odd_tree(self):
        self.assertEqual(model.tables(163, 3, 110), [3]*108+[2,2])
        tree = [payload for kind, payload in model.schedule(163, 3, 110) if kind=="reduction"]
        self.assertEqual(tree, [55,27,14,7,3,2,1])
        self.assertEqual(sum(tree),109)

    def test_n233_regression_not_two_separate_scalars(self):
        # 468 original bits, 117 strided lanes: 117 four-bit tables, not 118.
        self.assertEqual(model.tables(233,4,117),[4]*117)
        self.assertEqual(model.window_indices(233,4,117,0,0),(0,117,234,351))
        self.assertEqual(Counter(model.tables(233,4,118)),{4:114,3:4})
        self.assertEqual(model.predicted_depth(233,4,117),280)
        self.assertNotIn('tail',[kind for kind,_ in model.schedule(233,4,117)])

    def test_partial_round_keeps_only_active_lanes(self):
        self.assertEqual(model.window_batches(163,3,109),[(3,)*109,(1,)])
        self.assertIn((163,'tail',1),model.required_layers(163,3,109))

    def test_manuscript_index_partition_all_fields(self):
        # Only geometry is checked for larger n: no large gate streams run.
        for n in model.FIELDS:
            N=2*n+2
            for s in (2,3,4,7,18):
                for p in sorted({1,2,3,N//s,(N+s-1)//s,N-1,N}):
                    batches=model.window_batches(n,s,p)
                    seen=[]
                    for q,bits_per_lane in enumerate(batches):
                        for j,bits in enumerate(bits_per_lane):
                            indices=model.window_indices(n,s,p,j,q)
                            self.assertEqual(len(indices),bits)
                            self.assertEqual(indices,tuple((q*s+t)*p+j for t in range(bits)))
                            seen.extend(indices)
                    self.assertEqual(sorted(seen),list(range(N)))
                    self.assertEqual(len(batches)-1,(N+p*s-1)//(p*s)-1)
                    stages=model.schedule(n,s,p)
                    tree=[k for kind,k in stages if kind=='reduction']
                    self.assertEqual(len(tree),model.ceil_log2(p))
                    total_additions=sum(len(bits) for bits in batches[1:])+sum(tree)
                    self.assertEqual(total_additions,sum(map(len,batches))-1)

    def test_search_includes_full_parallelism_not_old_window_bound(self):
        self.assertIn((18,328),list(model.grid(163)))
        self.assertEqual(model.tables(163,18,328),[1]*328)
        self.assertEqual(model.qrom(163,1)['toffoli'],0)

    def test_missing_layers_never_emit_resources(self):
        row = model.evaluate(163,3,110,{})
        self.assertFalse(row["available"])
        self.assertNotIn("width",row)
        self.assertNotIn("toffoli_depth",row)
        self.assertEqual(row["predicted_toffoli_depth"],248)

    def test_retained_workspace_is_not_clean_workspace(self):
        reduction=model.liveness(163,3,"reduction")
        tail=model.liveness(163,3,"tail")
        self.assertEqual(reduction["retained_new_qubits"]-tail["retained_new_qubits"],2*163*3)
        self.assertEqual(reduction["clean_qubits"],tail["clean_qubits"])
        self.assertEqual(tail["qubits"],tail["input_qubits"]+tail["clean_qubits"]+tail["retained_new_qubits"])
        self.assertEqual(tail["output_alias_input_qubits"],163*3)

    def test_foreign_cache_rejected(self):
        with self.assertRaises(ValueError):
            model.validate_record(dict(model_version="ac-balanced-clean-v1"),163,"reduction",1)

    def test_parameter_bounds(self):
        for s,p in ((1,1),(19,1),(3,0),(3,329)):
            with self.assertRaises(ValueError): model.schedule(163,s,p)


if __name__=="__main__": unittest.main()
