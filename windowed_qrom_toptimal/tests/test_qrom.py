import unittest
from support.qrom import QromRegisters,emit_unary_qrom
from support.scheduler import PerLineAsap,ClassicalGateSimulator
from model import qrom


class QromTests(unittest.TestCase):
    def test_lookup_upper_formula_n163(self):
        n=163
        for b in range(1,7):
            bus=2*n
            regs=QromRegisters(tuple(range(b)),tuple(range(b,2*b-1)),
                tuple(range(2*b-1,2*b-1+bus)),tuple(range(2*b-1+bus,2*b-2+2*bus)))
            sink=PerLineAsap(2*b-2+2*bus)
            emit_unary_qrom(sink,regs,[(1<<bus)-1]*(1<<b))
            actual=sink.stats(); expected=qrom(n,b)
            self.assertEqual(actual.toffoli,expected["toffoli"])
            self.assertEqual(actual.toffoli_depth,expected["toffoli_depth"])
            self.assertEqual(actual.cnot,expected["cnot"])
            self.assertLessEqual(actual.full_depth,expected["full_depth"])

    def test_early_unlookup_and_inverse(self):
        # Independent computational-basis audit of all toy addresses.
        regs=QromRegisters((0,1),(2,),(3,4,5,6),(15,16,17))
        words=(1,6,9,15)
        for address,word in enumerate(words):
            initial=address | (1<<11) | (1<<14)
            sink=ClassicalGateSimulator(18,initial)
            emit_unary_qrom(sink,regs,words)
            for i in range(2):
                sink.cnot(11+i,7+i); sink.cnot(3+i,7+i)
                sink.cnot(13+i,9+i); sink.cnot(5+i,9+i)
            emit_unary_qrom(sink,regs,words,inverse=True)
            self.assertEqual((sink.state>>3)&15,0)
            self.assertEqual((sink.state>>2)&1,0)
            self.assertEqual(sink.state>>15,0)
            self.assertEqual((sink.state>>7)&3,1^(word&3))
            self.assertEqual((sink.state>>9)&3,2^(word>>2))


if __name__=="__main__": unittest.main()
