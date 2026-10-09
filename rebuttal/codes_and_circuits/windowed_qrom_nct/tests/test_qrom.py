from __future__ import annotations

import unittest

from windowed_qrom_nct.qrom import QromRegisters, emit_unary_qrom
from windowed_qrom_nct.scheduler import ClassicalGateSimulator, PerLineAsap


class QromTests(unittest.TestCase):
    WORDS = (0b000001, 0b010010, 0b101100, 0b111111)
    REGS = QromRegisters(address=(0, 1), selector_ancillas=(2,), bus=(3, 4, 5, 6, 7, 8))

    def test_all_addresses_and_literal_inverse(self) -> None:
        for address, word in enumerate(self.WORDS):
            sim = ClassicalGateSimulator(9, address)
            emit_unary_qrom(sim, self.REGS, self.WORDS)
            self.assertEqual((sim.state >> 3) & 0x3F, word)
            self.assertEqual(sim.state & 0b111, address)
            emit_unary_qrom(sim, self.REGS, self.WORDS, inverse=True)
            self.assertEqual(sim.state, address)

    def test_coherent_gate_counts(self) -> None:
        scheduler = PerLineAsap(9)
        emit_unary_qrom(scheduler, self.REGS, self.WORDS)
        stats = scheduler.stats()
        self.assertEqual(stats.toffoli, 2 * (len(self.WORDS) - 2))
        self.assertEqual(stats.cnot, len(self.WORDS) - 2 + sum(x.bit_count() for x in self.WORDS))

    def test_forward_then_inverse_doubles_non_x_counts(self) -> None:
        scheduler = PerLineAsap(9)
        emit_unary_qrom(scheduler, self.REGS, self.WORDS)
        first = scheduler.stats()
        emit_unary_qrom(scheduler, self.REGS, self.WORDS, inverse=True)
        both = scheduler.stats()
        self.assertEqual(both.toffoli, 2 * first.toffoli)
        self.assertEqual(both.cnot, 2 * first.cnot)

    def test_balanced_data_fanout_reduces_depth_and_has_exact_count(self) -> None:
        regs = QromRegisters(
            address=(0, 1),
            selector_ancillas=(2,),
            bus=(3, 4, 5, 6, 7, 8),
            data_fanout=(9, 10, 11, 12, 13),
        )
        scheduler = PerLineAsap(14)
        emit_unary_qrom(scheduler, regs, self.WORDS)
        stats = scheduler.stats()
        nonzero = sum(word != 0 for word in self.WORDS)
        weight = sum(word.bit_count() for word in self.WORDS)
        expected_data_cnot = 3 * weight - 2 * nonzero
        self.assertEqual(stats.cnot, len(self.WORDS) - 2 + expected_data_cnot)

        for address, word in enumerate(self.WORDS):
            sim = ClassicalGateSimulator(14, address)
            emit_unary_qrom(sim, regs, self.WORDS)
            self.assertEqual((sim.state >> 3) & 0x3F, word)
            emit_unary_qrom(sim, regs, self.WORDS, inverse=True)
            self.assertEqual(sim.state, address)

    def test_early_unlookup_copies_clean_z_v_and_clears_cache(self) -> None:
        # Two-bit x and two-bit y toy coordinates packed as x | (y << 2).
        words = (0b0000, 0b0110, 0b1001, 0b1111)
        regs = QromRegisters(
            address=(0, 1),
            selector_ancillas=(2,),
            bus=(3, 4, 5, 6),
            data_fanout=(15, 16, 17),
        )
        # z=(7,8), v=(9,10), R_x=01, R_y=10.
        for address, word in enumerate(words):
            initial = address | (1 << 11) | (1 << 14)
            sim = ClassicalGateSimulator(18, initial)
            emit_unary_qrom(sim, regs, words)
            for bit in range(2):
                sim.cnot(11 + bit, 7 + bit)
                sim.cnot(3 + bit, 7 + bit)
                sim.cnot(13 + bit, 9 + bit)
                sim.cnot(5 + bit, 9 + bit)
            emit_unary_qrom(sim, regs, words, inverse=True)
            self.assertEqual((sim.state >> 3) & 0b1111, 0)
            self.assertEqual((sim.state >> 2) & 1, 0)  # selector is clean
            self.assertEqual((sim.state >> 15) & 0b111, 0)  # fanout is clean
            expected_z = 0b01 ^ (word & 0b11)
            expected_v = 0b10 ^ ((word >> 2) & 0b11)
            self.assertEqual((sim.state >> 7) & 0b11, expected_z)
            self.assertEqual((sim.state >> 9) & 0b11, expected_v)


if __name__ == "__main__":
    unittest.main()
