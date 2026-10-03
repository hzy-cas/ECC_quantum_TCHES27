# algorithms/implementations.py
from .base_shor import BaseShor

# ==========================================
# GF(163) Implementation
# ==========================================
class Shor163(BaseShor):
    def __init__(self, data_root):
        super().__init__(n=163, block_size=743, nums2_const=9, data_root=data_root)

    def get_mul_total_size(self):
        return 906

    def inversion_logic(self, a, sqr, ancilla, ancilla2, count):
        sqr1, sqr2, sqr3, sqr4 = sqr[:163], sqr[163:326], sqr[326:489], sqr[489:652]
        anc_main = ancilla 
        anc_side = ancilla2
        
        n = 163
        mul_size = 906
        
        # Helper aliases
        def S(x, p): return self.Square(x, 2*n, p)
        def M(x, y, c_off): 
            res = self.mul(x, y, self.Toffoli_qubits[c_off:], 0)
            return res
        def C(x, y): self.CNOT_N(x, y)

        # 1 -> 2
        a, A_1_1 = S(a+sqr1, 1)
        A_2_0 = M(a+anc_main[:743], A_1_1+anc_main[743:], count)
        count += mul_size
        
        # Chain 1 (2->4)
        A_2_0, A_2_2 = S(A_2_0+sqr2, 2)
        a, A_1_1 = S(a+A_1_1, 1) # clear sqr1
        sqr1 = A_1_1
        A_4_0 = M(A_2_2+anc_main[:743], A_2_0+anc_main[743:], count)
        count += mul_size
        
        # 4 -> 8
        A_2_0, A_2_2 = S(A_2_0+A_2_2, 2)
        sqr2 = A_2_2
        A_4_0, A_4_4 = S(A_4_0+sqr1, 4)
        A_8_0 = M(A_4_0+anc_main[:743], A_4_4+anc_main[743:], count)
        count += mul_size

        # 8 -> 16
        A_4_0, A_4_4 = S(A_4_0+A_4_4, 4)
        sqr1 = A_4_4
        A_8_0, A_8_8 = S(A_8_0+sqr2, 8)
        A_16_0 = M(A_8_0+anc_main[:743], A_8_8+anc_main[743:], count)
        count += mul_size

        # 16 -> 32
        A_8_0, A_8_8 = S(A_8_0+A_8_8, 8)
        sqr2 = A_8_8
        A_16_0, A_16_16 = S(A_16_0+sqr1, 16)
        A_32_0 = M(A_16_0+anc_main[:743], A_16_16+anc_main[743:], count)
        count += mul_size
        
        # S - Chain 2 (Prep 2->34 side branch)
        C(A_2_0, sqr3)
        sqr3, A_2_32 = S(sqr3+sqr4, 32)
        C(A_2_0, sqr3)

        # Chain 1 (32 -> 64)
        A_16_0, A_16_16 = S(A_16_0+A_16_16, 16)
        sqr1 = A_16_16
        A_32_0, A_32_32 = S(A_32_0+sqr2, 32)
        
        # Merge (32 + 2 = 34)
        C(A_32_0, sqr3)
        # Main
        A_64_0 = M(A_32_0+anc_main[:743], A_32_32+anc_main[743:], count)
        count += mul_size
        # Side (A_2_32 * sqr3) -> A_34_0
        A_34_0 = M(A_2_32+anc_side[:743], sqr3+anc_side[743:], count)
        count += mul_size
        C(A_32_0, sqr3)
        
        # 64 -> 128
        A_32_0, A_32_32 = S(A_32_0+A_32_32, 32)
        sqr2 = A_32_32
        A_64_0, A_64_64 = S(A_64_0+sqr1, 64)
        A_128_0 = M(A_64_0+anc_main[:743], A_64_64+anc_main[743:], count)
        count += mul_size

        # Chain 1 Cleanup
        A_64_0, A_64_64 = S(A_64_0+A_64_64, 64)
        sqr1 = A_64_64 # clear sqr1
        
        # Final Merge Prep (128->1, 34->129)
        A_128_0, A_128_1 = S(A_128_0+sqr2, 1)
        A_34_0, A_34_129 = S(A_34_0+sqr3, 129)
        
        # Final Mul
        A_162_1 = M(A_34_129+anc_side[:743], A_128_1+anc_side[743:], count)
        count += mul_size
        
        # Cleanup A_2_0
        A_2_0, sqr4 = S(A_2_0+A_2_32, 32)
        
        return A_162_1, ancilla, count


# ==========================================
# GF(233) Implementation
# ==========================================
class Shor233(BaseShor):
    def __init__(self, data_root):
        super().__init__(n=233, block_size=1108, nums2_const=10, data_root=data_root)

    def get_mul_total_size(self):
        return 1341
    
    def get_curve_a(self):
        return 667980384196973083239232503037367571039377840757074144733028669832664

    def inversion_logic(self, a, sqr, ancilla, ancilla2, count):
        sqr1, sqr2, sqr3, sqr4 = sqr[:233], sqr[233:466], sqr[466:699], sqr[699:932]
        anc_main = ancilla
        anc_side = ancilla2
        
        n = 233
        mul_size = 1341
        
        def S(x, p): return self.Square(x, 2*n, p)
        def M(x, y, c_off): 
            return self.mul(x, y, self.Toffoli_qubits[c_off:], 0)
        def C(x, y): self.CNOT_N(x, y)

        # 1 -> 2
        a, A_1_1 = S(a+sqr1, 1)
        A_2_0 = M(a+anc_main[:1108], A_1_1+anc_main[1108:], count)
        count += mul_size
        
        # 2 -> 4
        A_2_0, A_2_2 = S(A_2_0+sqr2, 2)
        a, A_1_1 = S(a+A_1_1, 1) 
        sqr1 = A_1_1
        A_4_0 = M(A_2_2+anc_main[:1108], A_2_0+anc_main[1108:], count)
        count += mul_size

        # 4 -> 8
        A_2_0, A_2_2 = S(A_2_0+A_2_2, 2)
        sqr2 = A_2_2
        A_4_0, A_4_4 = S(A_4_0+sqr1, 4)
        A_8_0 = M(A_4_0+anc_main[:1108], A_4_4+anc_main[1108:], count)
        count += mul_size

        # Chain 1 (8 -> 16)
        A_4_0, A_4_4 = S(A_4_0+A_4_4, 4)
        sqr1 = A_4_4
        A_8_0, A_8_8 = S(A_8_0+sqr2, 8)
        A_16_0 = M(A_8_0+anc_main[:1108], A_8_8+anc_main[1108:], count)
        count += mul_size

        # 16 -> 32
        A_8_0, A_8_8 = S(A_8_0+A_8_8, 8)
        sqr2 = A_8_8
        A_16_0, A_16_16 = S(A_16_0+sqr1, 16)
        A_32_0 = M(A_16_0+anc_main[:1108], A_16_16+anc_main[1108:], count)
        count += mul_size

        # S - Chain 2 (Prep 8 -> 40 side branch)
        C(A_8_0, sqr3)
        sqr3, A_8_32 = S(sqr3+sqr4, 32)
        C(A_8_0, sqr3)

        # Chain 1 (32 -> 64)
        A_16_0, A_16_16 = S(A_16_0+A_16_16, 16)
        sqr1 = A_16_16
        A_32_0, A_32_32 = S(A_32_0+sqr2, 32)
        
        # Merge (32 + 8 = 40)
        C(A_32_0, sqr3)
        # Main
        A_64_0 = M(A_32_0+anc_main[:1108], A_32_32+anc_main[1108:], count)
        count += mul_size
        # Side (A_8_32 * sqr3) -> A_40_0
        A_40_0 = M(A_8_32+anc_side[:1108], sqr3+anc_side[1108:], count)
        count += mul_size
        C(A_32_0, sqr3)

        # S - Chain 2 (Prep 40 -> 104 side branch)
        A_8_0, sqr4 = S(A_8_0+A_8_32, 32) # clear sqr4
        A_40_0, A_40_64 = S(A_40_0+sqr3, 64)
        
        # Chain 1 (64 -> 128)
        A_32_0, A_32_32 = S(A_32_0+A_32_32, 32)
        sqr2 = A_32_32
        A_64_0, A_64_64 = S(A_64_0+sqr1, 64)
        
        # Merge (64 + 40 = 104)
        C(A_64_0, sqr4)
        # Main
        A_128_0 = M(A_64_0+anc_main[:1108], A_64_64+anc_main[1108:], count)
        count += mul_size
        # Side (A_40_64 * sqr4) -> A_104_0
        A_104_0 = M(A_40_64+anc_side[:1108], sqr4+anc_side[1108:], count)
        count += mul_size
        C(A_64_0, sqr4)
        
        # Chain 1 cleanup & Final Merge Prep
        A_64_0, sqr1 = S(A_64_0+A_64_64, 64) # clear sqr1
        A_128_0, A_128_1 = S(A_128_0+sqr2, 1)
        A_104_0, A_104_129 = S(A_104_0+sqr4, 129)

        # Final Mul
        A_232_1 = M(A_104_129+anc_side[:1108], A_128_1+anc_side[1108:], count)
        count += mul_size
        
        A_40_0, sqr3 = S(A_40_0+A_40_64, 64) # clear sqr3
        
        return A_232_1, ancilla, count


# ==========================================
# GF(283) Implementation
# ==========================================
class Shor283(BaseShor):
    def __init__(self, data_root):
        super().__init__(n=283, block_size=1385, nums2_const=11, data_root=data_root)

    def get_mul_total_size(self):
        return 1668

    def inversion_logic(self, a, sqr, ancilla, ancilla2, count):
        sqr1, sqr2, sqr3, sqr4 = sqr[:283], sqr[283:566], sqr[566:849], sqr[849:1132]
        ANC_LIMIT = 1385
        anc_main = ancilla
        anc_side = ancilla2
        
        n = 283
        mul_size = 1668

        def S(x, p): return self.Square(x, 2*n, p)
        def M(x, y, c_off): 
            return self.mul(x, y, self.Toffoli_qubits[c_off:], 0)
        def C(x, y): self.CNOT_N(x, y)

        # 1 -> 2
        a, A_1_1 = S(a+sqr1, 1) # a^2
        A_2_0 = M(a+anc_main[:ANC_LIMIT], A_1_1+anc_main[ANC_LIMIT:], count)
        count += mul_size

        # Chain 1 (2 -> 4)
        A_2_0, A_2_2 = S(A_2_0+sqr2, 2)
        a, A_1_1 = S(a+A_1_1, 1) # update sqr1
        sqr1 = A_1_1
        A_4_0 = M(A_2_2+anc_main[:ANC_LIMIT], A_2_0+anc_main[ANC_LIMIT:], count)
        count += mul_size

        # 4 -> 8
        A_2_0, A_2_2 = S(A_2_0+A_2_2, 2)
        sqr2 = A_2_2 # update sqr2
        A_4_0, A_4_4 = S(A_4_0+sqr1, 4)
        A_8_0 = M(A_4_0+anc_main[:ANC_LIMIT], A_4_4+anc_main[ANC_LIMIT:], count)
        count += mul_size

        # S - Chain 2 (Prep Side Branch)
        C(A_2_0, sqr3)
        sqr3, A_2_8 = S(sqr3+sqr4, 8)
        C(A_2_0, sqr3)

        # Chain 1 (8 -> 16)
        A_4_0, A_4_4 = S(A_4_0+A_4_4, 4)
        A_8_0, A_8_8 = S(A_8_0+sqr2, 8)
        
        # Merge
        C(A_8_0, sqr3)
        # Main Branch
        A_16_0 = M(A_8_0+anc_main[:ANC_LIMIT], A_8_8+anc_main[ANC_LIMIT:], count)
        count += mul_size
        # Side Branch
        A_10_0 = M(A_2_8+anc_side[:ANC_LIMIT], sqr3+anc_side[ANC_LIMIT:], count)
        count += mul_size
        C(A_8_0, sqr3)

        # S - Chain 2
        A_2_0, sqr4 = S(A_2_0+A_2_8, 8) # clear sqr4
        A_10_0, A_10_16 = S(A_10_0+sqr3, 16)

        # Chain 1
        A_8_0, sqr2 = S(A_8_0+A_8_8, 8)
        A_16_0, A_16_16 = S(A_16_0+sqr1, 16)

        # Merge
        C(A_16_0, sqr4)
        # Main Branch (16->32)
        A_32_0 = M(A_16_0+anc_main[:ANC_LIMIT], A_16_16+anc_main[ANC_LIMIT:], count)
        count += mul_size
        # Side Branch (10->26)
        A_26_0 = M(A_10_16+anc_side[:ANC_LIMIT], sqr4+anc_side[ANC_LIMIT:], count)
        count += mul_size
        C(A_16_0, sqr4)

        # Chain 1 (32 -> 64)
        A_16_0, A_16_16 = S(A_16_0+A_16_16, 16)
        sqr1 = A_16_16
        A_32_0, A_32_32 = S(A_32_0+sqr2, 32)
        A_64_0 = M(A_32_0+anc_main[:ANC_LIMIT], A_32_32+anc_main[ANC_LIMIT:], count)
        count += mul_size

        # Chain 1 (64 -> 128)
        A_32_0, A_32_32 = S(A_32_0+A_32_32, 32)
        sqr2 = A_32_32
        A_64_0, A_64_64 = S(A_64_0+sqr1, 64)
        A_128_0 = M(A_64_0+anc_main[:ANC_LIMIT], A_64_64+anc_main[ANC_LIMIT:], count)
        count += mul_size

        # Chain 1 (128 -> 256)
        A_64_0, A_64_64 = S(A_64_0+A_64_64, 64)
        sqr1 = A_64_64
        A_128_0, A_128_128 = S(A_128_0+sqr2, 128)
        A_256_0 = M(A_128_0+anc_main[:ANC_LIMIT], A_128_128+anc_main[ANC_LIMIT:], count)
        count += mul_size

        # Final Prep
        A_128_0, sqr2 = S(A_128_0+A_128_128, 128) # clear sqr2
        A_256_0, A_256_1 = S(A_256_0+sqr1, 1)
        A_26_0, A_26_257 = S(A_26_0+sqr4, 257)

        # Final Mul (Combine 256 + 26 = 282)
        # Note: Logic file says recursive_karatsuba(A_26_257, A_256_1, ..., ancilla2)
        # Mapping ancilla2(logic) -> ancilla3(circuit) -> anc_side
        A_282_1 = M(A_26_257+anc_side[:ANC_LIMIT], A_256_1+anc_side[ANC_LIMIT:], count)
        count += mul_size

        # Cleanup sqr3
        A_10_0, sqr3 = S(A_10_0+A_10_16, 16)

        return A_282_1, ancilla, count


# ==========================================
# GF(571) Implementation
# ==========================================
class Shor571(BaseShor):
    def __init__(self, data_root):
        super().__init__(n=571, block_size=2998, nums2_const=13, data_root=data_root)

    def get_mul_total_size(self):
        return 3569
        
    def inversion_logic(self, a, sqr, ancilla, ancilla2, count):

        sqr1, sqr2, sqr3, sqr4 = sqr[:571], sqr[571:1142], sqr[1142:1713], sqr[1713:2284]
        ANC_LIMIT = 2998  
        anc_main = ancilla
        anc_side = ancilla2
        
        n = 571
        mul_size = 3569

        def S(x, p): return self.Square(x, 2*n, p)
        def M(x, y, c_off): 
            return self.mul(x, y, self.Toffoli_qubits[c_off:], 0)
        def C(x, y): self.CNOT_N(x, y)

        # 1 -> 2
        a, A_1_1 = S(a+sqr1, 1) 
        A_2_0 = M(a+anc_main[:ANC_LIMIT], A_1_1+anc_main[ANC_LIMIT:], count)
        count += mul_size

        # Chain 1 (2 -> 4)
        A_2_0, A_2_2 = S(A_2_0+sqr2, 2)
        a, sqr1 = S(a+A_1_1, 1) 
        A_4_0 = M(A_2_2+anc_main[:ANC_LIMIT], A_2_0+anc_main[ANC_LIMIT:], count)
        count += mul_size

        # 4 -> 8
        A_2_0, sqr2 = S(A_2_0+A_2_2, 2)
        A_4_0, A_4_4 = S(A_4_0+sqr1, 4)
        A_8_0 = M(A_4_0+anc_main[:ANC_LIMIT], A_4_4+anc_main[ANC_LIMIT:], count)
        count += mul_size

        # S - Chain 2 (Prep Side)
        C(A_2_0, sqr3)
        sqr3, A_2_8 = S(sqr3+sqr4, 8)
        C(A_2_0, sqr3)

        # Chain 1 (8 -> 16)
        A_4_0, A_4_4 = S(A_4_0+A_4_4, 4)
        A_8_0, A_8_8 = S(A_8_0+sqr2, 8)

        # Merge
        C(A_8_0, sqr3)
        # Main
        A_16_0 = M(A_8_0+anc_main[:ANC_LIMIT], A_8_8+anc_main[ANC_LIMIT:], count)
        count += mul_size
        # Side
        A_10_0 = M(A_2_8+anc_side[:ANC_LIMIT], sqr3+anc_side[ANC_LIMIT:], count)
        count += mul_size
        C(A_8_0, sqr3)

        # S - Chain 2
        A_2_0, sqr4 = S(A_2_0+A_2_8, 8) # clear sqr4
        A_10_0, A_10_16 = S(A_10_0+sqr3, 16)

        # Chain 1
        A_8_0, A_8_8 = S(A_8_0+A_8_8, 8)
        sqr2 = A_8_8
        A_16_0, A_16_16 = S(A_16_0+sqr1, 16)

        # Merge
        C(A_16_0, sqr4)
        # Main
        A_32_0 = M(A_16_0+anc_main[:ANC_LIMIT], A_16_16+anc_main[ANC_LIMIT:], count)
        count += mul_size
        # Side
        A_26_0 = M(A_10_16+anc_side[:ANC_LIMIT], sqr4+anc_side[ANC_LIMIT:], count)
        count += mul_size
        C(A_16_0, sqr4)

        # S - Chain 2
        A_26_0, A_26_32 = S(A_26_0+sqr4, 32)
        A_10_0, sqr3 = S(A_10_0+A_10_16, 16) # clear sqr3

        # Chain 1 (32 -> 64)
        A_16_0, sqr1 = S(A_16_0+A_16_16, 16)
        A_32_0, A_32_32 = S(A_32_0+sqr2, 32)

        # Merge
        C(A_32_0, sqr3)
        # Main
        A_64_0 = M(A_32_0+anc_main[:ANC_LIMIT], A_32_32+anc_main[ANC_LIMIT:], count)
        count += mul_size
        # Side (Note: A_26_32 + sqr3)
        A_58_0 = M(A_26_32+anc_side[:ANC_LIMIT], sqr3+anc_side[ANC_LIMIT:], count)
        count += mul_size
        C(A_32_0, sqr3)

        # Chain 1 (64 -> 128)
        A_32_0, sqr2 = S(A_32_0+A_32_32, 32)
        A_64_0, A_64_64 = S(A_64_0+sqr1, 64)
        A_128_0 = M(A_64_0+anc_main[:ANC_LIMIT], A_64_64+anc_main[ANC_LIMIT:], count)
        count += mul_size

        # Chain 1 (128 -> 256)
        A_64_0, sqr1 = S(A_64_0+A_64_64, 64)
        A_128_0, A_128_128 = S(A_128_0+sqr2, 128)
        A_256_0 = M(A_128_0+anc_main[:ANC_LIMIT], A_128_128+anc_main[ANC_LIMIT:], count)
        count += mul_size

        # Chain 1 (256 -> 512)
        A_128_0, sqr2 = S(A_128_0+A_128_128, 128)
        A_256_0, A_256_256 = S(A_256_0+sqr1, 256)
        A_512_0 = M(A_256_0+anc_main[:ANC_LIMIT], A_256_256+anc_main[ANC_LIMIT:], count)
        count += mul_size

        # Final Prep
        A_256_0, sqr1 = S(A_256_0+A_256_256, 256) # clear sqr1
        A_512_0, A_512_1 = S(A_512_0+sqr2, 1)
        A_58_0, A_58_513 = S(A_58_0+sqr3, 513)

        # Final Mul (Combine 58 + 512 = 570)
        A_570_1 = M(A_58_513+anc_side[:ANC_LIMIT], A_512_1+anc_side[ANC_LIMIT:], count)
        count += mul_size

        # Cleanup sqr4
        A_26_0, sqr4 = S(A_26_0+A_26_32, 32)

        return A_570_1, ancilla, count