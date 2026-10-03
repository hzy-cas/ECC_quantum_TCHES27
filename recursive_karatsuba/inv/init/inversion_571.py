from .basic_gates_571 import (
    Square, 
    recursive_karatsuba, 
    Reduction, 
    CNOT_n,
    GF_N
)

def Inverison_Itoh_Tsujii_based(a, n, sqr1, sqr2, sqr3, sqr4, count, ancilla1, ancilla2):
    n = GF_N
    a, A_1_1 = Square(a+sqr1, n, 1) # a^2
    count = 0
    A_2_0, count, ancilla1 = recursive_karatsuba(a, A_1_1 , n, count, ancilla1) 
    A_2_0 = Reduction(A_2_0)


    # chain1
    A_2_0, A_2_2 = Square(A_2_0+sqr2, n, 2)
    a, sqr1 = Square(a+A_1_1, n, 1) 
    count = 0
    A_4_0, count, ancilla1 = recursive_karatsuba(A_2_2, A_2_0, n, count, ancilla1) 
    A_4_0 = Reduction(A_4_0)

    A_2_0, sqr2 = Square(A_2_0+A_2_2, n, 2)
    A_4_0, A_4_4 = Square(A_4_0+sqr1, n, 4)
    count = 0
    A_8_0, count, ancilla1 = recursive_karatsuba(A_4_0, A_4_4, n, count, ancilla1) 
    A_8_0 = Reduction(A_8_0)
    
    # S
    # chain2
    CNOT_n(A_2_0, sqr3)
    sqr3, A_2_8 = Square(sqr3+sqr4, n, 8)
    CNOT_n(A_2_0, sqr3)
    # chain1
    A_4_0, A_4_4 = Square(A_4_0+A_4_4, n, 4)
    A_8_0, A_8_8 = Square(A_8_0+sqr2, n, 8)
    #M
    CNOT_n(A_8_0, sqr3)
    count = 0
    A_16_0, count, ancilla1 = recursive_karatsuba(A_8_0, A_8_8, n, count, ancilla1) 
    A_16_0 = Reduction(A_16_0)
    count = 0
    A_10_0, count, ancilla2 = recursive_karatsuba(A_2_8, sqr3, n, count, ancilla2) 
    A_10_0 = Reduction(A_10_0)
    CNOT_n(A_8_0, sqr3)


    # S
    # chain2
    A_2_0,sqr4 = Square(A_2_0+A_2_8, n, 8) # clear sqr4
    A_10_0, A_10_16 = Square(A_10_0+sqr3, n, 16)
    # chain1    
    A_8_0, A_8_8 = Square(A_8_0+A_8_8, n, 8)
    sqr2 = A_8_8
    A_16_0, A_16_16 = Square(A_16_0+sqr1, n, 16)
    #M
    CNOT_n(A_16_0, sqr4)
    count = 0
    A_32_0, count, ancilla1 = recursive_karatsuba(A_16_0, A_16_16, n, count, ancilla1) 
    A_32_0 = Reduction(A_32_0)
    count = 0
    A_26_0, count, ancilla2 = recursive_karatsuba(A_10_16, sqr4, n, count, ancilla2) 
    A_26_0 = Reduction(A_26_0)
    CNOT_n(A_16_0, sqr4)

    # S
    # chain2
    A_26_0,A_26_32= Square(A_26_0+sqr4, n, 32) 
    A_10_0, sqr3 = Square(A_10_0+A_10_16, n, 16)# clear sqr3
    # chain1
    A_16_0, sqr1 = Square(A_16_0+A_16_16, n, 16)
    A_32_0, A_32_32 = Square(A_32_0+sqr2, n, 32)
    #M
    CNOT_n(A_32_0, sqr3)
    count = 0
    A_64_0, count, ancilla1 = recursive_karatsuba(A_32_0, A_32_32, n, count, ancilla1) 
    A_64_0 = Reduction(A_64_0)
    count = 0
    A_58_0, count, ancilla2 = recursive_karatsuba(A_26_32, sqr3, n, count, ancilla2) 
    A_58_0 = Reduction(A_58_0)
    CNOT_n(A_32_0, sqr3)

    A_32_0, sqr2 = Square(A_32_0+A_32_32, n, 32)
    A_64_0, A_64_64 = Square(A_64_0+sqr1, n, 64)
    count = 0
    A_128_0, count, ancilla1 = recursive_karatsuba(A_64_0, A_64_64, n, count, ancilla1)  
    A_128_0 = Reduction(A_128_0)

    A_64_0, sqr1 = Square(A_64_0+A_64_64, n, 64)
    A_128_0, A_128_128 = Square(A_128_0+sqr2, n, 128)
    count = 0
    A_256_0, count, ancilla1 = recursive_karatsuba(A_128_0, A_128_128, n, count, ancilla1)  
    A_256_0 = Reduction(A_256_0)

    A_128_0, sqr2 = Square(A_128_0+A_128_128, n, 128)
    A_256_0, A_256_256 = Square(A_256_0+sqr1, n, 256)
    count = 0
    A_512_0, count, ancilla1 = recursive_karatsuba(A_256_0, A_256_256, n, count, ancilla1)  
    A_512_0 = Reduction(A_512_0)


    A_256_0, sqr1 = Square(A_256_0+A_256_256, n, 256)# clear sqr1
    A_512_0, A_512_1 = Square(A_512_0+sqr2, n, 1)
    A_58_0, A_58_513 = Square(A_58_0+sqr3, n, 513)

    count = 0
    A_570_1, count, ancilla2 = recursive_karatsuba(A_58_513, A_512_1, n, count, ancilla2) 
    A_570_1 = Reduction(A_570_1)

    A_26_0,sqr4= Square(A_26_0+A_26_32, n, 32) 

    ## sqr4 sqr1 clear
    return A_570_1, ancilla1, ancilla2