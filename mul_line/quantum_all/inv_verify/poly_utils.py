class GF2n:
    def __init__(self, n, irreducible_poly_degrees):
        self.n = n
        self.modulus = 0
        for deg in irreducible_poly_degrees:
            self.modulus |= (1 << deg)

    def multiply(self, a, b):
        p = 0
        while b > 0:
            if b & 1:
                p ^= a
            a <<= 1
            if a & (1 << self.n):
                a ^= self.modulus
            b >>= 1
        return p & ((1 << self.n) - 1)

    def square(self, a):
        return self.multiply(a, a)

    def power(self, a, exponent):
        res = 1
        base = a
        while exponent > 0:
            if exponent & 1:
                res = self.multiply(res, base)
            base = self.square(base)
            exponent >>= 1
        return res

    def inverse(self, a):
        return self.power(a, (1 << self.n) - 2)