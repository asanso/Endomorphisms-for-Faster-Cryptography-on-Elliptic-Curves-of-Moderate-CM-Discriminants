from sage.all import ZZ, FiniteField
class WeierstrassCurve():
    def __init__(self, p, a, b, r, cofactor):
        self.p = p
        self.Fp = FiniteField(p)
        self.a = self.Fp(a)
        self.b = self.Fp(b)
        self.r = r
        self.cofactor = cofactor

    def __str__(self):
        a = ZZ(self.a)
        b = ZZ(self.b)
        p = ZZ(self.p)
        return f"Elliptic curve in Weierstrass form y² = x³ + {a}x + {b} over GF({p})"

    def random_point(self):
        x = self.Fp.random_element()
        rhs = x**3 + self.a * x + self.b
        while not rhs.is_square():
            x = self.Fp.random_element()
            rhs = x**3 + self.a * x + self.b
        y = rhs.sqrt()
        return PointWeierstrass(x, y, self)

    def point_of_order_r(self):
        P = self.random_point().clear_cofactor()
        while P.is_zero():
            P = self.random_point().clear_cofactor()
        assert P.scalar_mul(self.r).is_zero()
        return P

class PointWeierstrass():
    def __init__(self, x, y, curve):
        self.x = x
        self.y = y
        self.curve = curve

    def __str__(self):
        return 'x:' + str(self.x) + '\n' + \
            'y:' + str(self.y)

    def __eq__(self, other):
        if self.is_zero() and other.is_zero():
            return True  # Both are the point at infinity
        if self.is_zero() or other.is_zero():
            return False  # One is infinity, the other is not
        return self.x == other.x and self.y == other.y

    def is_zero(self):
        return self.y.is_zero() and self.x.is_zero()

    def neg(self):
        return PointWeierstrass(self.x, -self.y, self.curve)

    def on_curve(self):
        x, y = self.x, self.y
        a, b = self.curve.a, self.curve.b
        return y**2 == x**3 + a*x + b

    def add(self, other):
        if self.is_zero():
            return other
        if other.is_zero():
            return self
        if self.x == other.x and self.y != other.y:
            return PointWeierstrass(self.curve.Fp(0), self.curve.Fp(0), self.curve)  # Point at infinity

        if self == other:
            return self.double()

        x1, y1 = self.x, self.y
        x2, y2 = other.x, other.y
        p = self.curve.p

        # Compute the slope
        m = (y2 - y1) / (x2 - x1)
        x3 = m**2 - x1 - x2
        y3 = m * (x1 - x3) - y1

        return PointWeierstrass(x3, y3, self.curve)

    def double(self):
        if self.is_zero():
            return self

        x, y = self.x, self.y
        a = self.curve.a

        # Compute the slope for doubling
        m = (3 * x**2 + a) / (2 * y)
        x3 = m**2 - 2 * x
        y3 = m * (x - x3) - y

        return PointWeierstrass(x3, y3, self.curve)

    def scalar_mul(self, n):
        if n<0:
            n = -n
            P = self.neg()
        else:
            P = self
        R = P
        for b in ZZ(n).bits()[-2::-1]:
            R = R.double()
            if b == 1:
                R = R.add(P)
        return R

    def multi_scalar_mul(self, k1, other, k2):
        P = self
        if k1<0:
            k1=-k1
            P = P.neg()
        if k2<0:
            k2=-k2
            other = other.neg()
        PplusOther = P.add(other)
        bits_k1 = ZZ(k1).bits()
        bits_k2 = ZZ(k2).bits()
        while len(bits_k1) < len(bits_k2):
            bits_k1.append(0)
        while len(bits_k2) < len(bits_k1):
            bits_k2.append(0)
        R = PointWeierstrass(self.curve.Fp(0), self.curve.Fp(0), self.curve)
        for i in range(len(bits_k1)-1,-1,-1):
            R = R.double()
            if bits_k1[i] == 1 and bits_k2[i] == 0:
                R = R.add(self)
            if bits_k1[i] == 0 and bits_k2[i] == 1:
                R = R.add(other)
            if bits_k1[i] == 1 and bits_k2[i] == 1:
                R = R.add(PplusOther)
        return R

    def clear_cofactor(self):
        return self.scalar_mul(self.curve.cofactor)
