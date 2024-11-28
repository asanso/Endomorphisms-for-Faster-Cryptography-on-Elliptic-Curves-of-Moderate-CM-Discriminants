from sage.all import ZZ, Matrix, FiniteField, vector
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
            return PointWeierstrass(0, 0, self.curve)  # Point at infinity

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
        if n < 0:
            n = -n
            P = self.neg()
        else:
            P = self
        R = PointWeierstrass(self.curve.Fp(0), self.curve.Fp(0), self.curve)  # Point at infinity
        for b in ZZ(n).bits():
            R = R.double()
            if b == 1:
                R = R.add(P)
        return R

    def clear_cofactor(self):
        return self.scalar_mul(self.curve.cofactor)
