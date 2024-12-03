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
        z = self.Fp.random_element()
        x *= z
        y *= z
        return PointWeierstrass(x, y, z,self)

    def point_of_order_r(self):
        P = self.random_point().clear_cofactor()
        while P.is_zero():
            P = self.random_point().clear_cofactor()
        assert P.scalar_mul(self.r).is_zero()
        return P

class PointWeierstrass():
    def __init__(self, X, Y, Z, curve):
        self.X = X
        self.Y = Y
        self.Z = Z
        self.curve = curve

    def __str__(self):
        if self.is_zero():
            return "Point at infinity"
        x = self.X 
        y = self.Y
        z = self.Z
        return f"x:{x}\ny:{y}\nz:{z}"

    def __eq__(self, other):
        if self.is_zero() and other.is_zero():
            return True  # Both are the point at infinity
        if self.is_zero() or other.is_zero():
            return False  # One is infinity, the other is not
        return self.X * other.Z == other.X * self.Z and self.Y * other.Z == other.Y * self.Z

    def is_zero(self):
        return self.Z.is_zero()

    def neg(self):
        return PointWeierstrass(self.X, -self.Y, self.Z, self.curve)

    def on_curve(self):
        X, Y, Z = self.X, self.Y, self.Z
        a, b = self.curve.a, self.curve.b
        return Y**2 * Z == X**3 + a * X * Z**2 + b * Z**3

    def add(self, other):
        if self.is_zero():
            return other
        if other.is_zero():
            return self
        if self.X * other.Z == other.X * self.Z and self.Y * other.Z == -other.Y * self.Z:
            return PointWeierstrass(self.curve.Fp(0), self.curve.Fp(1), self.curve.Fp(0), self.curve)  # Point at infinity

        X1, Y1, Z1 = self.X, self.Y, self.Z
        X2, Y2, Z2 = other.X, other.Y, other.Z

        # Projective point addition formulas
        U1 = Y2 * Z1
        U2 = Y1 * Z2
        V1 = X2 * Z1
        V2 = X1 * Z2

        if V1 == V2:
            if U1 == U2:
                return self.double()  # Same point
            else:
                return PointWeierstrass(self.curve.Fp(0), self.curve.Fp(1), self.curve.Fp(0), self.curve)  # Point at infinity

        U = U1 - U2
        V = V1 - V2
        W = Z1 * Z2
        A = U**2 * W - V**3 - 2 * V**2 * V2
        X3 = V * A
        Y3 = U * (V**2 * V2 - A) - V**3 * U2
        Z3 = V**3 * W

        return PointWeierstrass(X3, Y3, Z3, self.curve)

    def double(self):
        if self.is_zero():
            return self

        X, Y, Z = self.X, self.Y, self.Z
        a = self.curve.a

        # Projective point doubling formulas
        W = 3 * X**2 + a * Z**2
        S = Y * Z
        B = X * Y * S
        H = W**2 - 8 * B
        X3 = 2 * H * S
        Y3 = W * (4 * B - H) - 8 * Y**2 * S**2
        Z3 = 8 * S**3

        return PointWeierstrass(X3, Y3, Z3, self.curve)

    def scalar_mul(self, n):
        if n < 0:
            n = -n
            P = self.neg()
        else:
            P = self
        R = PointWeierstrass(self.curve.Fp(0), self.curve.Fp(1), self.curve.Fp(0), self.curve)  # Point at infinity
        for b in ZZ(n).bits()[-2::-1]:
            R = R.double()
            if b == 1:
                R = R.add(P)
        return R

    def clear_cofactor(self):
        return self.scalar_mul(self.curve.cofactor)
