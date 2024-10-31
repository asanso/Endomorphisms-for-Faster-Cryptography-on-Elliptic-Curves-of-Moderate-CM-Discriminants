p = 2^255 +3225
a =  p -3
b = 28091019353058090096996979000309560759124368558014865957655842872397301267595
E0 = EllipticCurve(GF(p), [a,b])
t = E.trace_of_frobenius()
disc = t^2 - 4*p
factor(disc)
print()
#list of curves on the crater
Fp = GF(p)
FpX = Fp['X']
X = FpX.gen()
D = fundamental_discriminant(-619)
H = FpX(hilbert_class_polynomial(D))
print(H.roots())
print()

# computing the cycle
assert E0.j_invariant() == H.roots()[2][0]
phi0 =  E0.isogenies_prime_degree(11)[4]
E1 = phi0.codomain()
assert E1.j_invariant() == H.roots()[1][0]