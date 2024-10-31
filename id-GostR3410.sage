p = 2^255 +3225
a =  p -3
b = 28091019353058090096996979000309560759124368558014865957655842872397301267595
E0 = EllipticCurve(GF(p), [a,b])
t = E0.trace_of_frobenius()
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
assert E0.j_invariant() == H.roots()[4][0]
phi0 =  E0.isogenies_prime_degree(5)[0]
E1 = phi0.codomain()
assert E1.j_invariant() == H.roots()[3][0]

phi1 = E1.isogenies_prime_degree(5)[1]
E2 = phi1.codomain()
assert E2.j_invariant() == H.roots()[1][0]

phi2 = E2.isogenies_prime_degree(5)[1]
E3 = phi2.codomain()
assert E3.j_invariant() == H.roots()[2][0]

phi3 = E3.isogenies_prime_degree(5)[1]
E4 = phi3.codomain()
assert E4.j_invariant() == H.roots()[0][0]

phi4 = E4.isogenies_prime_degree(5)[0]
assert phi4.codomain().j_invariant() == H.roots()[4][0]