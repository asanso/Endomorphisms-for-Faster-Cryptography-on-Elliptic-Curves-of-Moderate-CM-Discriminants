p = 1910157204347957325700187962480217512925138482090399484362397
aboldhat=73275333332267847499581501376863252276520692179021512625126;
bboldhat=1538008641579097707704221968675032141849999412179326013460607;
E0 = EllipticCurve(GF(p), [aboldhat,bboldhat])
t =E0.trace_of_frobenius()
disc = t^2 - 4*p
print(factor(disc))
print()
#list of curves on the crater
Fp = GF(p)
FpX = Fp['X']
X = FpX.gen()
D = fundamental_discriminant(-547)
H = FpX(hilbert_class_polynomial(D))
print(H.roots())
print()

# computing the cycle

assert E0.j_invariant() == H.roots()[2][0]
phi0 =  E0.isogenies_prime_degree(11)[4]
E1 = phi0.codomain()
assert E1.j_invariant() == H.roots()[1][0]

phi1 = E1.isogenies_prime_degree(11)[0]
E2 = phi1.codomain()
assert E2.j_invariant() == H.roots()[0][0]

phi2 = E2.isogenies_prime_degree(11)[4]
assert phi2.codomain().j_invariant() == H.roots()[2][0]