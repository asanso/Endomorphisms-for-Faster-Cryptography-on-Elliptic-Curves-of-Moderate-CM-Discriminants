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

# Computing eigenvalue 
end =(phi4*phi3*phi2*phi1*phi0)
iso =end.codomain().isomorphism_to(E0)
full_end = (iso*end)
trace = full_end.trace()

norm = 5^5
r = E0.order()
Fr = GF(r)
R.<x> = PolynomialRing(Fr)
poly = x^2 - trace*x + norm
roots = poly.roots()
P = E0.random_point()
Q = full_end(P)
eigen = roots[1][0]

assert Q == eigen*P

# GLV

M = Matrix([[int(-eigen),1], [int(r),0]])
#print(M)
N = M.LLL()
N_inv = N**-1

def multi_scalar_mul(P, k1, endo, k2):
    return k1*P + k2*endo(P)

def fast_scalar_mul(n,P):
    beta = vector([n,0])*N_inv
    b = vector([int(beta[0]), int(beta[1])]) * N
    k1 = n-b[0]
    k2 = -b[1]
    #print(len(k1.str(2)))
    #print(len(k2.str(2)))
    return  multi_scalar_mul(P,k1, full_end, k2)

n = ZZ.random_element(r)
#print(len(n.str(2)))
S1 = n*P
S2 = fast_scalar_mul(n,P)
assert S1 == S2
