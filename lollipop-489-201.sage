import timeit

def end_composition(P):
    x0 = P[0]
    y0 = P[1]
    x1 = x_num0(x0)
    y1 = y_num0(x0,y0)
    x2 = x_num1(x1)
    y2 = y_num1(x1,y1)
    x3 = x_num2(x2)
    y3 = y_num2(x2,y2)
    x_fin = x_num_iso(x3)
    y_fin = y_num_iso(x3,y3)

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

# Computing eigenvalue 
end =(phi2*phi1*phi0)
iso =end.codomain().isomorphism_to(E0)
full_end = (iso*end)
trace = full_end.trace()
norm = 11^3
r = E0.order()
Fr = GF(r)
R.<x> = PolynomialRing(Fr)
poly = x^2 - trace*x + norm
roots = poly.roots()
P = E0.random_point()
Q = full_end(P)
eigen = roots[1][0]
 
assert Q == eigen*P

# endomorphism rational maps

x_num0 = phi0.x_rational_map().numerator()
x_num1 = phi1.x_rational_map().numerator()
x_num2 = phi2.x_rational_map().numerator()
x_num_iso = iso.x_rational_map().numerator()

x_denom0 = phi0.x_rational_map().denominator()
x_denom1 = phi1.x_rational_map().denominator()
x_denom2 = phi2.x_rational_map().denominator()

y_num0 = phi0.rational_maps()[1].numerator()
y_num1 = phi1.rational_maps()[1].numerator()
y_num2 = phi2.rational_maps()[1].numerator()
y_num_iso = iso.rational_maps()[1].numerator()

y_denom0 = phi0.rational_maps()[1].denominator()
y_denom1 = phi1.rational_maps()[1].denominator()
y_denom2 = phi2.rational_maps()[1].denominator()

# GLV

M = Matrix([[int(-eigen),1], [int(r),0]])
#print(M)
N = M.LLL()
N_inv = N**-1

n = ZZ.random_element(r)
#print(len(n.str(2)))
S1 = n*P
S2 = fast_scalar_mul(n,P)
assert S1 == S2