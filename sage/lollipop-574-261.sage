def multi_scalar_mul(P, k1, endo, k2):
    return k1*P + k2*endo(P)

def fast_scalar_mul(n,P):
    beta = vector([n,0])*N_inv
    b = vector([int(beta[0]), int(beta[1])]) * N
    k1 = n-b[0]
    k2 = -b[1]
    return  multi_scalar_mul(P,k1, full_end, k2)

def projective_maps_optimized(phi,Fp):
    rX,sXY = phi
    Fpx = Fp['x']
    x = Fpx.gen()
    FpX = Fpx.fraction_field()
    X = FpX.gen()
    Fpxz = Fp['x', 'z']
    FpXZ = FractionField(Fpxz)
    X, Z = FpXZ.gens()
    
    psi1 = rX.numerator()
    psi3 = -rX.denominator().sqrt()
    sX = sXY(y=1)
    assert psi3^3 == (sX.denominator()/2176782336)
    psi2 = sX.numerator()/2176782336
    psi1XZ = psi1(x=X/Z)
    psi2XZ = psi2(x=X/Z)
    psi3XZ = psi3(x=X/Z)
    a = psi1XZ*psi3XZ *Z^7*2176782336
    b = psi2XZ *Z^6*2176782336
    c = psi3XZ^3 *Z^6*2176782336
    return a,b,c



p = 2163160611951109656514578155686584121369935567340260530653195705987773102423967
aboldhat=1600224790493789893768878688270410785959510786167424698577648903928123990501736
bboldhat=442173884976592923728839681693322529291546262090176625621658004289727200834915
E0 = EllipticCurve(GF(p), [aboldhat,bboldhat])
t =E0.trace_of_frobenius()
disc = t^2 - 4*p
print(factor(disc))
print()
#list of curves on the crater
Fp = GF(p)
FpX = Fp['X']
X = FpX.gen()
D = fundamental_discriminant(-3019)
H = FpX(hilbert_class_polynomial(D))
print(H.roots())
print()

# computing the cycle
assert E0.j_invariant() == H.roots()[6][0]
phi0 =  E0.isogenies_prime_degree(5)[0]
E1 = phi0.codomain()
assert E1.j_invariant() == H.roots()[3][0]

phi1 = E1.isogenies_prime_degree(5)[0]
E2 = phi1.codomain()
assert E2.j_invariant() == H.roots()[0][0]

phi2 = E2.isogenies_prime_degree(5)[1]
E3 = phi2.codomain()
assert E3.j_invariant() == H.roots()[5][0]

phi3 = E3.isogenies_prime_degree(5)[0]
E4 = phi3.codomain()
assert E4.j_invariant() == H.roots()[2][0]

phi4 = E4.isogenies_prime_degree(5)[1]
E5 = phi4.codomain()
assert E5.j_invariant() == H.roots()[1][0]

phi5 = E5.isogenies_prime_degree(5)[0]
E6 = phi5.codomain()
assert E6.j_invariant() == H.roots()[4][0]

phi6 = E6.isogenies_prime_degree(5)[1]
assert phi6.codomain().j_invariant() == H.roots()[6][0]

# Computing eigenvalue 
end =(phi6*phi5*phi4*phi3*phi2*phi1*phi0)
iso =end.codomain().isomorphism_to(E0)
full_end = (iso*end)
trace = full_end.trace()

norm = 5^7
r = E0.order()
Fr = GF(r)
R.<x> = PolynomialRing(Fr)
poly = x^2 - trace*x + norm
roots = poly.roots()
P = E0.random_point()
Q = full_end(P)
eigen = roots[0][0]

assert Q == eigen*P

# GLV

M = Matrix([[int(-eigen),1], [int(r),0]])
#print(M)
N = M.LLL()
N_inv = N**-1

n = ZZ.random_element(r)
S1 = n*P
S2 = fast_scalar_mul(n,P)
assert S1 == S2

a0,b0,c0 = projective_maps_optimized(phi0,Fp)