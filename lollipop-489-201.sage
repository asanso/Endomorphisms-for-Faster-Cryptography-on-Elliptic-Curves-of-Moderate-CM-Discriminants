import timeit

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


def projective_maps(phi,Fp):
    rX,sXY = phi
    Fpx = Fp['x']
    x = Fpx.gen()
    FpX = Fpx.fraction_field()
    X = FpX.gen()
    Fpxz = Fp['x', 'z']
    FpXZ = FractionField(Fpxz)
    X, Z = FpXZ.gens()
    sX = sXY(y=1)
    rXZ = rX(x=X/Z)
    sXZ = sX(x=X/Z)
    # (x:y:z) -> (z * a(x,z) : y * b(x,z) : z * c(x,z))
    a = rXZ.numerator() * sXZ.denominator()
    b = sXZ.numerator() * rXZ.denominator()
    c = rXZ.denominator() * sXZ.denominator()
    return a,b,c


def end_composition(P):
    x0 = P[0]
    y0 = P[1]
    z0 = 1 

    #1st isogeny
    x1 = z0*a0(x0,z0)  
    y1 = y0 *b0(x0,z0)
    z1 = z0* c0(x0,z0)
    #2nd isogeny
    x2 = z1*a1(x1,z1)  
    y2 = y1 *b1(x1,z1)
    z2 = z1* c1(x1,z1)
    #3rd isogeny
    x3 = z2*a2(x2,z2)  
    y3 = y2 *b2(x2,z2)
    z3 = z2* c2(x2,z2)
    return  isoX(x3,1), isoY(y3,1),z3
    

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

end =(phi2*phi1*phi0)
iso =end.codomain().isomorphism_to(E0)
full_end = (iso*end)

# Computing eigenvalue 
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
a0,b0,c0 = projective_maps(phi0,Fp)
a1,b1,c1 = projective_maps(phi1,Fp)
a2,b2,c2 = projective_maps(phi2,Fp)
isoX = iso.rational_maps()[0]
isoY = iso.rational_maps()[1]

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