import random
from curve import WeierstrassCurve

p = 1910157204347957325700187962480217512925138482090399484362397
a=73275333332267847499581501376863252276520692179021512625126;
b=1538008641579097707704221968675032141849999412179326013460607;
r = 1910157204347957325700187962477453761504460135807772883063193
E = WeierstrassCurve(p, a, b, r, 1)

P = E.random_point()
Q = E.random_point()

# add
R = P.add(Q)
assert R.on_curve()
assert R == Q.add(P)

# double
R = P.double()
assert R.on_curve()
assert R  == P.add(P)


# scalar mul
n = 5
R = P.scalar_mul(n)
assert R.on_curve()
P.scalar_mul(E.cofactor * E.r).is_zero()


# multi scalar mul
#for [k1,k2] in [[12,13], [-12,13], [12,-13], [-12,-13]]:
    #P1 = E.random_point()
    #P2 = E.random_point()
    #assert P1.scalar_mul(k1).add(P2.scalar_mul(k2)) == P1.multi_scalar_mul(k1,P2,k2)
