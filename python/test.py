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