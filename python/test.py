import random
from curve import WeierstrassCurve
from lollipop489201 import Lollipop489201, Lollipop489201Point

p = 1910157204347957325700187962480217512925138482090399484362397
a=73275333332267847499581501376863252276520692179021512625126;
b=1538008641579097707704221968675032141849999412179326013460607;
r = 1910157204347957325700187962477453761504460135807772883063193
eigen = 500517855408530850485984845457693413895246088040882972578379

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
n = random.randint(0,r)
R = P.scalar_mul(n)
assert R.on_curve()
P.scalar_mul(E.cofactor * E.r).is_zero()