import random
from curve import WeierstrassCurve
from lollipop489201 import Lollipop489201
from lollipop574261 import Lollipop574261
from gc256c import GC256C

# Lollipop489201

p = 1910157204347957325700187962480217512925138482090399484362397
a = 73275333332267847499581501376863252276520692179021512625126
b = 1538008641579097707704221968675032141849999412179326013460607
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

# clear cofactor
assert P.clear_cofactor().scalar_mul(E.r).is_zero()

# multi scalar mul
for [k1,k2] in [[12,13], [-12,13], [12,-13], [-12,-13]]:
    P1 = E.random_point()
    P2 = E.random_point()
    assert P1.scalar_mul(k1).add(P2.scalar_mul(k2)) == P1.multi_scalar_mul(k1,P2,k2)

E = Lollipop489201(p, a, b, r, 1, eigen)
P = E.random_point()
R1 = P.scalar_mul(n)
R2 = P.fast_scalar_mul(n)
assert R1 == R2

# Lollipop574261

p = 2163160611951109656514578155686584121369935567340260530653195705987773102423967
a = 1600224790493789893768878688270410785959510786167424698577648903928123990501736
b = 442173884976592923728839681693322529291546262090176625621658004289727200834915
r = 2163160611951109656514578155686584121370430948758200237999126133244251338506461
eigen = 1158066151478108453899950316130401430242909130290171218593900211556999875716624

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

# clear cofactor
assert P.clear_cofactor().scalar_mul(E.r).is_zero()

# multi scalar mul
for [k1,k2] in [[12,13], [-12,13], [12,-13], [-12,-13]]:
    P1 = E.random_point()
    P2 = E.random_point()
    assert P1.scalar_mul(k1).add(P2.scalar_mul(k2)) == P1.multi_scalar_mul(k1,P2,k2)

E =  Lollipop574261(p, a, b, r, 1, eigen)
P = E.random_point()
R1 = P.scalar_mul(n)
R2 = P.fast_scalar_mul(n)
#assert R1 == R2

#GC256C

p = 2**255 + 3225
a =  p -3
b = 28091019353058090096996979000309560759124368558014865957655842872397301267595
r = 57896044618658097711785492504343953927102133160255826820068844496087732066703
eigen = 25745384436187293773257686313148572361251961181400080445531962069440031481800


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

# clear cofactor
assert P.clear_cofactor().scalar_mul(E.r).is_zero()

# multi scalar mul
for [k1,k2] in [[12,13], [-12,13], [12,-13], [-12,-13]]:
    P1 = E.random_point()
    P2 = E.random_point()
    assert P1.scalar_mul(k1).add(P2.scalar_mul(k2)) == P1.multi_scalar_mul(k1,P2,k2)

E = GC256C(p, a, b, r, 1, eigen)
P = E.random_point()
R1 = P.scalar_mul(n)
R2 = P.fast_scalar_mul(n)
assert R1 == R2
