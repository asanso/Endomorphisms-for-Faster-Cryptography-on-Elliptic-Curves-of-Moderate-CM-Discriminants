import random
from sage.misc.misc import cputime
# Bandersnatch curve
from lollipop489201 import Lollipop489201, Lollipop489201Point

p = 1910157204347957325700187962480217512925138482090399484362397
a=73275333332267847499581501376863252276520692179021512625126;
b=1538008641579097707704221968675032141849999412179326013460607;
r = 1910157204347957325700187962477453761504460135807772883063193
eigen = 500517855408530850485984845457693413895246088040882972578379

E = Lollipop489201(p, a, b, r, 1, eigen)
P = E.random_point()
n = random.randint(0,r)
print(n)

t = cputime()
P.scalar_mul(n)
t = cputime(t)
print("scalar multiplication:\t\t\t{:5.3f}ms".format(t))

P1 = E.random_point()
P2 = E.random_point()

k1 = random.randint(0,r)
k2 = random.randint(0,r)

t = cputime()
P.scalar_mul(n)
t = cputime(t)
print("scalar multiplication:\t\t\t{:5.3f}ms".format(t))

t2 = cputime()
P.fast_scalar_mul(n)
t2 = cputime(t2)
print('GLV:\t\t{:5.3f}ms ({:2.0f}% faster)'.format(t2, 100*(t-t2)/t2))