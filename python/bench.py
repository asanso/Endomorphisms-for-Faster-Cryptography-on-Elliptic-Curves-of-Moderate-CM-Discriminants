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
l = 101

t = cputime()
for i in range(1000):
    P.scalar_mul(2**l)
t = cputime(t)
print("scalar multiplication:\t\t\t{:5.3f}ms".format(t))

t2 = cputime()
for i in range(1000):
    P.psi()
t2 = cputime(t2)
print('endomorphism:\t\t{:5.3f}ms ({:2.0f}% faster)'.format(t2, 100*(t-t2)/t2))

print()

t = cputime()
for i in range(1000):
    P.scalar_mul(2**l)
t = cputime(t)
print("scalar multiplication:\t\t\t{:5.3f}ms".format(t))

t2 = cputime()
for i in range(1000):
    P.psi()
t2 = cputime(t2)
print('endomorphism:\t\t{:5.3f}ms ({:2.0f}% faster)'.format(t2, 100*(t-t2)/t2))