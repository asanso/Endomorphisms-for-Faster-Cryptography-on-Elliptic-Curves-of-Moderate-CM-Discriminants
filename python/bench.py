import random
from sage.misc.misc import cputime
from lollipop489201 import Lollipop489201
from lollipop574261 import Lollipop574261
from gc256c import GC256C
from util import *

# Lollipop489201

print_info("Lollipop489201")

p = 1910157204347957325700187962480217512925138482090399484362397
a=73275333332267847499581501376863252276520692179021512625126
b=1538008641579097707704221968675032141849999412179326013460607
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
    P.scalar_mul(n)
t = cputime(t)
print("scalar multiplication:\t\t\t{:5.3f}ms".format(t))

t2 = cputime()
for i in range(1000):
    P.fast_scalar_mul(n)
t2 = cputime(t2)
print('GLV:\t\t{:5.3f}ms ({:2.0f}% faster)'.format(t2, 100*(t-t2)/t2))

# Lollipop574261

print_info("Lollipop574261")

p = 2163160611951109656514578155686584121369935567340260530653195705987773102423967
a = 1600224790493789893768878688270410785959510786167424698577648903928123990501736
b = 442173884976592923728839681693322529291546262090176625621658004289727200834915
r = 2163160611951109656514578155686584121370430948758200237999126133244251338506461
eigen = 1158066151478108453899950316130401430242909130290171218593900211556999875716624

E = Lollipop574261(p, a, b, r, 1, eigen)
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
    P.scalar_mul(n)
t = cputime(t)
print("scalar multiplication:\t\t\t{:5.3f}ms".format(t))

t2 = cputime()
for i in range(1000):
    P.fast_scalar_mul(n)
t2 = cputime(t2)
print('GLV:\t\t{:5.3f}ms ({:2.0f}% faster)'.format(t2, 100*(t-t2)/t2))

#GC256C
print_info("GC256C")

p = 2**255 + 3225
a =  p -3
b = 28091019353058090096996979000309560759124368558014865957655842872397301267595
r = 57896044618658097711785492504343953927102133160255826820068844496087732066703
eigen = 25745384436187293773257686313148572361251961181400080445531962069440031481800

E = GC256C(p, a, b, r, 1, eigen)
E = Lollipop489201(p, a, b, r, 1, eigen)
P = E.random_point()
n = random.randint(0,r)
l = 128


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
    P.scalar_mul(n)
t = cputime(t)
print("scalar multiplication:\t\t\t{:5.3f}ms".format(t))

t2 = cputime()
for i in range(1000):
    P.fast_scalar_mul(n)
t2 = cputime(t2)
print('GLV:\t\t{:5.3f}ms ({:2.0f}% faster)'.format(t2, 100*(t-t2)/t2))
