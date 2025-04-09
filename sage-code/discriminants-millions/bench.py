import random
import sys
import os
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parent_dir)

from sage.misc.misc import cputime
from  mnt6992 import MNT6992
from  mnt4992 import MNT4992
from lollipop956451 import Lollipop956451
from common.util import *

#mnt6-992
print_info("MNT6-992")

p = 25641744522180213985074042382512419032582205964092241104233192501069961626908990804689046062985807367763157388208148142783024211337097607475328728764340929123496640236440585428440814378281395763310908490919915070930671222454079808350406894268997985005951152391219105800051892376272420163532656250001
a = -3
b = 1642739491069868239864206518917328724337439725254750097543357143990096951788153282690860861421586090025840253466387680079357191502263632255953680449527905010857985902362206848357179784304415938666315990692980659025732260124353177900313122450061158787371088871707038627436604708188993662841042819334
r = 25641744522180213985074042382512419032582205964092241104233192501069961626908990804689046062985807367763157388208148142783024211337097607475328728764501059521991679940707134911388454490566722006414134338905237592345300668055898078034132669792309029396443022768100968479265087506725202295216329962501
eigen = 8755393583322235624319443723664913402874200942755340958163213971984906486122241114422935779373476910592110677400288398529023760704760112304369251845295565360428092853190053203262375244370427985275769122677973709446606120471880416662042575031699787390164434036694557784674658579937262190468706225568

E = MNT6992(p, a, b, r, 1, eigen)
P = E.random_point()
n = random.randint(0,r)
l = 496

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

#mnt4-992
print_info("MNT4-992")

p = 25641744522180213985074042382512419032582205964092241104233192501069961626908990804689046062985807367763157388208148142783024211337097607475328728764501059521991679940707134911388454490566722006414134338905237592345300668055898078034132669792309029396443022768100968479265087506725202295216329962501
a = -3
b = 4630025964795094055640335377828829116104647285415921057280481574803823968907844433725986809743446857347156975219184199865771148305416646193332038712269279350807907061256345533099160012744717289879975555577660031297353421101603686130897543020890317371517012884189187869991863601436037916413727384042
r = 25641744522180213985074042382512419032582205964092241104233192501069961626908990804689046062985807367763157388208148142783024211337097607475328728764340929123496640236440585428440814378281395763310908490919915070930671222454079808350406894268997985005951152391219105800051892376272420163532656250001
eigen = 17198569052751224804696743053088666217728203453423791031198203236527434056515615959555990921179642139177634032804218270656023986020928859889848991661550386815656843735007111202968777541321356404760681716004717117625495332026061638555511595773881997921089002756511880920416528347470121239406208218808


E = MNT4992(p, a, b, r, 1, eigen)
P = E.random_point()
n = random.randint(0,r)
l = 496

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

#lollipop-956451
print_info("lollipop-956451")

p = 346970567564890728424891183049337962035601978622541415965400614418353641616109413196732123410548601214270492880973832688318759029906877448237427055904522531717582919192550166662067828746656443559096610485904608891275458843426979796895475532274125193125711684703883234686434449369116690283
a = 327109485262376726136172643805932252878964212240064991046803990546570692113955204504295853421320565971168014448160518214189382593729022922934619191479697131508401422054977766988346990626462981083698419529403033388974770759484439350971607246476622466769281101195803851961334206387433757192
b = 279050255428736347835787201657572365478299307911044416076073175912778177975625581983202449094627882457913470267629624709093466687587530683922179493820089554667794629254883792790600467080376496365458891552310274342079659991659178694655626276630944082371547597372703391340583858884507859049
r = 4554474620279221025376588013428792331658085875800439965576470444142136125074436360422816311318800204172284882642774195680976231192452677
eigen = 1338228162335440229158233933458297877654311818170011752647462023091210257786969948374341757427185085939177917420605544163794803852157834
h = 76182347360104295691389198837779736649246768246453788622413701157110735036917351906437908803769129014598039561668077363397378012546375286033150471993610

E = Lollipop956451(p, a, b, r, h, eigen)
P = E.random_point()
n = random.randint(0,r)
l = 226

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