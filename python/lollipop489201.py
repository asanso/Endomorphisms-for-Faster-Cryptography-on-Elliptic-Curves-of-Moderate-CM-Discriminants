import sage.all
from sage.all import FractionField
from sage.rings.integer_ring import Z as ZZ
from sage.matrix.constructor import Matrix
from sage.modules.free_module_element import free_module_element as vector
from curve import WeierstrassCurve
from curve import PointWeierstrass

class  Lollipop489201(WeierstrassCurve):

       def __init__(self, p, a, b, r, cofactor,L):
            super().__init__(p, a,b, r, cofactor)
            self.D = -547
            self.L = L
            self.cofactor = cofactor
            self.r = r
            M = Matrix([[-L,1], [r,0]])
            self.N = M.LLL()
            self.N_inv = self.N**-1
            self.Fpxy = self.Fp['x', 'y']
            self.FpXY = FractionField(self.Fpxy)
            x, y = self.FpXY.gens()
            self.x0_map = (x**11 + 344959091081700922598753234080208529238293783251065431047762*x**10 + 618656784706780202005835730745508498961394889899594573602133*x**9 + 345046775957625697537452490799924294023258841238858303398288*x**8 + 329474990218736431418333824967260937394458529551615570946491*x**7 + 691728297863296108633527868364274276540466958086067387095638*x**6 + 1357990828765921024869300783866802391624650241598296038729547*x**5 + 1685216347350882591326020390968131766626485131724987747810528*x**4 + 1731680257507500767922993247383625512641276166785391715703265*x**3 + 206258873338213133803988372492280260862220393203645177382245*x**2 + 196909889892527505195484041589754696597614298956718028744136*x + 73521355122320838166117586475381237055438418773047977432418)/(x**10 + 344959091081700922598753234080208529238293783251065431047762*x**9 + 352206763105908364536525348823514096225834345932995725877219*x**8 + 1695504367068528763531374852931338650240409735184673879453190*x**7 + 1554766224953940578561535963601254196374826673521628393510673*x**6 + 713756646364889744922619364491054959936101349245024037005080*x**5 + 175213194816461071782468817950244523852037975793395744465459*x**4 + 458041619840079561526774390340896744813640059447166396713262*x**3 + 434818279101297896759676143702454744440592746076302229341350*x**2 + 495373531680873851138512740297031175502347015531087253549186*x + 516081654464090886604019685633066413036969723657457805496148)
            self.y0_map = (x**15*y + 517438636622551383898129851120312793857440674876598146571643*x**14*y + 1422805883596174815897152555190191305231430178199889242864234*x**13*y + 1235830203396332037342231702584403457060140993768197357746200*x**12*y + 1167158140851704127267172514520448125703155271915952985411840*x**11*y + 456048028442012048627378561237883653662609805799893458863180*x**10*y + 1534336943728983425292386694750981729344544471154362735643661*x**9*y + 711772070484181506589990348964463770754736827247877029827955*x**8*y + 1683757980840900324661482267629579663984837987483223179178481*x**7*y + 912226166872639571476828563993659630303059718942536037154585*x**6*y + 812804831647521380827713151385824707104134522078992566134775*x**5*y + 1269070611714168808805831756997274848818655773786457230641707*x**4*y + 631191777072740644729159478367012612095051457616463979190485*x**3*y + 726052255759378714070536380160393013950206642992737462675985*x**2*y + 1435966611309991001234369342698537076900972744407002246922282*x*y + 704184042616084674633071720690470858279640440021930108424115*y)/(x**15 + 517438636622551383898129851120312793857440674876598146571643*x**14 + 1689255905197046653366462937112185707966990722166488090589148*x**13 + 34184037023971003773513558157493389356190589711706974158149*x**12 + 194939472451628947553677133271819755265242139331575826191037*x**11 + 746860776195678464781702814837297514039620557481234570265595*x**10 + 245689221103848566439197527045061433678041811274024299542282*x**9 + 791698389953916123455561934848960593509116352014720763520985*x**8 + 1714039539474161075262684821036280740359688180554080022946109*x**7 + 974758460415018235105122873304820941738705495115246742896025*x**6 + 157593574622784075069994530838464719574945254591144557657770*x**5 + 1437124048537544763886933410010839295836389506437768612605644*x**4 + 1322867069493606242684309427454639040503730401885301513206879*x**3 + 658098423054303618211936394842433054629632148355581864298661*x**2 + 180767545670412292249655637241650859681554831261352253591473*x + 328727618918465900241876749242248711072761081420664765183639)
            self.x1_map = (x**11 + 749861857783311000014356152513809480380392700806514784754805*x**10 + 997291177211325205792785862936610252284658084024334301745370*x**9 + 810066573108720648530602943396138386481411011131483509441084*x**8 + 1054780613612252450032968808635709512834630733226824386828276*x**7 + 577692782843913339163583507087732602684890927438020408063199*x**6 + 945056594389552854914012177255801422868994665115724641248671*x**5 + 1227626843823511688209156821595022515275132821672484684922895*x**4 + 1550970379751372840827846019897581621656189794415350027614838*x**3 + 666996815352161636346872221816251913304351982520128255054223*x**2 + 873542576599292704141659019413582171582354377655559645064800*x + 660780789182856670738519540858284909653672248553371938267705)/(x**10 + 749861857783311000014356152513809480380392700806514784754805*x**9 + 1481472285823047983558513007074644713658356024207260839152791*x**8 + 394196683414760660123801058313157078893978216676595445546271*x**7 + 509110282304387659918500389726951220407244364452713684595528*x**6 + 433507335960871305095262156337016280812237868461663910427091*x**5 + 426994628681265832536850626666777832515735985846498734791310*x**4 + 1497278124599783278271294051624976103004089500892420548906896*x**3 + 1089144409304109311180596060827714780091476011751015231492402*x**2 + 1789715397434387899048223087068856453649121904903595233961854*x + 755170933323745962233074790393075359734761695185524492460588)
            self_y1_map = (x**15*y + 169714184500987837171440247530605464108019810164572434951009*x**14*y + 316526920063841588397385125819545739756400771206564869373286*x**13*y + 1621748525962821278244075682294049013731372771435883562918298*x**12*y + 1020604364678510187750094270772210460720794459682032557804251*x**11*y + 1769966919800560752973066299897586338452115743193252624033742*x**10*y + 522016255930273182727649989872394632484034118064938662728209*x**9*y + 1846052372446903973156404524974757989581031680226873681438897*x**8*y + 835903026492167665490860919650150008081472190360366227971309*x**7*y + 471407388196409080974211990263505404494943237189020722440138*x**6*y + 1489926204200302820134220093171232056006692261366655998565222*x**5*y + 363905601934455950976925103353037313325598776931687528440497*x**4*y + 99400210835168699769006621862787257514505818992994209599985*x**3*y + 699499579134806698630466360035284374524405593600176146418627*x**2*y + 651151254369888966351319266638044161252114955090809986456442*x*y + 1851537510780650617746242121464108055752524819273314818642073*y)/(x**15 + 169714184500987837171440247530605464108019810164572434951009*x**14 + 1742503015800076136331845944161728791307841313114037816328262*x**13 + 711680147034458979158052988175889350204011499326185009556440*x**12 + 1738885441951506731201190707875585146527572055387597591835977*x**11 + 223732728328134780947631657688890615520519500271365846565634*x**10 + 636631239210979590062161636346746623773714857267214622673445*x**9 + 1785278976854569441933648921208594033191681841468644419093652*x**8 + 1170101243343293969494756228373977720417185587773506218148504*x**7 + 407462923705238090437393384624573952338953445710929626383721*x**6 + 1685578154686000451271452047232040088421853028668941500286678*x**5 + 1519803727649360341027186592045515188041446554356386947963469*x**4 + 710733985312765847394099059022165290464746180390483229176598*x**3 + 691832373221706282676972249266298394100233935057397529673308*x**2 + 923742251077608245860858207477852941473179789842337294543876*x + 1413271646844725027533739430800669600223244366756053426600893)

class Lollipop489201Point(PointWeierstrass):
      
        def fast_scalar_mul(self, n):
            psiP = self.psi()
            beta = vector([n,0]) * self.curve.N_inv
            b = vector([int(beta[0]), int(beta[1])]) * self.curve.N
            k1 = n-b[0]
            k2 = -b[1]
            return self.multi_scalar_mul(k1, psiP, k2)