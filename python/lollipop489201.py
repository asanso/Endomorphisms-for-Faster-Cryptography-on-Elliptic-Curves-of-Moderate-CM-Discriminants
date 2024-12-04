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
        self.a = a
        self.b = b
        self.D = -547
        self.L = L
        self.cofactor = cofactor
        self.r = r
        M = Matrix([[-L,1], [r,0]])
        self.N = M.LLL()
        self.N_inv = self.N**-1
        self.Fpxz = self.Fp['x', 'z']
        self.FpXZ = FractionField(self.Fpxz)
        x, z = self.FpXZ.gens()
        
        self.a0_map = x**16 + 517438636622551383898129851120312793857440674876598146571643*x**15*z + 45548722449961165135585356553962597777412784042687453951665*x**14*z**2 + 1004771838759537490760840541150595414258965276111549598306288*x**13*z**3 + 1802196064355541574277856509986400465783152521716485431707157*x**12*z**4 + 171394251038162832245290564422946125007474434422483245271934*x**11*z**5 + 1766869533051641042659913569656177082454393263090678097158146*x**10*z**6 + 669791093723921388876352212349937218928107312052467956586541*x**9*z**7 + 1611812048340327603707858737917294143177814170729864950224756*x**8*z**8 + 292140100130670505835712820807802200028929972166633187148159*x**7*z**9 + 1723878874686292508980675111944388275406384272060486097798484*x**6*z**10 + 1133077027993103686431749500288479289997885486174502352337203*x**5*z**11 + 596082167637905227071725876381115039405406907212206354633101*x**4*z**12 + 476012972870516618243827403164107892260507939360409723419564*x**3*z**13 + 1286863496611669913418526005933126942968694800583081878036525*x**2*z**14 + 707734358177489618446932539101928861059544642930875191744900*x*z**15 + 472188838773936911376041873658187425992126159503843950370949*z**16
        self.b0_map = x**15 + 517438636622551383898129851120312793857440674876598146571643*x**14*z + 1422805883596174815897152555190191305231430178199889242864234*x**13*z**2 + 1235830203396332037342231702584403457060140993768197357746200*x**12*z**3 + 1167158140851704127267172514520448125703155271915952985411840*x**11*z**4 + 456048028442012048627378561237883653662609805799893458863180*x**10*z**5 + 1534336943728983425292386694750981729344544471154362735643661*x**9*z**6 + 711772070484181506589990348964463770754736827247877029827955*x**8*z**7 + 1683757980840900324661482267629579663984837987483223179178481*x**7*z**8 + 912226166872639571476828563993659630303059718942536037154585*x**6*z**9 + 812804831647521380827713151385824707104134522078992566134775*x**5*z**10 + 1269070611714168808805831756997274848818655773786457230641707*x**4*z**11 + 631191777072740644729159478367012612095051457616463979190485*x**3*z**12 + 726052255759378714070536380160393013950206642992737462675985*x**2*z**13 + 1435966611309991001234369342698537076900972744407002246922282*x*z**14 + 704184042616084674633071720690470858279640440021930108424115*z**15
        self.c0_map = x**15 + 517438636622551383898129851120312793857440674876598146571643*x**14*z + 1689255905197046653366462937112185707966990722166488090589148*x**13*z**2 + 34184037023971003773513558157493389356190589711706974158149*x**12*z**3 + 194939472451628947553677133271819755265242139331575826191037*x**11*z**4 + 746860776195678464781702814837297514039620557481234570265595*x**10*z**5 + 245689221103848566439197527045061433678041811274024299542282*x**9*z**6 + 791698389953916123455561934848960593509116352014720763520985*x**8*z**7 + 1714039539474161075262684821036280740359688180554080022946109*x**7*z**8 + 974758460415018235105122873304820941738705495115246742896025*x**6*z**9 + 157593574622784075069994530838464719574945254591144557657770*x**5*z**10 + 1437124048537544763886933410010839295836389506437768612605644*x**4*z**11 + 1322867069493606242684309427454639040503730401885301513206879*x**3*z**12 + 658098423054303618211936394842433054629632148355581864298661*x**2*z**13 + 180767545670412292249655637241650859681554831261352253591473*x*z**14 + 328727618918465900241876749242248711072761081420664765183639*z**15

        self.a1_map = x**16 + 169714184500987837171440247530605464108019810164572434951009*x**15*z + 1258321907188353358566118800023694329934143372931111278920841*x**14*z**2 + 959200990696743917764293375062775423568532672710148270240340*x**13*z**3 + 1459288548386100856751346242083650278921011403280076948667214*x**12*z**4 + 1221444656052679258569022724875817397051118512787965159227345*x**11*z**5 + 233868830000414568252156415597688031219570001693680530151948*x**10*z**6 + 292110529114196314192941152604518402384926068350862969939065*x**9*z**7 + 1896788176091141740812394784081948467510354311132087627746074*x**8*z**8 + 154006570746017283064665376978825879215925437330009283229147*x**7*z**9 + 342001349241338896389959480198412214946877879124345235088729*x**6*z**10 + 1756885581559708330285451739999582441225484492461434608234645*x**5*z**11 + 1389482338735724714866503021906927571650435595450174128908462*x**4*z**12 + 910727020050103194785937317895820973436669392844426071887777*x**3*z**13 + 1758333466106987618743909711502302623369418465232133023344959*x**2*z**14 + 234218957324430311613232941363280549391988680090564532274889*x*z**15 + 1156165718476500662284098733393459561165691222571434936767925*z**16
        self.b1_map = x**15 + 169714184500987837171440247530605464108019810164572434951009*x**14*z + 316526920063841588397385125819545739756400771206564869373286*x**13*z**2 + 1621748525962821278244075682294049013731372771435883562918298*x**12*z**3 + 1020604364678510187750094270772210460720794459682032557804251*x**11*z**4 + 1769966919800560752973066299897586338452115743193252624033742*x**10*z**5 + 522016255930273182727649989872394632484034118064938662728209*x**9*z**6 + 1846052372446903973156404524974757989581031680226873681438897*x**8*z**7 + 835903026492167665490860919650150008081472190360366227971309*x**7*z**8 + 471407388196409080974211990263505404494943237189020722440138*x**6*z**9 + 1489926204200302820134220093171232056006692261366655998565222*x**5*z**10 + 363905601934455950976925103353037313325598776931687528440497*x**4*z**11 + 99400210835168699769006621862787257514505818992994209599985*x**3*z**12 + 699499579134806698630466360035284374524405593600176146418627*x**2*z**13 + 651151254369888966351319266638044161252114955090809986456442*x*z**14 + 1851537510780650617746242121464108055752524819273314818642073*z**15
        self.c1_map = x**15 + 169714184500987837171440247530605464108019810164572434951009*x**14*z + 1742503015800076136331845944161728791307841313114037816328262*x**13*z**2 + 711680147034458979158052988175889350204011499326185009556440*x**12*z**3 + 1738885441951506731201190707875585146527572055387597591835977*x**11*z**4 + 223732728328134780947631657688890615520519500271365846565634*x**10*z**5 + 636631239210979590062161636346746623773714857267214622673445*x**9*z**6 + 1785278976854569441933648921208594033191681841468644419093652*x**8*z**7 + 1170101243343293969494756228373977720417185587773506218148504*x**7*z**8 + 407462923705238090437393384624573952338953445710929626383721*x**6*z**9 + 1685578154686000451271452047232040088421853028668941500286678*x**5*z**10 + 1519803727649360341027186592045515188041446554356386947963469*x**4*z**11 + 710733985312765847394099059022165290464746180390483229176598*x**3*z**12 + 691832373221706282676972249266298394100233935057397529673308*x**2*z**13 + 923742251077608245860858207477852941473179789842337294543876*x*z**14 + 1413271646844725027533739430800669600223244366756053426600893*z**15

        self.a2_map = x**16 + 249412976004351119321546162401705573618734615646130015438663*x**15*z + 331586010354288757593815988330979866551347781590323759813063*x**14*z**2 + 189493671609550927101538751995202510433674731790506477035820*x**13*z**3 + 423493981085860762879278009647485952463906869324804812246997*x**12*z**4 + 910943213355410029737284549502320841872943381765507823532423*x**11*z**5 + 876090887633442224452397984143883009879568347183264345836008*x**10*z**6 + 1006826018526071345920730667795961090516611185448796217947641*x**9*z**7 + 1815860317046454549085567022955354622945022423307526899611584*x**8*z**8 + 1430441365535157702499606748092150821225009053082384704236158*x**7*z**9 + 949309344355634515639425813172517804106187246114581178059862*x**6*z**10 + 140840236073205946956469089411799418798821525300992431407092*x**5*z**11 + 1566590983266380684984011251595378489786022763547768678275206*x**4*z**12 + 1766183926763734121106923054068582473933814437311998201417097*x**3*z**13 + 869773147092731717284723641670006979617756581517945560780021*x**2*z**14 + 43814905145794841824853879671867046097939560633324832249850*x*z**15 + 1363240216739756793966996039902707518383662609834952341410940*z**16
        self.b2_map = x**15 + 249412976004351119321546162401705573618734615646130015438663*x**14*z + 1384829807186116028194758208377613837577045136949046726192020*x**13*z**2 + 1669970933976859231005510236665805921400495534188185851252026*x**12*z**3 + 794890973360162428414987632315972344016820540669556741229124*x**11*z**4 + 1164417450400411974623164806128702779941378045096401437000697*x**10*z**5 + 551831035483681211011764789535463288424217754710852509327005*x**9*z**6 + 1885737894326071143269784986072018146362326287530457095928904*x**8*z**7 + 1590922590099200992405113920664146203981328707557282384201178*x**7*z**8 + 867742116421271536787888916969681222580122510846044064891353*x**6*z**9 + 744409225891123561974923259717318337604606620337222120384812*x**5*z**10 + 1206895556990867769185235914140807234199121602002569596101045*x**4*z**11 + 1731983330995967981289354227865318361194636445533149450203467*x**3*z**12 + 1737056527971484937617597109637790025504150093869326964767799*x**2*z**13 + 1018044081109230946373749874799231292102607804766975906751099*x*z**14 + 474288560065795157389500767946407971425462752751174093722746*z**15
        self.c2_map = x**15 + 249412976004351119321546162401705573618734615646130015438663*x**14*z + 1813286510944181055744381079594405608526765700314884985183740*x**13*z**2 + 1499881329535285102607130415844196244467728270804278375715389*x**12*z**3 + 1519576638283851453594183806213673143457505961576993747138285*x**11*z**4 + 1160705524494879067407165518213407500636075761687346515359473*x**10*z**5 + 168289226021318157423925338284968952448556615905064702492552*x**9*z**6 + 1567172386010113319667824428725146094924025785040220539790973*x**8*z**7 + 1310958242764537731145210729145101620952796488477196222458734*x**7*z**8 + 139711277389620670125373681944854353861220199740102356131869*x**6*z**9 + 1150857506775430358811264708976325613298796444800826749897342*x**5*z**10 + 1661095751470907232121236741776204593341884599851029629939829*x**4*z**11 + 1376732737741033707625712959871528652916497115945013045837396*x**3*z**12 + 1584971691504606307458577682272465104509087252715544159060104*x**2*z**13 + 1363573205342400802718856090135950303393594404518014812358818*x*z**14 + 1784423724414520149961100580555710058712897684772207396181465*z**15
        
        self.iso_x = 210189393351005824025698136604114443811292381716765236542912
        self.iso_y = 686904640397606125397390333239731317785401578633329391615079
        
    def random_point(self):
        P = super().random_point()
        return Lollipop489201Point(P.X, P.Y,P.Z, self)
       
class Lollipop489201Point(PointWeierstrass):
      
        def __init__(self, x, y, z, curve):
            super().__init__(x, y, z, curve)

        def fast_scalar_mul(self, n):
            psiP = self.psi()
            beta = vector([n,0]) * self.curve.N_inv
            b = vector([int(beta[0]), int(beta[1])]) * self.curve.N
            k1 = n-b[0]
            k2 = -b[1]
            return self.multi_scalar_mul(k1, psiP, k2)

        def psi(self):
            x0,y0,z0 = self.X, self.Y, self.Z
            #1st isogeny
            x1 = self.curve.a0_map(x0,z0) 
            y1 =  y0* self.curve.b0_map(x0,z0) 
            z1 =  z0* self.curve.c0_map(x0,z0) 
            #2nd isogeny
            x2 = self.curve.a1_map(x1,z1) 
            y2 =  y1* self.curve.b1_map(x1,z1) 
            z2 =  z1* self.curve.c1_map(x1,z1) 
            #3rd isogeny

            x3 = self.curve.a2_map(x2,z2) 
            y3 =  y2* self.curve.b2_map(x2,z2) 
            z3 =  z2* self.curve.c2_map(x2,z2) 
            return Lollipop489201Point( self.curve.iso_x*x3,  self.curve.iso_y*y3, z3, self.curve)
