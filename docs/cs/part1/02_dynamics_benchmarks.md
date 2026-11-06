# Dynamika: dráhy, sekulární vývoj a stabilita {#dynamics}

Tato kapitola ověřuje, jak SWINGERsim pohybuje hvězdami: N-částicový backend s relativitou,
sekulární backendy, detekci nestability a reprodukci čtyř klasických studií hierarchických soustav.
Ověřuje také hierarchie s více než čtyřmi hvězdami, rozptylové experimenty a vliv okolí (plyn,
slapy Galaxie, procházející hvězdy, gravitační vlny).

V celé zprávě se opakují dva jevy. **Apsidální stáčení** je pomalé otáčení eliptické dráhy v její
vlastní rovině (argument periastra $\omega$ roste); jeho příčinami jsou ostatní hvězdy, obecná
relativita a slapové a rotační výdutě hvězd. **Precese uzlů** je otáčení roviny dráhy kolem
celkového momentu hybnosti soustavy; mění sklon k rovině oblohy, a tím způsobuje objevování
a mizení zákrytů. Je-li vzájemný sklon vnitřní a vnější dráhy velký (zhruba 40°–140°), vnější
těleso si vyměňuje moment hybnosti s vnitřní dráhou a periodicky zvyšuje a snižuje její
excentricitu: to jsou **cykly von Zeipelovy–Lidovovy–Kozaiovy (ZLK)** (Kozai 1962; Lidov 1962;
Naoz 2016). Charakteristickou dobou těchto cyklů je **Kozaiova doba** $t_K$, řádově
$P_\mathrm{out}^2/P_\mathrm{in}$ krát poměr hmotností.

## Relativita, cykly ZLK a precese uzlů v N-částicovém backendu {#nbody-basics}

**Reference a vstup.** N-částicový backend ověřují tři analytické limity. (i) Stáčení perihelia
Merkuru v prvním postnewtonovském řádu (1PN), 42,98″ za století pro skutečnou dráhu
($a$ = 0,387 au, $e$ = 0,206), jeden z klasických testů obecné relativity (Will 2014). (ii) Cykly
ZLK testovací částice ve vzdálenosti 1 au od hvězdy o hmotnosti 1 M☉, rušené hvězdou o hmotnosti
1 M☉ na kruhové dráze ve vzdálenosti 20 au. V této kvadrupólové limitě testovací částice dává teorie
$e_\mathrm{max} = \sqrt{1 - \tfrac53\cos^2 i_0}$ pro počáteční vzájemný sklon $i_0$, vzájemný
sklon 39,23° v maximu excentricity, žádný růst excentricity pod 39,2° a přesně zachovávanou
veličinu $\sqrt{1-e^2}\cos i$ (Kozai 1962; Lidov 1962; Naoz 2016). (iii) Perioda precese uzlů
téhož uspořádání při $i$ = 20° podle kvadrupólového vzorce Söderhjelma (1975). Sekulární backendy
se ověřují ve stejných limitách a navíc na relativistické apsidální periodě osamocené dvojhvězdy
($a$ = 0,05 au, $e$ = 0,3).

**Tabulka 3.** Relativita, cykly ZLK a precese: N-částicový backend (třída R) a sekulární
backendy (třída own). Sloupce jsou definovány v oddílu [Jak číst tuto zprávu](#how-to-read).

| Test | Třída | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|---|
| `gr-mercury` | R | 42,98″/století | 42,9803″/století | +0,00 % | 0,5 % | PASS |
| `zlk-emax-65` | R | 0,838047 | 0,84023 | +0,26 % | 1 % | PASS |
| `zlk-emax-80` | R | 0,974548 | 0,976155 | +0,16 % | 1 % | PASS |
| `zlk-imin-65` | R | 39,2315° | 39,3441° | +0,113° | 0,5° | PASS |
| `zlk-below-critical` | R | 0 | 0,00173037 | +0,00173 | 0,01 | PASS |
| `zlk-kozai-constant` | R | 0 | 5,16811·10⁻⁵ | +5,17·10⁻⁵ | 0,001 | PASS |
| `nodal-precession` | R | 11351,4 roku | 11339,8 roku | −0,10 % | 1 % | PASS |
| `zlk-period-nbody` | own | 0 | 0,00660493 | +0,0066 | 0,02 | PASS |
| `quad-tp-emax` | own | 0,838047 | 0,838082 | +0,00 % | 0,1 % | PASS |
| `secular-kozai-emax` | own | 0,838047 | 0,838082 | +0,00 % | 0,5 % | PASS |
| `secular-3p1-kozai-emax` | own | 0,838047 | 0,838082 | +0,00 % | 0,5 % | PASS |
| `secular-nodal-precession` | own | 11351,4 roku | 11351,4 roku | −0,00 % | 0,5 % | PASS |
| `secular-gr-apsidal` | own | 7113,77 roku | 7113,77 roku | −0,00 % | 0,5 % | PASS |
| `brouwer-law` | R | 0,5 | 0,59301 | +0,093 | 0,15 | PASS |

**Výsledek.** Tabulka 3 uvádí srovnání. N-částicový backend dává stáčení perihelia Merkuru
s přesností lepší než 10⁻⁴; s vypnutou relativitou dává 0″, takže test odhalí chybějící fyziku.
Maxima ZLK z N-částicových běhů leží o 0,16–0,26 % nad vzorcem. Tento malý přebytek je očekávaný
a má správné znaménko: vzorec zanedbává členy řádu $a_\mathrm{in}/a_\mathrm{out}$ (zde 1/20)
a N-částicová dráha začíná s $e$ = 0,01, nikoli přesně s nulou. Veličina $\sqrt{1-e^2}\cos i$,
středovaná přes první a poslední čtvrtinu běhu dlouhého 20 000 let, se posune o 5·10⁻⁵, zatímco
$e$ kmitá mezi 0 a 0,84; její okamžitá hodnota kolísá asi o 0,5 %, protože elementy jsou
oskulační (spočtené z okamžitých poloh a rychlostí). Pod kritickým sklonem zůstává excentricita
pod 0,002. Perioda ZLK z N-částicového běhu je o 0,66 % delší než kvadrupólová sekulární perioda.
Perioda precese uzlů souhlasí se Söderhjelmovým vzorcem na 0,10 %.

Tři sekulární backendy souhlasí se stejnými analytickými limitami s přesností lepší než 10⁻⁴
(rozdíly 0,00 % ve sloupci „Rozdíl“). To je očekávané a je to silný test jejich implementace:
řeší přesně středované rovnice těchto limit, bez krátkoperiodických členů, kvůli nimž se
N-částicové hodnoty mírně odchylují. Hodnotu $e_\mathrm{max}$ testovací částice reprodukuje kód pro
trojhvězdy, větev 3+1 vnořeného kódu i kód 2+2, v němž je rušící těleso nahrazeno dvojhvězdou
zanedbatelné velikosti.

**Numerická kvalita integrátoru.** Chyba energie integrátoru vysokého řádu se zaokrouhlovacími
chybami roste jako náhodná procházka, $\propto t^{1/2}$ (Brouwerův zákon; Rein & Spiegel 2015).
Šest realizací přes 3·10⁴ oběhů dává exponent růstu 0,59 (`brouwer-law`), v rámci tolerance 0,15
dané rozptylem šesti realizací. Jde o test integrátoru IAS15 tak, jak jej SWINGERsim používá,
nikoli o test fyziky SWINGERsimu.

## Stabilita: změny hierarchie a indikátory chaosu {#stability-basics}

Hierarchická konfigurace je **dynamicky stabilní**, jestliže si velmi dlouho uchová svou hierarchii
(kdo obíhá kolem koho). **Kritérium Mardlinga–Aarsetha** (MA01; Mardling & Aarseth 2001) udává
nejmenší poměr vzdálenosti vnějšího periastra k velké poloose vnitřní dráhy, při němž je
trojhvězda stabilní. **MEGNO** (Mean Exponential Growth factor of Nearby Orbits; Cincotta & Simó
2000) měří, jak rychle se rozbíhají sousední trajektorie: pro pravidelný (kvaziperiodický) pohyb
směřuje k 2 a pro chaotický pohyb roste neomezeně. SWINGERsim počítá MEGNO pomocí variačních rovnic
knihovny REBOUND.

**Tabulka 4.** Detekce nestability a chaosu (třída R).

| Test | Třída | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|---|
| `stability-unstable` | R | 1 | 1 | +0 | 0 | PASS |
| `stability-stable` | R | 0 | 0 | +0 | 0 | PASS |
| `megno-regular` | R | 2 | 1,99959 | −0,00041 | 0,05 | PASS |
| `megno-chaotic` | R | 1 | 1 | +0 | 0 | PASS |

Koplanární čtyřhvězda 2+2 se stejnými hmotnostmi se umístí na polovinu a na 2,5násobek meze MA01
a integruje se po dobu 100 vnějších period. První musí změnit hierarchii (1 = ano) a také ji
změní, když asi po dvou vnějších obězích vymrští hvězdu; druhá ji změnit nesmí (0 = beze změny)
a nezmění. Rezerva dvojnásobku na obě strany je nutná, protože MA01 bylo odvozeno pro trojhvězdy,
nikoli pro čtyřhvězdy. MEGNO pravidelné hierarchické trojhvězdy konverguje k 1,998 (očekáváno 2)
a čtyřhvězda 2+2 uvnitř nestabilní oblasti dosáhne MEGNO > 4 (1 = ano). Nástroje pro stabilitu
jako celek se ověřují proti publikovanému neuronovému klasifikátoru v [kapitole 9](#stability).

## Reprodukce Naoz et al. (2013) {#naoz2013}

Naoz et al. (2013) ukázali, že členy oktupólového řádu – další řád za kvadrupólem v rozvoji podle
$a_\mathrm{in}/a_\mathrm{out}$, důležitý při excentrické vnější dráze – mohou překlopit vnitřní
dráhu z prográdní na retrográdní. Dva jejich obrázky slouží jako srovnávací testy úrovně B se
vstupními hodnotami převzatými z popisků a z textu.

**Vstup.** *Obr. 3* (excentrické překlopení ZLK): hvězda o hmotnosti 1 M☉ s planetou o hmotnosti
Jupiteru ve vzdálenosti 6 au a hnědým trpaslíkem o 40 hmotnostech Jupiteru ve vzdálenosti 100 au
s $e$ = 0,6, počáteční vzájemný sklon 65°. *Obr. 4* (uspořádání PSR B1620−26): hmotná vnitřní
dvojhvězda, která se překlápí v kvadrupólovém řádu. Text udává rozsah vzájemného sklonu
57,5°–106,7°. Popisek uvádí počáteční sklony obou drah k celkovému momentu hybnosti 6,75°
a 63,25°; tomu odpovídá poměr momentů hybnosti
$G_1/G_2 = \sin 63{,}25°/\sin 6{,}75° = 7{,}60$, kdežto uvedené hmotnosti a $a_2$ = 50 au dávají
6,79. Nezávislý kvadrupólový výpočet reprodukuje publikovaný rozsah s poměrem 7,60
(57,53°–106,70°), nikoli s 6,79 (57,3°–105,0°). Test proto používá $a_2$ = 39,9 au, což dává
publikovaný poměr.

**Tabulka 5.** Naoz et al. (2013), obr. 3 a 4: N-částicový backend (třída R) a sekulární backendy
(own).

| Test | Třída | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|---|
| `naoz13-fig4-imin` | R | 57,5° | 57,2391° | −0,261° | 0,5° | PASS |
| `naoz13-fig4-imax` | R | 106,7° | 106,93° | +0,23° | 0,5° | PASS |
| `quad-naoz13-fig4-imin` | own | 57,5° | 57,5376° | +0,0376° | 0,1° | PASS |
| `quad-naoz13-fig4-imax` | own | 106,7° | 106,721° | +0,0207° | 0,1° | PASS |
| `naoz13-fig3-flip-time` | R | 3,8 Myr | 3,9868 Myr | +4,92 % | 10 % | PASS |
| `naoz13-fig3-imax` | R | 142,6° | 143,01° | +0,408° | 2° | PASS |
| `secular-naoz13-fig3-flip-time` | own | 3,8 Myr | 4,04981 Myr | +6,57 % | 10 % | PASS |
| `secular-naoz13-fig3-imax` | own | 142,6° | 142,556° | −0,0435° | 2° | PASS |

![**Obrázek 2.** SWINGERsim reprodukuje Naoz et al. (2013). Oba panely ukazují vzájemný sklon vnitřní a vnější dráhy (svislá osa, stupně) v závislosti na čase; vodorovná černá čára vyznačuje 90°, hranici mezi prográdním a retrográdním pohybem. *Vlevo:* jejich obr. 3, excentrické překlopení ZLK, N-částicový běh přes 5 Myr (modrá čára). Čárkovaná červená svislá čára je první překlopení odečtené z publikovaného obrázku (3,8 Myr), tečkovaná červená vodorovná čára publikované maximum po překlopení (142,6°). *Vpravo:* jejich obr. 4, uspořádání PSR B1620−26 přes 500 kyr; modře s poměrem momentů hybnosti z článku ($G_1/G_2$ = 7,60), šedě s uvedenou hodnotou $a_2$ = 50 au ($G_1/G_2$ = 6,79); tečkované červené čáry vyznačují publikovaný rozsah 57,5°–106,7°.](../../figures/part1/naoz2013.png)

**Výsledek.** Obrázek 2 ukazuje obě uspořádání. V levém panelu si všimněte minim ZLK před
překlopením: leží na 38° a jejich obálka pomalu stoupá k 90°, stejně jako v publikovaném obrázku.
Dráha poprvé překročí 90° v čase 3,99 Myr v N-částicovém běhu a 4,05 Myr v sekulárním backendu,
proti 3,8 Myr odečteným z publikovaného obrázku (+4,9 % a +6,6 %; tolerance 10 % odráží odečítání
času z obrázku a citlivost okamžiku překlopení na fázi cyklů). Po překlopení dosáhne sklon 143,0°
(N-částicově) a 142,6° (sekulárně), proti 142,6°. V pravém panelu kmitá modrá křivka přesně mezi
oběma tečkovanými čarami: sekulární kvadrupólový kód reprodukuje publikovaný rozsah na 0,04°
(57,54°–106,72°), N-částicový backend na 0,26° (57,24°–106,93°); rozdíl tvoří krátkoperiodické
členy, které středované rovnice vynechávají. Šedá křivka ukazuje, že uvedená hodnota $a_2$ = 50 au
by dala užší rozsah; publikovaná křivka patří k poměru z popisku.

## Reprodukce Vokrouhlického (2016) {#vokrouhlicky2016}

Vokrouhlický (2016) vypracoval sekulární teorii **čtyřhvězdy 2+2** – dvou dvojhvězd A a B na
společné vnější dráze – a použil ji na soustavu podobnou V994 Her. Veličinou, která jej zajímá, je
$J_{AB}$, úhel mezi rovinami drah A a B: určuje, zda mohou obě dvojhvězdy zakrývat současně.
Odvodil analytickou horní mez $J_{AB}$ a znázornil čtyři uspořádání (své obr. 2, 3, 4 a 7).

**Vstup.** Soustava podobná V994 Her z článku, sestavená z elementů vztažených k rovině oblohy:
obr. 2 s kruhovými drahami A a B a uzly $\Omega_A$ = 10°, $\Omega_B$ = 15°; obr. 3 s pozorovanými
$e_A$ = 0,031, $e_B$ = 0,125; obr. 4 s $\Omega_A$ = 35°, $\Omega_B$ = 40°; obr. 7 jako obr. 4 plus
slapová kvadrupólová deformace hvězd. Pro obr. 7 článek ladí jednu deformovatelnost pro každou
dvojhvězdu na pozorované apsidální periody asi 116 let (A) a 111 let (B), s poloměry
$R\propto M^{0{,}8}$ (2,4, 2,1, 1,6 a 1,7 R☉). SWINGERsim používá stejné poloměry, slapy bez
rotace a disipace a apsidální konstanty $k_2$ = 0,00754 (A) a 0,00333 (B), pro něž Sterneův vzorec
dává přesně tyto apsidální periody osamocených dvojhvězd. Čtyři uspořádání byla zintegrována
N-částicovým backendem (každé 25–60 minut) a sekulárním backendem 2+2 (každé sekundy).

**Tabulka 6.** Vokrouhlický (2016): úhly a pravděpodobnosti zákrytů (vše třída own). Sloupce jsou definovány v oddílu [Jak číst tuto zprávu](#how-to-read).

| Test | Třída | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|---|
| `vok16-v994-JA` | own | 11,6° | 11,6462° | +0,0462° | 0,1° | PASS |
| `vok16-v994-JB` | own | 15,8° | 15,7932° | −0,00678° | 0,1° | PASS |
| `vok16-jab-max-analytic` | own | 5,2° | 5,20074° | +0,000736° | 0,05° | PASS |
| `vok16-fig2-jab-min` | own | 3,5° | 3,47152° | −0,0285° | 0,15° | PASS |
| `vok16-fig2-jab-max` | own | 5,2° | 5,20069° | +0,000687° | 0,15° | PASS |
| `vok16-fig2-jab-1myr` | own | 5,2° | 5,2006° | +0,000633° | 0,15° | PASS |
| `vok16-fig3-jab-max` | own | 6° | 6,20874° | +0,209° | 0,5° | PASS |
| `vok16-eclipse-prob-A` | own | 0,366 | 0,369265 | +0,89 % | 2 % | PASS |
| `vok16-eclipse-prob-B` | own | 0,397 | 0,399075 | +0,52 % | 2 % | PASS |

První dva řádky tabulky 6 ověřují převod elementů vztažených k rovině oblohy na sklony $J_A$, $J_B$
dvojhvězd k vnější dráze ve SWINGERsimu – nezávislá kontrola konvencí úhlů. Analytická mez
$J_{AB}$ je 5,201° (článek 5,2°). Sekulární backend reprodukuje rozsahy $J_{AB}$ uvedené v textu
pro obr. 2 a 3 na 0,03° a 0,21° a během 1 Myr zůstává maximum na analytické mezi (článek zjišťuje
totéž během 100 Myr). Kritérium zákrytu SWINGERsimu použité na izotropně orientované dráhy dává
pravděpodobnosti, že A a B zakrývají, 0,369 a 0,399, proti publikovaným 0,366 a 0,397.

**Tabulka 7.** Vokrouhlický (2016), obr. 2, 3, 4 a 7: publikované chování proti N-částicovému
backendu (třída R) a sekulárnímu backendu 2+2 (třída own). Každý řádek je jedno uspořádání; druhý
sloupec uvádí, co ukazuje článek, třetí a čtvrtý rozsahy, které spočte SWINGERsim.

| Uspořádání | Článek | SWINGERsim N-částicově | SWINGERsim sekulárně |
|---|---|---|---|
| Obr. 2: kruhové A a B | $J_{AB}$ 3,5°–5,2°, $e_A$, $e_B$ blízko nuly | $J_{AB}$ 3,41°–5,19° během 50 kyr, $e$ < 0,005 | $J_{AB}$ 3,47°–5,20°, $e_A$ < 0,0051 |
| Obr. 3: pozorované $e_A$, $e_B$ | $J_{AB}$ ≈ 1°–6°, malé oscilace $e$ | $J_{AB}$ 0,94°–6,30°, $e_A$ 0,029–0,036, $e_B$ 0,123–0,140 | $J_{AB}$ 0,78°–6,21°, $e_A$ 0,029–0,036, $e_B$ 0,123–0,138 |
| Obr. 4: $\Omega_A$ = 35°, $\Omega_B$ = 40° | $J_A$ ≈ 35°, $J_B$ ≈ 40°, $e_B$ nad 0,6, $J_{AB}$ od 0° do více než 70° | $J_A$ = 35,4°, $J_B$ = 40,3°, $e_B$ až 0,631, $J_{AB}$ 3,5°–73,4° | $e_B$ až 0,610, $J_{AB}$ 3,2°–71,1° |
| Obr. 7: obr. 4 + slapová deformace | obnovený řád: malé oscilace $e$ a $i$ | $e_B$ 0,124–0,136, $e_A$ 0,029–0,032, $J_{AB}$ 0,4°–5,8°, apsidální periody 113 a 108 let | $e_B$ 0,123–0,134, $J_{AB}$ 0,29°–5,91°, apsidální periody 113,1 a 108,4 roku |

![**Obrázek 3.** SWINGERsim reprodukuje Vokrouhlického (2016), obr. 2, 3, 4 a 7 (jeden sloupec pro každé uspořádání). *Horní řada:* úhel $J_{AB}$ mezi rovinami drah obou dvojhvězd (stupně) v závislosti na čase (kyr); černě: N-částicový backend; červeně čárkovaně: sekulární backend 2+2; šedě čárkovaně: analytická mez 5,2° z článku. *Dolní řada:* excentricity A (modře) a B (oranžově; červeně čárkovaně pro sekulární běh) v závislosti na čase.](../../figures/vokrouhlicky2016.png)

**Výsledek.** Obrázek 3 a tabulka 7 srovnávají čtyři uspořádání. V prvním sloupci kmitá $J_{AB}$
mezi 3,4° a 5,2° a nikdy nepřekročí šedou analytickou mez; černá a červená křivka se překrývají.
Ve druhém sloupci dovolují pozorované excentricity, aby $J_{AB}$ kmital mezi asi 1° a 6°, jako
v článku. Zajímavý je třetí sloupec: při větších sklonech vstoupí B do cyklů ZLK a dosáhne
$e_B$ = 0,63, $J_{AB}$ putuje až k 73° a N-částicová a sekulární křivka souhlasí asi 4 kyr a pak
se rozejdou. Článek popisuje přesně totéž pro své dva integrátory („they start to significantly
deviate after ≃5 kyr“): pohyb je chaotický a dvě správné integrace se rozbíhají. V tomto
uspořádání dosáhne nejbližší periastrum B 1,06násobku součtu poloměrů hvězd, což podporuje závěr
článku, že model hmotných bodů je zde neúplný. Čtvrtý sloupec přidává slapové výdutě: rychlé
apsidální stáčení, které způsobují (113 a 108 let), rozladí cykly ZLK dvojhvězdy B, její
excentricita zůstane na 0,124–0,136 místo aby dosáhla 0,63, a roviny A a B zůstanou vůči sobě
v rozmezí 6° – řád, který článek obnovuje ve svém obr. 7. Energie se v N-částicových bězích
zachovává na 2,5·10⁻¹². Kvantitativně jsou publikované rozsahy reprodukovány na 0,1–0,3° všude,
kde je pohyb pravidelný.

## Reprodukce Pejchy et al. (2013): ansámbl {#pejcha2013}

Pejcha et al. (2013) N-částicovým kódem zintegrovali tisíce čtyřhvězd 2+2 se stejnými hmotnostmi
a počítali, jak často dvojhvězda dosáhne extrémní excentricity, $1-e<10^{-3}$, během 300
Kozaiových dob dvojhvězdy A. Sklony A a B k vnější dráze nebyly losovány izotropně, ale nastaveny
na mřížce v $(\cos i_A, \cos i_B)$, se 100 náhodnými nastaveními ostatních úhlů v každém bodě
mřížky; SWINGERsim je losuje rovnoměrně v $\cos i$ (izotropně), což je spojitý protějšek
rovnoměrné mřížky. Pro svou sadu F6 ($a_A$ = 0,573 au, $e_A$ = 0,3, $a_B$ = 0,351 au,
$e_B$ = 0,1, $a_{AB}$ = 6 au, $e_{AB}$ = 0,3, $t_K$ = 64,9 roku) zjistili 35,8 % pro A a 53,0 %
pro B – 4,8krát a 10krát více než u odpovídajících trojhvězd. Toto „zesílení“ ve čtyřhvězdách 2+2
je hlavním výsledkem článku. Vokrouhlický (2016, oddíl 4) zopakoval experiment se svým sekulárním
modelem ($a_{AB}$ = 5 au, 11 kyr) a zjistil 31,0 % a 48,1 %.

**Ansámbl** je ve SWINGERsimu množina běhů (**členů**), jejichž neznámé parametry jsou losovány
z **apriorního rozdělení** (rozdělení přijatelných hodnot); výsledkem je podíl členů s daným
výsledkem, se 68% Wilsonovým intervalem jako statistickou chybou.

**Tabulka 8.** Pejcha et al. (2013) a Vokrouhlický (2016): podíly dvojhvězd, které dosáhnou
$1-e<10^{-3}$. N-částicové testy jsou ruční (třída R), sekulární pomalé (třída own). Sloupce jsou definovány v oddílu [Jak číst tuto zprávu](#how-to-read).

| Test | Třída | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|---|
| `pejcha13-F6-fA` | R | 0,358 | 0,24 | −0,118 | 0,14 | PASS |
| `pejcha13-F6-fB` | R | 0,530 | 0,41 | −0,12 | 0,15 | PASS |
| `vok16-pejcha-fA` | own | 0,31 | 0,21 | −0,10 | 0,07 | **EXPLAINED** |
| `vok16-pejcha-fB` | own | 0,481 | 0,4025 | −0,0785 | 0,07 | **EXPLAINED** |

**Tabulka 9.** Ansámbly F6 podrobně. Sloupec 2 je počet členů, sloupce 3–4 podíly A a B, které
dosáhnou $1-e<10^{-3}$, s jejich statistickými chybami 1σ, sloupec 5 publikované hodnoty.

| Uspořádání | Členy | Podíl A | Podíl B | Publikováno |
|---|---|---|---|---|
| F6, N-částicově, 300 $t_K$ | 100 | 0,24 ± 0,04 | 0,41 ± 0,05 | 0,358, 0,530 (N-částicově) |
| F6, sekulárně, 300 $t_K$ | 400 | 0,215 ± 0,021 | 0,358 ± 0,024 | – |
| F6, sekulárně, 900 $t_K$ | 200 | 0,325 ± 0,033 | 0,505 ± 0,035 | – |
| uspořádání Vokrouhlického, sekulárně, 11 kyr (≈ 300 $t_K$) | 400 | 0,21 | 0,40 | 0,310, 0,481 (sekulárně) |
| totéž, sekulárně, 50 kyr | 400 | 0,36 | 0,57 | 0,310, 0,481 (sekulárně) |
| F6 s B nahrazenou hmotným bodem (trojhvězda), sekulárně | 400 | 0,020 ± 0,007 | – | asi 0,07 (N-částicově), 0,03 analyticky |
| tatáž trojhvězda, N-částicově | 100 | 0,06 ± 0,02 | – | asi 0,07 |

**Co souhlasí** (tabulky 8 a 9; v tabulce 9 je každý řádek jeden ansámbl a podíly jsou uvedeny
s chybami 1σ tam, kde byl ansámbl dost velký, aby je bylo možné uvést). (i) *Trojhvězda.*
S dvojhvězdou B nahrazenou hmotným bodem dává N-částicový backend 6 ± 2 % členů, které dosáhnou
$1-e<10^{-3}$, stejně jako Pejcha et al. (asi 7 %). Sekulární backend dává 2 % (analytický
kvadrupólový odhad je asi 3 %): v trojhvězdě s takto mírnou hierarchií zvyšují krátkoperiodické
členy vrcholy excentricity (Luo, Katz & Dong 2016). (ii) *Zesílení.* N-částicově dosáhne
dvojhvězda A ve čtyřhvězdě meze 4,0krát častěji než v trojhvězdě (článek: 4,8), B častěji než A
a žádný člen se tam nedostane mimo okno ZLK – hlavní výsledek obou článků. (iii) *Oba backendy
SWINGERsimu souhlasí navzájem* v rámci statistických chyb: pro týchž 100 členů 0,24 proti 0,24 (A)
a 0,41 proti 0,39 (B), přestože jednotlivé chaotické členy končí odlišně.

**Co se liší a proč.** Při nominálních 300 Kozaiových dobách zjišťuje SWINGERsim asi dvě třetiny
publikovaných podílů pro čtyřhvězdy, s oběma backendy a v obou uspořádáních, zatímco trojhvězda
souhlasí. Příčinu zužují čtyři testy. Nezávislý sekulární kód SecularMultiple (Hamers & Portegies
Zwart 2016) dává pro týchž 200 členů F6 hodnoty 0,19 (A) a 0,37 (B) při stejných multipólových
řádech a 0,225 a 0,385 se všemi řády, proti 0,23 a 0,37 ze SWINGERsimu pro tyto členy. Sekulární
rovnice SWINGERsimu integrované pevným krokem 0,1 roku metodou Rungeho–Kutty čtvrtého řádu, jako
u Vokrouhlického (2016), dávají stejné podíly jako adaptivní integrátor (0,275 a 0,30 pro 40
členů). A podíly s délkou běhu dále rostou (obrázek 4).

![**Obrázek 4.** Sada F6 Pejchy et al. (2013): podíl členů, které alespoň jednou dosáhly $1-e<10^{-3}$ (svislá osa), v závislosti na čase v Kozaiových dobách dvojhvězdy A (vodorovná osa, logaritmická); sekulární backend SWINGERsimu, 200 členů. Modře: dvojhvězda A; oranžově: dvojhvězda B. Tečkované vodorovné čáry jsou publikované podíly; čárkovaná svislá čára vyznačuje 300 $t_K$, uvedenou délku publikovaných běhů.](../../figures/part1/pejcha_time_curve.png)

Obrázek 4 ukazuje kumulativní podíl v závislosti na čase. Křivky jen stoupají, protože člen, který
mez dosáhl, zůstává započten. Nenasycují se: rostou zhruba logaritmicky, o 0,2 (A) a 0,3 (B) na
každý desetinásobek času, a publikovaných hodnot dosáhnou asi při 1500 $t_K$ (A) a 1200 $t_K$ (B).
Fang, Thompson & Hirata (2018) zjistili stejné nenasycení v širších soustavách 2+2 (asi 10 % na
desetinásobek času). Při takovém růstu závisejí publikovaná čísla silně na přesné délce běhů. Běh
sady Vokrouhlického (2016) do 50 kyr dává 0,36 a 0,57: publikované 0,31 a 0,48 leží mezi
SWINGERsimem při 11 kyr a při 50 kyr (obrázek 6 vpravo dole v [kapitole 4](#discs)). Tři nezávislé
integrace (sekulární SWINGERsim, N-částicový SWINGERsim s REBOUND, SecularMultiple) souhlasí
navzájem a integrační metoda sekulárního článku je jako příčina vyloučena. Rozdíl tedy spočívá
v efektivní délce publikovaných běhů (publikované podíly odpovídají 4–5násobku uvedené délky)
a možná ve vzorkování jejich mřížky sklonů, o níž Pejcha et al. poznamenávají, že „might not be
distributed uniformly“. N-částicový test vyhovuje na úrovni tří směrodatných odchylek; sekulární
test je vysvětlen (EXPLAINED) délkou publikovaných běhů. Tvar oblasti vysokých excentricit v rovině
$(\cos J_A, \cos J_B)$ souhlasí s obr. 11 Vokrouhlického (2016): diagonální pás pro A a široký pás
$|\cos J_B|\lesssim 0{,}6$ pro B.

## Srovnání s MSE {#mse}

**Reference.** MSE (Multiple Stellar Evolution; Hamers et al. 2021) je nezávislý program pro
hierarchické vícenásobné hvězdy. Byl sestaven z veřejného repozitáře a spuštěn na stejných
čtyřhvězdách 3+1 a 2+2 jako sekulární backendy SWINGERsimu, se stejným hamiltoniánem (členy MSE
typu binary–triplet, hexadekupólové a dotriakontupólové jsou vypnuty). Třetím srovnáním je
*osud*: čtyřhvězda 3+1, jejíž vnitřní dvojhvězdu (0,3 au) cykly ZLK dovedou k přetečení Rocheova
laloku – událost MSE „RLOF start“ proti prvnímu vyplnění Rocheova laloku v periastru ve SWINGERsimu,
se stejnými poloměry.

**Tabulka 10.** Sekulární backendy SWINGERsimu proti MSE (úroveň C, třída own). Sloupce jsou definovány v oddílu [Jak číst tuto zprávu](#how-to-read).

| Test | Třída | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|---|
| `mse-3p1-inner-emax` | own | 0,9707 | 0,970714 | +1,37·10⁻⁵ | 0,002 | PASS |
| `mse-3p1-inner-t-emax` | own | 49 325 let | 49 322,8 roku | −0,00 % | 0,2 % | PASS |
| `mse-3p1-mid-emax` | own | 0,2245 | 0,224461 | −3,86·10⁻⁵ | 0,002 | PASS |
| `mse-2p2-A-emax` | own | 0,9494 | 0,949492 | +9,2·10⁻⁵ | 0,002 | PASS |
| `mse-2p2-B-emax` | own | 0,7758 | 0,776227 | +0,000427 | 0,002 | PASS |
| `mse-3p1-rlof-time` | own | 380,03 roku | 380,1 roku | +0,02 % | 0,2 % | PASS |

**Výsledek (tabulka 10).** Maxima excentricity souhlasí na 10⁻⁴ nebo lépe a okamžik maxima vnitřní
dvojhvězdy čtyřhvězdy 3+1 na 2 roky z 49 kyr. U čtyřhvězdy 2+2 souhlasí první tři cykly ZLK
dvojhvězdy A a všechny cykly B, poté se oba běhy rozejdou ve fázi, jak se v chaotickém režimu
očekává. Dodatečné multipólové členy MSE mění $e_\mathrm{max}$ v těchto soustavách o méně než 0,001.
Počátek přetečení Rocheova laloku souhlasí na 0,07 roku z 380 let. Stojí za to znát dva praktické
rozdíly: MSE přechází na přímou N-částicovou integraci, jakmile se soustava stane dynamicky
nestabilní, kdežto sekulární backendy SWINGERsimu jen hlásí diagnostiku stability (uživatel pak
přejde na N-částicový backend); a MSE ztrácí argument periastra dráhy se sklonem přesně 0°, takže
srovnání používá 0,5–2°.

## Rozptyl a statistická teorie úniku {#scattering}

SWINGERsim umí sestavit setkání dvojhvězd a osamocených hvězd (rozptylové experimenty) a integrovat
je s REBOUND, dokud není výsledek zřejmý. Antognini & Thompson (2016, tabulka 3) uvádějí účinné
průřezy setkání dvou dvojhvězd (stejné hmotnosti, $a_1 = a_2$, kruhové dráhy, rychlost
v nekonečnu rovná kritické rychlosti $v_c$, při níž je celková energie nulová, srážkové parametry
do 7$a$). Stone & Leigh (2019) odvodili statistiku rozpadu rezonanční („demokratické“) soustavy tří
těles z ergodických úvah; SWINGERsim jejich rozdělení implementuje.

**Tabulka 11.** Rozptylové experimenty (třída R) a statistická teorie úniku (own). Sloupce jsou definovány v oddílu [Jak číst tuto zprávu](#how-to-read).

| Test | Třída | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|---|
| `scatter-hyperbola-deflection` | own | 173,548° | 173,548° | +1,4·10⁻¹²° | 10⁻⁶° | PASS |
| `at16-binbin-single-ionization` | R | 2,645 | 2,6873 | +1,60 % | 10 % | PASS |
| `at16-binbin-exchange-ionization` | R | 3,259 | 3,2305 | −0,88 % | 10 % | PASS |
| `at16-binbin-exchange` | R | 0,302 | 0,37165 | +23,06 % | 30 % | PASS |
| `escape-theory-energy-slope-L0` | own | −3,5 | −3,5026 | −0,0026 | 0,05 | PASS |
| `escape-theory-energy-slope-highL` | own | −4 | −4,02479 | −0,0248 | 0,1 | PASS |
| `escape-theory-mean-e-L0` | own | 0,8133 | 0,812817 | −0,000483 | 0,01 | PASS |
| `escape-theory-mean-e-highL` | own | 0,7 | 0,710776 | +0,0108 | 0,02 | PASS |
| `escape-theory-vs-nbody-light-escaper` | own | 0,769 | 0,75369 | −0,0153 | 0,07 | PASS |

**Výsledek (tabulka 11).** Nastavení jednoho hyperbolického průletu dává přesný úhel odklonu
problému dvou těles. Účinné průřezy (v jednotkách $\pi a^2$) pocházejí z 12 000 setkání, proti 10⁶
v článku. Prostá ionizace a výměna + ionizace souhlasí na 1,6 % a 0,9 %. Článek uvádí normování
$\pi(a_1^2+a_2^2)$, jeho čísla však souhlasí s $\pi a^2$ jedné dvojhvězdy; testy používají to druhé
(nejednoznačnost faktoru dva v článku). Výměny vycházejí o 23 % častější, než je publikováno; při 91
výměnných událostech to představuje 1,8 směrodatné odchylky a tolerance 30 % je statistickou chybou
takto malého počtu událostí. Stabilní trojhvězdy jsou vzácné v obou případech (0,008 proti
0,020 ± 0,002). Dvě z 12 000 setkání neměla v rámci integrační meze konečný výsledek (0,02 %)
a z účinných průřezů jsou vynechána. Sklony rozdělení energie unikající dvojhvězdy a její střední
excentricita reprodukují uzavřené tvary Stonea & Leigha (2019) na 0,0026–0,025 (sklony)
a 0,0005–0,011 (střední $e$). V rezonančních setkáních dvojhvězdy a osamocené hvězdy (dvojhvězda
1 + 0,6 M☉, hvězda 0,8 M☉) integrovaných N-částicově uniká nejlehčí hvězda v 75,4 % případů proti
76,9 % podle teorie. Hodnota 76,9 % není publikované číslo: jde o vlastní vyčíslení rozdělení
úniku podle Stonea & Leigha (2019) pro tyto hmotnosti ve SWINGERsimu, takže tento test srovnává
N-částicové integrace SWINGERsimu s jeho vlastní implementací teorie.

## Hierarchie s více než čtyřmi hvězdami {#larger-hierarchies}

SWINGERsim popisuje libovolnou hierarchii vnořených dvojic (2+2+1, 3+1+1, 2+(2+1), 2+2+2, 3+3
nebo vlastní vnoření) a je ověřen až do deseti hvězd. Větší hierarchie pokrývají dva druhy důkazů.

**Sekulárně proti N-částicově pro šestihvězdu (úroveň C).** Šestihvězda 2+2+2, jejíž dvojhvězda A
je skloněna o 60° k dráze čtyřhvězdy AB a jejíž dráha AB je skloněna o 50° k vnější dráze, se
integruje N-částicovým backendem (IAS15 + 1PN, 400 let, 98 s) a vnořeným sekulárním backendem
(páry dvojhvězd do oktupólového řádu; Hamers et al. 2015; Hamers & Portegies Zwart 2016; 1 s).
Pětihvězda 2+2+1 ověřuje cykly ZLK dráhy obou dvojhvězd buzené pátou hvězdou.

**V994 Her jako šestihvězda (úroveň D).** V994 Her je trojnásobně zákrytová šestihvězda (Zasche
et al. 2023): dvě dvojhvězdy A a B na vzájemné dráze a třetí dvojhvězda C na široké dráze kolem
nich. Tato dráha není známa: medián rozdělení periody, který Zasche et al. odhadují, je blízko
3000 let, kdežto jejich fotodynamický model používá pevnou kruhovou dráhu viděnou shora s periodou
asi 30 kyr. SWINGERsim pro ni bere 3000 let ($e$ = 0,3). S publikovanými hmotnostmi, poloměry
a ostatními drahami (všechny koplanární) dává N-částicový běh SWINGERsimu bez relativity a slapů
apsidální stáčení, které v každé dvojhvězdě způsobují ostatní hvězdy; s relativitou dává
relativistický podíl. Obojí je vyjádřeno jako podíl pozorovaných rychlostí (1,86, 3,66
a 1,68 °/rok).

**Tabulka 12.** Větší hierarchie. Sloupce jsou definovány v oddílu [Jak číst tuto zprávu](#how-to-read).

| Test | Třída | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|---|
| `sextuple-secular-emax` | own | 0,75207 | 0,751821 | −0,03 % | 2 % | PASS |
| `sextuple-secular-zlk` | own | 70,2 roku | 71,04 roku | +1,20 % | 5 % | PASS |
| `sextuple-nbody-emax` | R | 0,75207 | 0,75207 | +0,00 % | 0,2 % | PASS |
| `quint-out-zlk-emax` | R | 0,81272 | 0,8152 | +0,31 % | 1 % | PASS |
| `v994s-dyn-A` | R | 0,12 | 0,0999365 | −16,72 % | 25 % | PASS |
| `v994s-dyn-B` | R | 0,05 | 0,0519361 | +3,87 % | 25 % | PASS |
| `v994s-dyn-C` | R | 0 | 1,23·10⁻⁶ | +1,23·10⁻⁶ | 0,01 | PASS |
| `v994s-gr-A` | R | 0,10 | 0,0938886 | −6,11 % | 10 % | PASS |
| `v994s-gr-B` | R | 0,08 | 0,0748582 | −6,43 % | 10 % | PASS |
| `v994s-gr-C` | R | 0,06 | 0,0811201 | +35,20 % | 10 % | **EXPLAINED** |

**Výsledek (tabulka 12).** Vnořený sekulární backend dává největší excentricitu A s přesností
0,03 % a její periodu ZLK s přesností 1,2 % vůči N-částicovému běhu, za setinu výpočetní ceny.
Souhlasí i rozsahy ostatních drah (AB: $e$ = 0,190–0,327 sekulárně, 0,179–0,324 N-částicově;
vzájemný sklon B 6,6°–23,5° a 6,6°–23,2°). Samotný N-částicový výsledek je reprodukovatelný do
poslední číslice (`sextuple-nbody-emax`). Pro V994 Her jsou dynamické podíly 10,0 % (A), 5,2 % (B)
a 0 (C), proti publikovaným 12 %, 5 % a „zanedbatelný“. Rozdíl pro A leží v rámci vlivu nejistoty
vnější excentricity ($e_{AB}$ = 0,687 +0,050/−0,037), která mění rychlost až o 25 % – odtud
tolerance. Relativistické podíly jsou 9,4 % (A) a 7,5 % (B), proti uvedeným 10 % a 8 %. Pro C dává
SWINGERsim 8,1 %, článek 6 %. Tento rozdíl je vysvětlen (EXPLAINED): rychlost 1PN spočtená ručně
z vlastní tabulky 3 článku (1,81 + 1,077 M☉, $P$ = 1,960 d, $e$ = 0,189) je 0,136 °/rok, což je
8,1 % z pozorovaných 1,68 °/rok – hodnota SWINGERsimu. Hodnota 6 % je patrně přehlédnutím v článku
(samotná hmotnost primární složky dává 6 %); A a B, spočtené stejným způsobem, souhlasí.

Stabilita větších hierarchií se posuzuje algebraickými kritérii použitými na každou vnořenou
trojici a N-částicovým testem duchem ([kapitola 9](#stability)); neuronové klasifikátory byly
natrénovány jen na čtyřhvězdách a SWINGERsim je používá na každou dílčí čtyřhvězdu 2+2 a 3+1, což je
přiblížení mimo jejich trénovací množinu.

## Okolí: plyn, slapy Galaxie, procházející hvězdy a gravitační vlny {#environment}

Na vícenásobnou hvězdu působí zvenčí nebo během velmi dlouhých dob čtyři jevy.

- **Parametrická plynná fáze** (pro mladé soustavy) tlumí velkou poloosu, excentricitu a sklon
  k plynu s předepsanými časovými škálami, $\dot a = -a/\tau_a$, $\dot e = -e/\tau_e$, a stáčí
  normálu dráhy podle $\tan(i/2) = \tan(i_0/2)\,e^{-t/\tau_i}$ (exponenciální migrace, např.
  Papaloizou & Larwood 2000). Testy ověřují, že zintegrované dráhy těmto definicím odpovídají.
- **Slapy Galaxie** od disku (Heisler & Tremaine 1986) budí ZLK podobnou oscilaci velmi širokých
  drah; pro dráhu skloněnou o $i$ k disku a bez radiálních členů je maximální excentricita
  $e_\mathrm{max} = \sqrt{1 - \tfrac54\cos^2 i}$ (Grishin & Perets 2022, rov. 6).
- **Procházející hvězdy** udělují dvojhvězdě impulz; pro vzdálený parabolický průlet udává změnu
  vektoru excentricity v prvním řádu vzorec Hamerse & Samsinga (2019, rov. 9).
- **Gravitační vlny** zmenšují a zakulacují těsné dráhy rychlostmi středovanými přes oběh podle
  Peterse (1964). Referencí úrovně D je pozorované zkracování dráhy dvojitého pulzaru PSR B1913+16
  (Weisberg & Huang 2016).

**Tabulka 13.** Okolí a gravitační vlny. Sloupce jsou definovány v oddílu [Jak číst tuto zprávu](#how-to-read).

| Test | Třída | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|---|
| `gas-migration-a` | own | 5,48812 au | 5,48812 au | −0,00 % | 10⁻⁴ % | PASS |
| `gas-damping-e` | own | 0,111565 | 0,111565 | +0,00 % | 10⁻⁴ % | PASS |
| `gas-alignment-i` | own | 5,63999° | 5,63943° | −0,000561° | 0,01° | PASS |
| `galactic-tide-emax-60` | own | 0,829156 | 0,829181 | +2,5·10⁻⁵ | 0,01 | PASS |
| `galactic-tide-emax-60-nbody` | R | 0,82916 | 0,82908 | −0,01 % | 1 % | PASS |
| `flyby-first-order-kick` | own | 0 | 0,00594917 | +0,00595 | 0,03 | PASS |
| `flyby-nbody-backend` | R | 0 | 0,00738403 | +0,00738 | 0,05 | PASS |
| `gw-hulse-taylor-pdot` | own | −2,40263·10⁻¹² | −2,4021·10⁻¹² | +0,02 % | 0,1 % | PASS |
| `gw-secular-circular-decay` | own | 1 | 1 | −0,00 % | 0,01 % | PASS |
| `gw-nbody-decay` | R | 1 | 1,0003 | +0,03 % | 1 % | PASS |

**Výsledek (tabulka 13).** Dráha o velikosti 2·10⁴ au skloněná o 60° k disku dosáhne
$e$ = 0,82918 v sekulárním backendu a 0,82908 v N-částicovém backendu, proti analytické hodnotě
0,82916; energie středované soustavy se zachovává na 10⁻¹⁰. Vzorec impulzu prvního řádu souhlasí
s přímou integrací parabolického průletu integrátorem IAS15 ve vzdálenosti periastra 20 $a$ na
0,6 % (0,2 % při 40 $a$); tentýž průlet integrovaný jako částice v N-částicovém backendu s ním
souhlasí na 0,7 %. Petersovy rychlosti dávají zkracování dráhy PSR B1913+16 (hmotnosti 1,438
a 1,390 M☉, $P$ = 0,323 d, $e$ = 0,617) s přesností 0,02 % vůči hodnotě obecné relativity;
sekulární backend sleduje analytické zmenšování kruhové dvojhvězdy s přesností lepší než 10⁻⁴,
N-částicový backend (v němž SWINGERsim uplatňuje ztrátu středovanou přes oběh mezi integračními
kroky) rychlost s přesností 3·10⁻⁴.
