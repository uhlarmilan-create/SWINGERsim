# Přenos hmoty a vývoj {#mass-transfer-and-evolution}

Hvězdy mění své hmotnosti: ztrácejí hmotu větrem, předávají ji průvodci, když vyplní svůj Rocheův
lalok, sdílejí obálku v kontaktních dvojhvězdách, pohltí průvodce ve společné obálce nebo s ním
splynou. Každá taková změna působí na všechny dráhy soustavy, které hvězdu obsahují. Tato kapitola
popisuje, jak SWINGERsim tyto procesy modeluje a jak váže vývoj hvězd (kapitola [Hvězdy](#stars))
na dynamiku (kapitola [Dynamika](#dynamics)).

## Vazba hvězdného vývoje na dynamiku {#coupling}

Hvězdy se mění na jaderných časových škálách (Myr–Gyr), ve většině soustav mnohem pomaleji než
sekulární dynamika. S `physics.evolution: true` je vazbou **operátorové štěpení** jako v kódu MSE
(Hamers et al. 2021):

- **Sekulární backendy.** Běh se rozdělí na úseky, v nichž se žádný poloměr hvězdy, apsidální
  konstanta, moment setrvačnosti (a při ztrátě hmoty ani hmotnost) nezmění o více než
  `run.evolution_max_dr` (výchozí 1 %). Uvnitř úseku se stavba hvězdy (poloměr, $k_2$, $r_g^2$,
  slapové zpoždění) drží na hodnotě ve středu úseku; stav drah se přenáší z úseku do úseku a rotace
  zachovávají svůj moment hybnosti, když se moment setrvačnosti změní. Na tepelně pulzující AGB, kde
  poloměr při každém pulzu dosahuje maxima jen na několik staletí, obsahuje mřížka úseků uzly drah
  MIST a zjemňuje se půlením, dokud se žádný poloměr mezi dvěma body mřížky nezmění o více než
  `evolution_max_dr`. Běhy na miliony až miliardy let se zpracovávají proudově po částech
  a redukují se na své extrémy, takže se vejdou do paměti.
- **N-částicový backend.** Stavba hvězdy se aktualizuje v každém výstupním kroku (poloměr pro srážky
  a Rocheovy laloky, Loveovo číslo, moment setrvačnosti a zpoždění pro slapové síly, rotace při
  stálém rotačním momentu hybnosti) a drží se na hodnotě ve středu kroku. Výstupní kroky
  (`run.n_output`) proto musejí být krátké ve srovnání s vývojem.

Běh končí na konci vývojové dráhy kterékoli hvězdy (`stellar_evolution_end`: zapálení helia u SSE,
opuštění AGB nebo konec hoření uhlíku u MIST), zpět v čase při zrodu hvězd (`stellar_birth`)
a – pokud se nesleduje přenos hmoty – když hvězda vyplní svůj Rocheův lalok (`roche_overflow`).

## Větry {#winds}

S `physics.mass_loss` se hmotnosti během běhu mění, a to z jednoho ze dvou zdrojů:

- `source: track` – hmotnost každé hvězdy sleduje její vývojovou dráhu. Hmotnosti MIST obsahují
  předpisy pro vítr z modelů MESA, na nichž je MIST založen (Vink et al. 2001 pro horké hvězdy
  v plné síle, Nugis & Lamers 2000 pro Wolfovy–Rayetovy hvězdy, Reimers 1975 na větvi obrů
  a Blöcker 1995 na AGB, u rotujících hvězd zesílené faktorem $(1-\Omega/\Omega_\mathrm{crit})^{-0{,}43}$;
  Choi et al. 2016). Dráhy SSE jsou spočteny bez větru; pro ně SWINGERsim přidává vítr kódu SSE
  (Hurley, Pols & Tout 2000, odd. 7.1): větší z hodnot
  $9{,}6\times10^{-15}(Z/0{,}02)^{1/2}R^{0{,}81}L^{1{,}24}M^{0{,}16}$ M☉ yr⁻¹ (Nieuwenhuijzen & de Jager 1990)
  nad 4000 L☉ a $4\times10^{-13}\,\eta\,LR/M$ M☉ yr⁻¹ (Reimers 1975) od Hertzsprungovy mezery dál,
  s $\eta$ = `reimers_eta` (výchozí 0,5); $R$, $L$, $M$ ve slunečních jednotkách.
- `source: rates` – konstantní rychlosti $\dot M$ pro každou hvězdu, pro parametrické studie
  a srovnávací testy.

Vítr je izotropní v soustavě hvězdy a odnáší hybnost hvězdy (*Jeansův mód*; Jeans 1924;
Hadjidemetriou 1963). Jeho účinek na dráhu X s celkovou hmotností $M_X$ měří index ztráty hmoty
(Veras et al. 2011, rov. 15)

$$ \Psi = \frac{|\dot M_X|}{n_X M_X}, $$

tj. poměr oběžné periody k časové škále ztráty hmoty dělený $2\pi$ ($n_X$ je střední pohyb).

- **N-částicový backend.** Hmotnosti se během integrace mění operátorem REBOUNDx `modify_mass`
  ($\dot m = m/\tau$ při nezměněných rychlostech), přičemž $\tau$ každé hvězdy se nastaví v každém
  výstupním kroku tak, aby hmotnost na konci kroku dosáhla své hodnoty. To je přesné pro libovolné
  $\Psi$, od adiabatického režimu ($\Psi \ll 1$: $aM$ konstantní, $e$ v průměru konstantní) po
  impulzní ($\Psi \gg 1$: poloha a rychlost se zachovají, zatímco hmotnost klesne; kruhová
  dvojhvězda, která naráz ztratí polovinu hmotnosti, se rozváže; Blaauw 1961; Hills 1983).
- **Sekulární backendy.** Přes oběh středované rovnice potřebují adiabatickou odezvu: každá dráha
  zachovává vektor excentricity a orientaci a $a_X M_X$ je konstantní (Hadjidemetriou 1963). Běh se
  odmítne, když $\Psi$ kterékoli dráhy překročí `physics.mass_loss.max_psi` (výchozí 0,05; přechod
  do neadiabatického režimu leží při $\Psi \sim 0{,}01$–0,1, Veras et al. 2011); ukládá se největší
  $\Psi$ každé dráhy. Cykly ZLK testovací částice s rušitelem proměnné hmotnosti pak souhlasí
  s N-částicovým backendem v časování na 0,2 %.
- **Rotace.** Vítr odnáší rotační moment hybnosti rychlostí $\dot J = \tfrac{2}{3}\dot M R^2\Omega$
  (tenká kulová slupka; Hurley et al. 2002, rov. 11), takže $J \propto M^{2/(3 r_g^2)}$.

**Dynamická nestabilita vývoje trojhvězd.** Ztráta hmoty z vnitřní dvojhvězdy rozšíří vnitřní
dráhu o faktor $M_{\mathrm{in},0}/M_\mathrm{in}$ a vnější jen o menší faktor
$M_{\mathrm{tot},0}/M_\mathrm{tot}$, takže hierarchie slábne a stabilní trojhvězda může překročit
mez stability (Perets & Kratter 2012; Toonen et al. 2022). Ve srovnávacím testu začne vnitřní
excentricita růst s odchylkou do 7 % od okamžiku, kdy adiabaticky se rozšiřující dráhy překročí
mez Mardlinga & Aarsetha (2001); samotný rozpad pak trvá ještě mnoho vnějších oběhů, jak zjistili
Toonen et al.

## Přetečení Rocheova laloku a stabilní přenos hmoty {#roche-lobe-overflow}

Hvězda dvojhvězdy ze dvou hvězd, která vyplní svůj Rocheův lalok ([Hvězdy](#roche-lobes)), předává
hmotu svému průvodci. S `physics.mass_transfer` (vyžaduje `physics.evolution` s hvězdami SSE nebo
MIST) sleduje SWINGERsim stabilní přenos algoritmem vývoje dvojhvězd (BSE) podle Hurleyho, Touta
& Polse (2002, odd. 2.6). Každá dvojhvězda ze dvou hvězd se pojímá jako kruhová – slapy těsné
dvojhvězdy před přetečením zakruhují – a její **historie** (hmotnosti, vzdálenost složek, stavy
hvězd) se spočte před dynamikou z vývojových drah hvězd a počáteční vzdálenosti. Dynamika pak tuto
historii sleduje (oddíl [Historie v dynamice](#histories)).

![**Obrázek 9.** Přenos hmoty, schematicky. Dárce (oranžový, vlevo) vyplňuje svůj Rocheův lalok;
proud plynu (tmavě oranžová čára) jej opouští vnitřním Lagrangeovým bodem L1 a dopadá na příjemce
(modrý, vpravo), buď přímo, nebo přes akreční disk (světle oranžová elipsa). Šedá poznámka nahoře
připomíná, že část přenesené hmoty může dvojhvězdu opustit. Hmotnosti obou hvězd, dráha i stavba
obou hvězd se mění společně.](../../figures/concepts/mass_transfer.png)

Obrázek 9 pojmenovává složky, které níže uvedené rovnice kvantifikují: rychlost, s jakou dárce ztrácí
hmotu, podíl, který si příjemce ponechá, moment hybnosti odnášený hmotou, která odchází, a odezvu
obou hvězd.

**Rychlost** (Hurley et al. 2002, rov. 58–63). Dárce o poloměru $R$, který přeplňuje svůj lalok
$R_L$, přenáší

$$ \dot M = 3\times10^{-6}\,[\min(M, 5)]^2\,[\ln(R/R_L)]^3\ \mathrm{M_\odot\,yr^{-1}} $$

($M$ v M☉), takže zůstává těsně nad svým lalokem, zatímco přenos pohání jeho rozpínání nebo
zmenšování laloku. Rychlost je omezena tepelnou rychlostí $M/\tau_\mathrm{KH}$ dárců podobných
obrům (Hertzsprungova mezera a později; $\tau_\mathrm{KH} = 10^7\,M M'/(RL)$ yr s hmotností obálky
$M' = M - M_c$, kde $M_c$ je hmotnost jádra, ve slunečních jednotkách) nebo dynamickou rychlostí.

**Stabilita.** Přenos je dynamicky nestabilní a vede ke společné obálce, když poměr hmotností
$q = M_\mathrm{donor}/M_\mathrm{accretor}$ překročí kritickou hodnotu $q_\mathrm{crit}$
(tabulka 8).

**Tabulka 8.** Kritické poměry hmotností $q_\mathrm{crit} = M_\mathrm{donor}/M_\mathrm{accretor}$,
nad nimiž je přenos dynamický (Hurley et al. 2002, odd. 2.6). *Dárce* je vývojová fáze dárce;
$x = 0{,}3$ je exponent vztahu hmotnost–poloměr pro obry a $M_c/M$ podíl hmotnosti jádra.

| Dárce | $q_\mathrm{crit}$ | Poznámka |
|---|---|---|
| hlavní posloupnost pod 0,7 M☉ (hluboce konvektivní) | 0,695 | |
| hlavní posloupnost nad 0,7 M☉, hoření helia v jádře | 3,0 | `q_crit_ms`; např. 1,6 v binary_c (de Mink et al. 2007) |
| Hertzsprungova mezera | 4 | |
| větev obrů, AGB | $(1{,}67 - x + 2(M_c/M)^5)/2{,}13$ | rov. 57 |

S `stability: zeta` se dárci-obři posuzují místo toho podle adiabatické odezvy kondenzované
polytropy (Hjellming & Webbink 1987; Soberman et al. 1997, rov. 60–61),

$$\zeta_\mathrm{ad} = \frac{2}{3}\frac{m_c}{1-m_c} - \frac{1}{3}\frac{1-m_c}{1+2m_c} - 0{,}03\,m_c + \frac{0{,}2\,m_c}{1+(1-m_c)^{-6}},\qquad m_c = \frac{M_c}{M},$$

kde $\zeta_\mathrm{ad} = d\ln R/d\ln M$ dárce na dynamické časové škále; přenos je nestabilní, když
$\zeta_\mathrm{ad}$ klesne pod exponent Rocheova laloku $\zeta_L = d\ln R_L/d\ln M$ pro skutečné
podíly akrece.

**Akrece** (rov. 64–65). Příjemce na hlavní posloupnosti, v Hertzsprungově mezeře nebo hořící helium
v jádře přijme nejvýše $10\,\tau_{\dot M}/\tau_\mathrm{KH,2}$ přenesené hmoty
($\tau_{\dot M} = M_2/\dot M$, $\tau_\mathrm{KH,2}$ je jeho tepelná časová škála); zbytek odchází se
specifickým dráhovým momentem hybnosti dárce. Příjemci-obři přijmou vše. S `accretion: fractions`
jsou níže uvedené podíly `alpha` a `beta` konstantní.

**Dráha.** Z hmoty $|dM_d|$ ztracené dárcem odchází podíl $\alpha$ z okolí dárce, podíl $\beta$
z okolí příjemce (izotropní reemise) a $1-\alpha-\beta$ je akreováno. S dráhovým momentem hybnosti
$J = M_d M_a\sqrt{Ga/M}$ kruhové dvojhvězdy ($M_d$, $M_a$ jsou hmotnosti dárce a příjemce,
$M = M_d + M_a$) platí

$$\frac{d\ln a}{dM_d} = 2\,\frac{\alpha M_a^2 + \beta M_d^2}{M M_d M_a} + \frac{\alpha+\beta}{M}
- \frac{2}{M_d} + \frac{2(1-\alpha-\beta)}{M_a}.$$

To reprodukuje uzavřené tvary Sobermana et al. (1997) a Postnova & Jungelsona (2014) na $10^{-12}$:
konzervativní přenos zachovává $a M_d^2 M_a^2$, Jeansův mód $aM$.

**Odezva dárců a příjemců.** Každá hvězda nese svou hmotnost, svůj věk na vývojové dráze, kterou
sleduje, a hmotnost této dráhy (Hurley et al. 2000, odd. 7.1; 2002, odd. 2.6.6). Hvězda hlavní
posloupnosti nebo hvězda Hertzsprungovy mezery, jejíž nová hmotnost by nezmenšila její heliové
jádro, sleduje dráhu své současné hmotnosti ve stejném vývojovém stavu (EEP u MIST, zlomek fáze
u SSE): dárce se se svou hmotností zmenšuje a mezeru překonává na delší časové škále své nové
hmotnosti, příjemce se zjasňuje. S `rejuvenation: hurley` (výchozí) promíchá příjemce
s konvektivním jádrem (nad 1,25 M☉) nebo plně konvektivní příjemce (pod 0,35 M☉) čerstvý vodík do své
oblasti hoření a omladí se dále, $t' = (M/M')\,(t'_\mathrm{MS}/t_\mathrm{MS})\,t$, kde $t$ a $t'$
jsou věky před a po, $M$, $M'$ hmotnosti a $t_\mathrm{MS}$, $t'_\mathrm{MS}$ doby setrvání na hlavní
posloupnosti (Hurley et al. 2000). Příjemce proto vypadá mladší a méně vyvinutý než osamocená hvězda
jeho hmotnosti. Obři si ponechávají dráhu hmotnosti, kterou měli při opuštění mezery (jejich stavba
sleduje jádro), a jejich poloměr sleduje hmotnost obálky s $\zeta_\mathrm{eq} = -x + 2(M_c/M)^5$,
$x = 0{,}3$ (rov. 56). Heliové jádro se nikdy nezmenšuje.

**Události**, které ukončují historii: `common_envelope` (pokud se společná obálka nesleduje),
`contact` (příjemce vyplní svůj vlastní lalok, pokud se kontaktní dvojhvězdy nesledují),
`donor_stripped` (dárce je obnažen na své heliové jádro; heliové hvězdy se nesledují).

**Kontrola.** Pro příklad algolidy z BSE (2,9 + 0,9 M☉, 2,7 d) končí SWINGERsim s hmotností dárce
podle BSE (0,404 M☉; BSE v COSMIC, Breivik et al. 2020: 0,401; článek 0,42); akreovaná hmotnost je
nejméně jistým číslem tohoto příkladu a publikovaná čísla nejsou navzájem konzistentní (viz
*Testování, část I*).

## Historie v dynamice {#histories}

V dynamice sledují hmotnosti obou hvězd dvojhvězdy s přenosem hmoty její historii. Vzdálenost složek
se v každém kroku operátorového štěpení (sekulární backendy) nebo výstupním kroku (N-částicový:
relativní dráha $\mathbf{r} \to g\mathbf{r}$, $\mathbf{v} \to \mathbf{v}/\sqrt g$ při pevném vektoru
excentricity a fázi) škáluje o změnu nad rámec adiabatické odezvy $aM$ = konst., kterou backendy
uplatňují samy. Vnější dráhy reagují na hmotu opouštějící dvojhvězdu adiabaticky. Během přenosu je
poloměr dárce roven poloměru jeho Rocheova laloku (pro zákryty, slapy a diagnostiky) a přetečení
hvězdy s přenosem běh neukončí. Přemodelované hvězdy přenášejí svůj poloměr, zářivý výkon a teplotu
do dynamiky a do světelných křivek. Větry a přenos se kombinují: hvězdy ztrácejí hmotu větrem před
přenosem, během něj i po něm a dvojhvězda se pro větrnou část rozšiřuje podle $aM$ = konst.

### Přenos spuštěný dynamikou {#dynamic-onset}

Historie se počítají pro kruhové dvojhvězdy s nakonfigurovanými vzdálenostmi složek. Dynamika však
může přes lalok vyhnat hvězdu, s níž žádná historie nepočítá: cykly ZLK zvyšují excentricitu vnitřní
dvojhvězdy, dokud se hvězda v periastru nedotkne svého laloku (Eggleton & Kisseleva-Eggleton 2006;
Fabrycky & Tremaine 2007). Slapy pak dvojhvězdu rychle zakruhují; kódy pro trojhvězdy předpokládají,
že dráha je při začátku přenosu kruhová (Toonen et al. 2016). S
`physics.mass_transfer.dynamic_onset` ukončí SWINGERsim běh při dotyku – při excentricitě $e_c$, pro
niž $a(1-e_c)\,r_L = R$, kde $r_L$ je lalok v jednotkách vzdálenosti složek – a pokračuje se
zakruhovanou dvojhvězdou: při stálém dráhovém momentu hybnosti $a_c = a(1-e_c^2)$ (volba
`angular_momentum`; Hut 1981), nebo v periastru $a_c = a(1-e_c)$ (volba `periastron`, dárce zůstává
na svém laloku). Dvojice si zachová těžiště, rovinu a fázi, její hvězdy rotují synchronně a odtud
začíná nová historie přenosu hmoty. V trojhvězdě vyhnané na $e_c = 0{,}89$ najdou N-částicový
i sekulární backend dotyk s odchylkou do 1,7 % v čase a stejnou kruhovou dráhu.

## Polodotykové a kontaktní dvojhvězdy {#semidetached-contact}

**Předepsaná změna periody.** Algolidy a dvojhvězdy typu W UMa mění své periody přenosem hmoty
a ztrátou momentu hybnosti; jejich minima pak leží na parabole, $O - C = \tfrac12 \dot P P E^2$ ($E$
je číslo cyklu; Kwee 1958; Hilditch 2001). Když příčinu není třeba modelovat, předepíše
`orbits.<pair>.Pdot` dvojhvězdě ze dvou hvězd uvnitř N-částicového běhu konstantní $\dot P$
(d yr⁻¹), $P(t) = P + \dot P\,t$. Časová osa se rozdělí na intervaly délky $\Delta t$ a tam, kde
integrace překročí hranici (kterýmkoli směrem), se relativní dráha přeškáluje na periodu nového
intervalu, $a \to a(P_\mathrm{new}/P_\mathrm{old})^{2/3}$, při pevném vektoru excentricity a fázi
a s nezměněným těžištěm dvojice. Okamžiky zákrytů se pak od paraboly odchylují nejvýše
o $(\dot P/P)\Delta t^2/8$, což se drží pod 0,05 s. `Pdot_mode: conservative` mění hmotnosti
s periodou jako při konzervativním přenosu, $\dot P/P = 3\dot M_1 (M_1 - M_2)/(M_1 M_2)$. Všechny
pozorovatelné veličiny i třetí těleso vidí tutéž driftující dráhu, takže parabolu a efekt světelného
času lze prokládat společně ([Prokládání](#free-parameters)). Kvadratický koeficient spočtených minim
je roven $\tfrac12 P\dot P$ s přesností $5\times10^{-5}$. Tato volba existuje jen v N-částicovém
backendu.

**Kontaktní dvojhvězdy.** Když obě hvězdy vyplní své laloky, nechá je
`physics.mass_transfer.contact` sdílet společnou obálku (Lucy 1968), místo aby historii ukončil: méně
hmotná hvězda zůstává na svém laloku a přes obálku (konzervativně) zásobuje primární složku, zatímco
magnetizovaný vítr obálky odnáší moment hybnosti,

$$\dot J = -4{,}9\times10^{41}\,\frac{R_1^2 M_1 + R_2^2 M_2}{P}\ \mathrm{g\,cm^2\,s^{-1}\,yr^{-1}}$$

($R$, $M$ ve slunečních jednotkách, $P$ ve dnech; Gazeas & Stępień 2008, rov. 19). Dráha se zmenšuje
a poměr hmotností klesá; dvojhvězda splyne, když dosáhne Darwinovy nestability ($q < q_\mathrm{min}$,
výchozí 0,09; Rasio 1995) nebo když přeplní svou vnější kritickou plochu přes L2 (události
`contact_merger_darwin`, `contact_merger_l2`). Oscilace tepelné relaxace a přenos energie mezi
složkami nad rámec jedné společné teploty se nemodelují.

## Přenos hmoty z třetí hvězdy {#tertiary-transfer}

Hvězda, která obíhá kolem dvojhvězdy ze dvou hvězd (například obří třetí hvězda v HD 181068)
a vyplní svůj Rocheův lalok, přenáší hmotu na tuto dvojhvězdu (`physics.mass_transfer.tertiary`,
standardně zapnuto, sekulární backend). SWINGERsim spočte historii tří hvězd před dynamikou, stejně
jako u dvojhvězdy ze dvou hvězd: dárce vidí vnitřní dvojhvězdu jako hmotný bod, takže jeho lalok je
lalokem dvojhvězdy s $q = M_3/(M_1+M_2)$ v periastru vnější dráhy (Eggleton 1983), a rychlost a její
meze se řídí BSE (Hurley et al. 2002; Hamers et al. 2022). Proud vytvoří cirkumbinární disk, když je
vnitřní vzdálenost složek menší než cirkularizační vzdálenost $a_c = a_\mathrm{peri}(1 - R_L)^2$, kde
$R_L$ je poloměr Rocheova laloku dárce v jednotkách vnější vzdálenosti (Eggleton 1983). Leigh et al.
(2020, rov. 5) píší $a_c = a_\mathrm{out}(1 - R_L)^2$ s vnější velkou poloosou; SWINGERsim ji
vyhodnocuje v periastru $a_\mathrm{peri}$ vnější dráhy, kde dárce přetéká. S diskem se akreuje 90 %
hmoty (dvě třetiny hmotnější vnitřní hvězdou), jinak 20 % (Hamers et al. 2022). Zbytek opouští
vnitřní dvojhvězdu izotropně; pro vnější dráhu jej reemituje příjemce (Soberman et al. 1997) a jeho
vymrštění stojí vnitřní dráhovou energii (de Vries et al. 2014):

$$\mathrm{d}\ln a_\mathrm{in} = \frac{\mathrm{d}M_1}{M_1} + \frac{\mathrm{d}M_2}{M_2} - \frac{2 M_\mathrm{in}\,\mathrm{d}M_\mathrm{ej}}{\alpha\lambda\, M_1 M_2},$$

kde $M_1$, $M_2$ jsou vnitřní hmotnosti, $M_\mathrm{in}$ jejich součet, $\mathrm{d}M_\mathrm{ej}$
vymrštěná hmota a $\alpha\lambda$ = 5 účinnost (de Vries et al. našli z hydrodynamiky 3–7);
`inner_orbit: angular_momentum` místo toho nechá akreovaný plyn přinést specifický moment hybnosti
$h\sqrt{G M_\mathrm{in} a_\mathrm{in}}$. Historie se zastaví při společné obálce trojhvězdy, když
vnitřní hvězda vyplní svůj vlastní lalok, když se trojhvězda stane dynamicky nestabilní (Mardling
& Aarseth 2001), když je dárce obnažen na jádro nebo na konci jeho vývojové dráhy. Nemodeluje se:
odpor proudu na vnitřní dráhu nad rámec energetického formalismu, změna vnější excentricity během
přenosu a epizody v každém periastru excentrické vnější dráhy (rychlost je středována přes oběh).
Mimo sekulární backend nebo bez této volby přetečení třetí hvězdy běh ukončí; ansámbly zaznamenávají,
která hvězda přetekla první a zda je jejím průvodcem hvězda, nebo dvojhvězda.

## Společná obálka {#common-envelope}

Když se přenos stane dynamicky nestabilním, obálka dárce pohltí obě hvězdy. S
`physics.mass_transfer.common_envelope` spočte SWINGERsim výsledek energetickým formalismem Hurleyho,
Touta & Polse (2002, odd. 2.7.1; podle Webbinka 1984 a Touta et al. 1997), místo aby běh zastavil:

- vazebná energie obálek (rov. 69),
  $E_\mathrm{bind} = -\frac{G}{\lambda}\left(\frac{M_1 M_\mathrm{env,1}}{R_1} + \frac{M_2 M_\mathrm{env,2}'}{R_2}\right)$,
  se strukturním parametrem $\lambda$ (`lam`, výchozí 0,5), přičemž průvodce na hlavní posloupnosti
  se počítá celý jako jádro ($M_{c,2}' = M_2$, $M'_\mathrm{env,2} = 0$);
- jádra se spirálovitě přibližují a předávají obálce dráhovou energii s účinností
  $\alpha_\mathrm{CE}$ (`alpha`, výchozí 1; rov. 70–72),
  $E_\mathrm{bind} = \alpha_\mathrm{CE}\left(-\frac{GM_{c,1}M_{c,2}'}{2a_f} + \frac{GM_{c,1}M_{c,2}'}{2a_i}\right)$,
  což dává výslednou vzdálenost $a_f$ z počáteční $a_i$;
- kdyby jádro při $a_f$ vyplnilo svůj Rocheův lalok, jádra splynou tam, kde k tomu dojde, a splynulá
  hvězda si ponechá tu část obálky, kterou uvolněná energie nedokázala rozvázat (rov. 73–77).
  Poloměry jader: poloměr bílého trpaslíka pro degenerovaná heliová jádra a jádra AGB, jinak poloměr
  heliové hvězdy nulového stáří (Hurley et al. 2000).

Obálka odchází izotropně během `duration_yr`. V N-částicovém backendu si těžiště dvojhvězdy zachová
polohu a rychlost, zatímco hvězdy získají nové hmotnosti a kruhovou dráhu se vzdáleností $a_f$;
vnější dráhy pocítí impulzní ztrátu hmoty a rozšíří se, stanou se excentrickými ($e' = f/(1-f)$ pro
kruhovou dráhu, která ztratí podíl $f$ své hmotnosti) nebo nevázanými ($f \ge 1/2$; Hadjidemetriou
1963; Veras et al. 2011). Sekulární backendy se na konci obálky zastaví, protože impulzní odezva
závisí na oběžné fázi. Po obálce je jádro hmotným bodem stálé hmotnosti (pozůstatky se nevyvíjejí);
průvodce pokračuje po dráze své současné hmotnosti. Výsledek je ověřen na publikovaných příkladech
BSE a na analytické impulzní odezvě (viz *Testování, část I*). Třetí hvězda, která vyplní svůj lalok
na excentrické vnější dráze, může vést také ke společné obálce trojhvězdy (Glanz & Perets 2021;
Hamers et al. 2022), která běh ukončí.

## Splynutí {#mergers}

S `physics.mergers` běh pokračuje i poté, co dvě hvězdy splynou – srážkou v N-částicovém backendu
nebo společnou obálkou či kontaktní dvojhvězdou, která skončí splynutím –, s jednou splynulou
hvězdou a zbývajícími drahami (ze čtyřhvězdy se stane trojhvězda). Srážky ve vázaných soustavách
několika těles jsou pomalé (relativní rychlost zhruba úniková rychlost), takže hvězdy spíše splynou,
než by se navzájem zničily. Splynutí vymrští podíl $\phi = 0{,}3\,q/(1+q)^2$ celkové hmotnosti,
$q = M_2/M_1 \le 1$ (čelní srážky hvězd hlavní posloupnosti, Glebbeek et al. 2013: 2–7,5 %; Lombardi
et al. 2002: nejvýše asi 8 %); splynutí ve společné obálce si ponechá hmotnost z energetického
formalismu. Produkt se nachází v těžišti obou hvězd s jejich celkovou hybností; hmota odchází
naráz, takže vnější dráhy pocítí impulzní ztrátu hmoty. Do jádra produktu se vmíchá čerstvý vodík:
jeho poměrný věk na hlavní posloupnosti je $f = r\,(M_1 f_1 + M_2 f_2)/(M_1 + M_2)$, s $r$ od 0,81
(M ≤ 6 M☉) do 0,64 (M ≥ 60 M☉) (Glebbeek & Pols 2008; Schneider et al. 2016) nebo 0,1 podle BSE
(Hurley et al. 2002). Produkt vzniká blízko meze rozpadu rotací, a proto začíná s podílem
`spin_fraction` (výchozí 0,3) kritické rychlosti rotace a rychle se brzdí (Schneider et al. 2019).

## Magnetické brzdění a aktivita {#magnetic}

**Magnetické brzdění.** S `physics.magnetic_braking` (vyžaduje `evolution` a `rotation`) se rotace
hvězd s konvektivními obálkami zpomaluje momentem síly podle Matta et al. (2015),

$$T = -T_0\left(\frac{\tau_\mathrm{cz}}{\tau_\mathrm{cz,\odot}}\right)^p\left(\frac{\Omega}{\Omega_\odot}\right)^{p+1}
\ \text{(nenasycený režim)},\qquad T = -T_0\,\chi^p\,\frac{\Omega}{\Omega_\odot}\ \text{(nasycený režim)},$$

s $T_0 = 9{,}5\times10^{30}\,\mathrm{erg}\,(R/\mathrm{R}_\odot)^{3{,}1}(M/\mathrm{M}_\odot)^{0{,}5}$,
$\chi = 10$, $p = 2$, $\Omega_\odot = 2{,}6\times10^{-6}$ s⁻¹, nasycením pod Rossbyho číslem
$\mathrm{Ro}_\odot/\chi$ ($\mathrm{Ro} = 1/(\Omega\tau_\mathrm{cz})$) a konvektivním obrátkovým časem
podle Cranmera & Saara (2011),
$\tau_\mathrm{cz} = 314{,}24\,\exp[-T_\mathrm{eff}/1952{,}5\,\mathrm{K} - (T_\mathrm{eff}/6250\,\mathrm{K})^{18}] + 0{,}002$ d
(12,9 d pro Slunce, téměř nula nad 6500 K). Uvnitř každého úseku je moment setrvačnosti konstantní
a $d\Omega/dt = T/I$ se integruje přesně. Se slapy je brzděná hvězda v těsné dvojhvězdě opět
roztáčena svým průvodcem a dvojhvězda ztrácí dráhový moment hybnosti: magnetické brzdění zmenšuje
těsné dvojhvězdy, což je jedna z cest ke kontaktním dvojhvězdám v trojhvězdách. Rotace Slunce ve
slunečním věku a Skumanichův zákon $\Omega \propto t^{-1/2}$ (Skumanich 1972) jsou reprodukovány; pro
hvězdy pod 0,6 M☉ je moment síly známě příliš slabý (Douglas et al. 2017).

**Applegateův mechanismus.** Magneticky aktivní hvězda, jejíž kvadrupólový moment se během cyklu
aktivity mění, moduluje periodu své dvojhvězdy (Applegate 1992). `swingersim applegate` porovná
energii, kterou modulace potřebuje, s energií, kterou hvězda vyzáří za jednu periodu modulace (model
tenké slupky podle Tiana et al. 2009 ve tvaru používaném Völschowem et al. 2016,
$\Delta P/P = 2\pi\,(O-C)/P_\mathrm{mod}$ s plnou amplitudou O−C), a udá největší poloamplitudu O−C,
kterou může aktivita vyvolat. Změna okamžiků minim větší než tato nemůže být způsobena aktivitou;
dvouzónové modely Völschowa et al. potřebují o jeden až dva řády více energie, takže mez modelu tenké
slupky je velkorysá. S `physics.activity` se cyklus stává součástí dopředného modelu: zákryty
dvojhvězdy se v syntetických okamžicích a minimech posunou o $K\sin(2\pi t/P_\mathrm{mod} + \varphi)$.

## Meze {#mass-transfer-limits}

- Historie přenosu hmoty předpokládají kruhové dvojhvězdy; excentrický přenos vstupuje jen přes
  začátek spuštěný dynamikou s okamžitou cirkularizací.
- Heliové hvězdy, bílí trpaslíci, neutronové hvězdy a černé díry se nevyvíjejí; jádro po společné
  obálce je hmotným bodem stálé hmotnosti.
- Větry jsou izotropní a opouštějí soustavu; akrece větru průvodcem se nemodeluje.
- Rychlosti a stabilita se řídí předpisy BSE a nesou jejich nejistoty (kritické poměry hmotností,
  $\alpha_\mathrm{CE}\lambda$, účinnost akrece).

## Implementace {#implementation-mass-transfer}

**Externí knihovna a data.** Změna hmotnosti v N-částicovém backendu používá operátor REBOUNDx
`modify_mass` (Tamayo et al. 2020); větry MIST jsou součástí vývojových drah MIST (Choi et al. 2016).
**Vlastní kód.** Vše ostatní bylo pro SWINGERsim napsáno podle citovaných prací: operátorové štěpení
hvězdného vývoje a dynamiky s adaptivními úseky a proudovým zpracováním; předpis větru SSE
integrovaný podél drah SSE; adiabatická odezva a kontrola $\Psi$ v sekulárních backendech; algoritmus
přenosu hmoty BSE (rychlosti, meze, kritické poměry hmotností, účinnost akrece, rovnice
$d\ln a/dM_d$, přemodelování dárců a příjemců, omlazení); kritérium stability $\zeta$; začátek
přenosu spuštěný dynamikou s cirkularizací; předepsaná změna periody; kontaktní dvojhvězdy se ztrátou
momentu hybnosti a jejich kritéria splynutí; přenos hmoty z třetí hvězdy; energetický formalismus
společné obálky; splynutí; magnetické brzdění; a energetická kontrola Applegateova mechanismu
a modulace aktivitou.
