# Slapy a rotace {#tides-and-rotation}

Hvězda v těsné dvojhvězdě není bod. Průvodce přitahuje její bližší stranu silněji než vzdálenější,
takže se hvězda protáhne podél spojnice středů (*slapová výduť*); rotuje-li, je navíc zploštělá
na pólech (*rotační výduť*). Obě deformace přidávají ke gravitaci hvězdy kvadrupól, který stáčí
periastrum (apsidální stáčení) a váže rotaci hvězdy na dráhu. Vnitřní tření způsobuje, že slapová
výduť zaostává za průvodcem; toto zpoždění přenáší moment hybnosti mezi rotací a dráhou a disipuje
energii, takže se rotace synchronizují a dráhy se zakruhovávají a zmenšují. Ve vícenásobné hvězdě
mají slapy ještě druhou úlohu: rychlé apsidální stáčení, které v těsných dvojhvězdách způsobují,
potlačuje cykly ZLK ([Hierarchie a dráhy](#zlk-cycles)).

![**Obrázek 8.** Slapy, schematicky. Vlevo: hvězda (oranžová) deformovaná svým průvodcem (modrý);
černá přímka je osa slapové výdutě, která kvůli vnitřnímu tření zaostává za spojnicí středů (šedá
tečkovaná). Vpravo: vývoj pod vlivem slapového tření v jednotkách slapové časové škály (vodorovná
osa). Oranžově: rotační frekvence hvězdy dělená oběžnou frekvencí, rostoucí z 0,3 na 1
(synchronizace); modře: excentricita dráhy, klesající z 0,4 k 0 (cirkularizace). Synchronizace je
dosaženo několikrát rychleji než cirkularizace, protože hvězda uchovává ve své rotaci mnohem méně
momentu hybnosti než dráha.](../../figures/concepts/tides.png)

Obrázek 8 shrnuje pořadí dějů. Pravý panel je ilustrativní (libovolná jednotka času), jeho sdělení
však platí i pro skutečné dvojhvězdy: nejdřív se přizpůsobí rotace, později dráha. Rovnice, na nichž
je založen, jsou uvedeny níže.

## Rovnovážné slapy v N-částicovém backendu {#equilibrium-tides}

S `physics.tides: true` se hvězdy stávají rozlehlými tělesy. SWINGERsim používá efekt REBOUNDx
`tides_spin` (Lu et al. 2023), který implementuje model rovnovážných slapů Eggletona, Kiselevové
& Huta (1998) pro libovolný počet těles:

- každá hvězda je deformována svou vlastní rotací a slapovým polem **všech** ostatních hvězd
  a deformované hvězdy se přitahují odpovídajícími kvadrupólovými silami; tato konzervativní část
  stáčí periastra a způsobuje precesi rovin drah kolem rotačních os hvězd;
- vektor rotace každé hvězdy se vyvíjí pod vlivem momentů sil všech ostatních hvězd (například
  precese rotační osy skloněné k dráze);
- s časovým zpožděním $\Delta t > 0$ slapová výduť zaostává, energie se disipuje, dráhy se
  zakruhovávají a zmenšují a rotace se blíží pseudosynchronní.

Sílu deformace určuje potenciálové Loveovo číslo $2k_2$ a $R^5$, vývoj rotace moment setrvačnosti
$I = r_g^2 M R^2$ ([Hvězdy](#structure-constants)). `physics.rotation` (vyžaduje `physics.tides`)
udělí hvězdám rotaci: standardně synchronní s jejich dvojhvězdou (`stars.<name>.P_rot: sync`),
pseudosynchronní (`pseudosync`) nebo se zadanou periodou ve dnech, se sklonem rotační osy
`obliquity` (úhlem mezi rotací a dráhou) a azimutem `spin_azimuth`. Bez `physics.rotation` začínají
hvězdy bez rotace, takže zbývá jen slapová deformace (model, který použil Vokrouhlický 2016, jeho
obr. 7).

**Konvence časového zpoždění.** SWINGERsim očekává časové zpoždění `lag_s` (v sekundách) v konvenci
Mignarda (1979), Huta (1981) a Leconta et al. (2010). Rozbor zdrojového kódu `tides_spin`
a srovnání jeho rychlostí s těmito vzorci ukázaly, že stejné číslo předané REBOUNDx jako jeho
parametr $\tau$ dává rychlosti přesně dvojnásobné; SWINGERsim proto předává $\tau = \Delta t/2$
a výsledek kontrolují srovnávací testy s Lecontem et al. (2010) (viz *Testování, část I*). Pro
orientaci: konstantní časové zpoždění odpovídá slapovému faktoru kvality
$Q \approx 1/(2\,|\Omega - n|\,\Delta t)$ na slapové frekvenci, kde $\Omega$ je úhlová rychlost
rotace a $n$ oběžná úhlová frekvence.

**Zachování.** Bez disipace se energie včetně rotační a slapové energie zachovává na
$|dE/E| \sim 10^{-13}$ (bez relativity). S disipací energie fyzikálně klesá a běh místo ní kontroluje
celkový moment hybnosti drah a rotací.

## Klasické apsidální stáčení dvojhvězdy {#classical-apsidal-motion}

Pro každou dvojhvězdu ze dvou hvězd udává SWINGERsim stáčení apsid očekávané od dvojhvězdy samotné,
klasický slapový a rotační člen podle Sterna (1939):

$$\frac{\dot\omega}{n} = \sum_{i} k_{2,i}\left(\frac{R_i}{a}\right)^5
\left[15\,\frac{m_j}{m_i}\,f(e) + \left(1 + \frac{m_j}{m_i}\right)
\left(\frac{\Omega_i}{n}\right)^2 g(e)\right],$$

$$ f(e) = \frac{1 + \tfrac32 e^2 + \tfrac18 e^4}{(1-e^2)^5}, \qquad g(e) = \frac{1}{(1-e^2)^2}, $$

kde se sčítá přes obě hvězdy $i$ (s průvodcem $j$), $k_{2,i}$ je apsidální konstanta, $R_i$ poloměr,
$m_i$ hmotnost, $\Omega_i$ úhlová frekvence rotace (osa rotace kolmá k dráze), $a$ velká poloosa, $e$
excentricita a $n$ střední pohyb. První člen v závorkách je slapový příspěvek, druhý rotační;
relativistický člen ([Hierarchie a dráhy](#apsidal-motion)) se přidává zvlášť. Strmá závislost na
$(R/a)^5$ dělá z apsidálního stáčení sondu poloměrů a prostřednictvím $k_2$ i vnitřní stavby hvězd
(Claret & Giménez 2010).

Ve vícenásobné soustavě obsahuje naměřené stáčení navíc poruchu od ostatních hvězd; souhrn běhu
uvádí obojí a `swingersim apsidal` rozloží celek na slapovou, rotační, relativistickou a dynamickou
část a udá faktor, o který by se $k_2$ obou hvězd musela změnit, aby bylo dosaženo cílové rychlosti.
Publikovaná fotodynamická řešení často prokládají rychlost stáčení apsid jako volný parametr, místo
aby ji počítala ze stavby hvězd; rozdíl oproti SWINGERsimu pak měří napětí mezi naměřenou rychlostí
a modely hvězd, nikoli rozdíl mezi dvěma kódy.

**Pseudosynchronizace.** Na excentrické dráze nemohou slapy synchronizovat rotaci s dráhou v každém
bodě; přivedou ji k *pseudosynchronní* rychlosti, při níž přes oběh středovaný slapový moment síly
vymizí (Hut 1981, rov. 42) a která je rychlejší než střední pohyb (1,24krát při $e = 0{,}2$).
SWINGERsim ji používá pro `P_rot: pseudosync` a uvádí pseudosynchronní periodu mezi slapovými
časovými škálami níže.

**Rotace a $k_2$.** Rotující hvězda je koncentrovanější ke středu než statická. S
`physics.k2_rotation` se modelové konstanty opraví o $\delta\log k_2 = -0{,}87\,\lambda_s$,
$\lambda_s = 2\Omega^2R^3/(3GM) = 2V^2/(3gR)$ (Claret 1999; Claret & Giménez 2010 tuto rotační
opravu používají u pozorovaných dvojhvězd), s úhlovou rychlostí $\Omega$, rovníkovou rychlostí $V$,
povrchovým tíhovým zrychlením $g$, poloměrem $R$ a hmotností $M$ každé hvězdy v daném okamžiku,
ve všech backendech.

## Slapový vývoj: model konstantního časového zpoždění {#tidal-evolution}

Sekulární backendy a analytické časové škály používají teorii konstantního časového zpoždění
(Mignard 1979; Hut 1981) ve tvaru podle Leconta et al. (2010, rov. 1–3) pro rotační osy kolmé
k dráze. S potenciálovým Loveovým číslem $2k_{2,i}$ a časovým zpožděním $\Delta t_i$ hvězdy $i$ je

$$ K_i = \tfrac32\,(2k_{2,i})\,\Delta t_i\,\frac{G m_j^2}{R_i}\left(\frac{R_i}{a}\right)^6 n^2 , $$

$$ \frac{da}{dt} = \frac{4a^2}{G m_1 m_2}\sum_i K_i\left[N(e)\frac{\Omega_i}{n} - N_{a1}(e)\right], \qquad
\frac{de}{dt} = \frac{11\,a\,e}{G m_1 m_2}\sum_i K_i\left[\Omega_e(e)\frac{\Omega_i}{n} - \frac{18}{11}N_{e1}(e)\right], $$

$$ \frac{d\Omega_i}{dt} = -\frac{2K_i}{C_i\,n}\left[\Omega(e)\frac{\Omega_i}{n} - N(e)\right], \qquad C_i = r_{g,i}^2 m_i R_i^2 , $$

kde $N$, $N_{a1}$, $N_{e1}$, $\Omega_e$ a $\Omega(e)$ jsou polynomiální funkce excentricity podle
Leconta et al. (2010) (tytéž funkce se objevují u Huta 1981, rov. 9–11), $C_i$ moment setrvačnosti
a ostatní symboly jako výše. Rovnice pro rotaci plyne z rovnováhy slapového momentu síly se změnou
dráhového momentu hybnosti; pro kruhovou dráhu přechází v rovnici Murraye & Dermotta (1999,
rov. 4.159). Sekulární backendy integrují tyto rovnice pro každou dráhu dvou osamocených hvězd spolu
s gravitační dynamikou.

Pro každou dvojhvězdu uvádí souhrn běhu časové škály na začátku: cirkularizace $e/|\dot e|$,
zmenšování dráhy $a/|\dot a|$, dobu, za niž se odchylka rotace každé hvězdy od pseudosynchronismu
zmenší $e$-krát, a pseudosynchronní periodu. Realistická zpoždění hvězd dávají časové škály milionů
let a delší, mimo dosah přímých N-částicových integrací; v N-částicových bězích je disipace užitečná
hlavně pro silně excentrické fáze (maxima ZLK, kde zmenšuje dráhy; Fabrycky & Tremaine 2007) a pro
testy.

## Disipace ze stavby hvězdy {#dissipation-from-structure}

Časové zpoždění je nejhůře známým vstupem slapového vývoje. S `physics.tidal_lag: hurley2002` je
SWINGERsim odvodí pro každou hvězdu dvojhvězdy ze dvou hvězd z její stavby podle Hurleyho, Touta
& Polse (2002, odd. 2.3), místo pevného `lag_s`:

- hvězda s **konvektivní obálkou** disipuje turbulentní viskozitou (rovnovážný slap; Zahn 1977;
  Rasio et al. 1996),
  $$ \left(\frac{k}{T}\right)_\mathrm{c} = \frac{2}{21}\,\frac{f_\mathrm{conv}}{\tau_\mathrm{conv}}\,\frac{M_\mathrm{env}}{M}, \qquad f_\mathrm{conv} = \min\left[1, \left(\frac{P_\mathrm{tid}}{2\tau_\mathrm{conv}}\right)^2\right], $$
  s hmotností obálky $M_\mathrm{env}$, obrátkovým časem vírů
  $\tau_\mathrm{conv} = 0{,}4311\,[M_\mathrm{env} R_\mathrm{env}(R - R_\mathrm{env}/2)/(3L)]^{1/3}$ yr
  ($R_\mathrm{env}$ je hloubka obálky, $L$ zářivý výkon, ve slunečních jednotkách) a periodou
  slapového buzení $P_\mathrm{tid}$, $1/P_\mathrm{tid} = |1/P_\mathrm{orb} - 1/P_\mathrm{spin}|$
  (Hurley et al. 2002, rov. 33; faktor 2 v $f_\mathrm{conv}$, jejich rov. 32, je tam proto, že slap
  během této periody proběhne dvakrát). Faktor $f_\mathrm{conv}$ snižuje viskozitu, když je slapová
  perioda kratší než obrátkový čas – standardně kvadraticky (Goldreich & Nicholson 1977; Hurley et
  al. 2002), nebo lineárně (Zahn 1966, 1989; Zahn & Bouchet 1989) s
  `physics.tidal_convective_reduction: linear`;
- hvězda se **zářivou obálkou** disipuje dynamickým slapem (Zahn 1975, 1977),
  $(k/T)_\mathrm{r} = 1{,}9782\times10^4\,(M R^2/a^5)(1+q)^{5/6} E_2$ yr⁻¹ s
  $E_2 = 1{,}592\times10^{-9} M^{2{,}84}$ ($q$ je poměr hmotností, ve slunečních jednotkách).

Veličina $(k/T)$ podle Huta (1981) odpovídá konstantnímu časovému zpoždění přes
$T = R^3/(G M \tau_\mathrm{Hut})$ a $\tau_\mathrm{Hut} = 2\,\Delta t$ (SWINGERsim používá Loveovo
číslo $2k_2$ podle Leconta et al. 2010; faktor byl ověřen numericky). Zpoždění se počítá
z počáteční dráhy a rotace a u vyvíjejících se hvězd znovu v každém úseku; `lag_s` zadané pro
hvězdu se zachová a `physics.tidal_efficiency` disipaci násobí.

**Kalibrace.** Tento předpis je široce používán (BSE, MSE, COSMIC), jeho kalibrace je však
nejslabším článkem dlouhodobého vývoje. Na samotné hlavní posloupnosti se dvojhvězdy slunečního typu
za 4 Gyr zakruhovají jen do periody asi 3,8 d, kdežto hvězdokupa M67 ukazuje 12,1 +1,0/−1,5 d (Meibom
& Mathieu 2005, tab. 3); jen pro hlavní posloupnost by byl potřeba `tidal_efficiency` asi 150 (100 dává
11,2 d, 200 dává 12,9 d).

### Inerciální vlny {#inertial-waves}

V rotující konvektivní obálce budí slapové působení s frekvencí pod $2\Omega$ *inerciální vlny*,
jejichž vratnou silou je Coriolisova síla. Jejich disipace – dynamický slap konvektivních zón – je
mnohem silnější než tlumení rovnovážného slapu turbulentní viskozitou, které je slabé, protože
slapová perioda je krátká ve srovnání s konvektivním obrátkovým časem (Ogilvie 2013; Mathis 2015;
Barker 2020, 2022). Po středování přes slapovou frekvenci závisí jen na stavbě a rotaci,

$$ \frac{1}{\langle Q'_\mathrm{IW}\rangle} = \epsilon^2\,C, \qquad \epsilon^2 = \frac{\Omega^2 R^3}{GM}, $$

kde $Q'$ je modifikovaný slapový faktor kvality, $\Omega$ úhlová rychlost rotace a $C$ strukturní
koeficient. S `physics.inertial_waves: true` (vyžaduje `tidal_lag: hurley2002`) počítá SWINGERsim $C$
pro hvězdu se zářivým jádrem o poloměru $\alpha R$ a hmotnosti $\beta M$ (z modelu obálky)
dvouvrstvým vzorcem Ogilvieho (2013) a pro plně konvektivní hvězdu (pod 0,35 M☉) pro polytropu
$n = 1{,}5$, $C = 1/130{,}9$, spočtené jednou metodou Barkera (2020, rov. 12–38; 2022, příl. A) pro
spojité profily. Tlumení excentricity $d\ln e/dt = -(225/16)(M_2/M)(R/a)^5\,n/Q'$ (Barker 2022,
rov. 1) se převede na ekvivalentní konstantní časové zpoždění tak, aby souhlasilo s rychlostí
tlumení modelu zpoždění pro synchronní hvězdu, $-(21/2)(k/T)\,q(1+q)(R/a)^8$ (Hut 1981), a přičte se
k rovnovážnému slapu. Inerciální vlny se vypnou, když hvězda rotuje pomaleji než polovinou oběžné
frekvence (žádná budicí frekvence pod $2\Omega$).

Numerická metoda reprodukuje dvouvrstvý vzorec ve tvaru Barkera (2020, rov. 39) pro vyhlazený
dvouvrstvý profil na 0,2 %, což je kontrola řešiče a jeho normování. Pro hvězdy hlavní posloupnosti
0,3 a 1 M☉ rotující s periodou 10 d dává $Q' = 0{,}9$ a $1{,}0\times10^7$, jak nachází Barker (2022)
s profily MESA. Pro realistické profily hvězd slunečního typu dvouvrstvý model podhodnocuje disipaci
asi dvakrát a pro hvězdy typu F nad 1,2 M☉ o řád i více (Barker 2020).

### Slapy před hlavní posloupností {#pms-tides}

Mladé hvězdy jsou velké a plně konvektivní a většina cirkularizace dvojhvězd slunečního typu
proběhne dřív, než dosáhnou hlavní posloupnosti (Zahn & Bouchet 1989). S hvězdami MIST ve fázi před
hlavní posloupností ([Hvězdy](#mist-model)) působí rovnovážný slap i inerciální vlny i tam. Dvojice
stejných hvězd 1 M☉ spuštěné v 1 Myr jsou v 50 Myr kruhové až do počáteční periody 7,9 d (Zahn
& Bouchet: 7,2–8,5 d). V periodě, kterou používají pozorovatelé – výsledné periodě právě
zakruhovaných dvojhvězd –, jsou kruhové do 6,6 d ve 3 Myr (dvojhvězdy před hlavní posloupností:
7,1 ± 1,2 d; Meibom & Mathieu 2005, tab. 3) a do 7,5 d ve 4 Gyr (M67: 12,1 +1,0/−1,5 d).
S inerciálními vlnami zakruhová samotná hlavní posloupnost dvojice hvězd slunečního typu do 7,5 d
za 4 Gyr (3,8 d jen s rovnovážným slapem). Fáze před hlavní posloupností je tedy reprodukována;
další růst na hlavní posloupnosti nikoli. Barker (2022) s disipací inerciálních vln středovanou
přes frekvenci, spočtenou z profilů MESA, a se synchronizovanými rotacemi nachází periody
cirkularizace, které s věkem na hlavní posloupnosti rostou, jak se pozoruje. Dvouvrstvý odhad
SWINGERsimu se stavbou SSE a MIST a vyvíjejícími se rotacemi to nereprodukuje; rozdíl spočívá ve
stavbě hvězd (dvouvrstvý model podhodnocuje disipaci realistických profilů slunečního typu, viz
výše) a v historii rotace a zůstává otevřenou otázkou. Rezonanční zachycení gravitačních módů, které
SWINGERsim nemodeluje, působí hlavně před hlavní posloupností: Zanazzi & Wu (2021) z něj dostávají
kruhové periody asi 4–6 d a zanedbatelný příspěvek hlavní posloupnosti. Výsledek závisí i na
začátku výpočtu: mladší, větší hvězda se zakruhuje více (při startu v 0,15 Myr 11,6 d již ve 3 Myr).
Před hlavní posloupností převládají inerciální vlny, takže lineární a kvadratické snížení
turbulentní viskozity se tam shodují do 1 %.

## Meze {#tides-limits}

- Model rovnovážných slapů předpokládá malé, kvazistatické deformace. Dynamické slapy hvězd se
  zářivou obálkou vstupují jen přes Zahnův předpis ve zpoždění; rezonanční zachycení a nelineární
  slapy se nemodelují.
- V sekulárních backendech zůstávají rotační osy kolmé k drahám; odklon rotace od dráhy a precesi
  rotační osy sleduje jen N-částicový backend.
- Rychlosti disipace nesou výše diskutovanou nejistotu kalibrace; periody cirkularizace starých
  hvězdokup nereprodukuje žádný model rovnovážných slapů, včetně tohoto.
- Bez `physics.tides` je růst excentricity v cyklech ZLK těsných dvojhvězd nadhodnocen
  (Vokrouhlický 2016).

## Implementace {#implementation-tides}

**Externí knihovna.** Slapové a rotační síly a vývoj rotace v N-částicovém backendu zajišťuje efekt
REBOUNDx `tides_spin` (Lu et al. 2023; Tamayo et al. 2020). **Vlastní kód.** Převod konvence
časového zpoždění ($\tau = \Delta t/2$, zjištěný rozborem zdrojového kódu REBOUNDx a ověřený proti
Lecontovi et al. 2010); počáteční rotace (synchronní, pseudosynchronní podle Huta 1981, zadané
periody a sklony rotačních os); Sterneovy (1939) rychlosti stáčení apsid a jejich rozklad na
příspěvky; rovnice konstantního časového zpoždění podle Leconta et al. (2010) v sekulárních
backendech a slapové časové škály; zpoždění ze stavby hvězdy podle Hurleyho, Touta & Polse (2002)
s oběma sníženími viskozity; rotační oprava $k_2$; a disipace inerciálních vln – dvouvrstvý vzorec
Ogilvieho (2013) a vlastní řešič metody spojitých profilů podle Barkera (2020, 2022) pro polytropy
a profily MESA.
