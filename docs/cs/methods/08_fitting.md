# Prokládání pozorovaných dat {#fitting}

Předchozí kapitoly spouštějí model dopředu: od drah a hvězd ke světelným křivkám, rychlostem
a okamžikům minim. Tato kapitola popisuje obrácený krok – přizpůsobení modelu naměřeným datům –, který
SWINGERsim provádí v posloupnosti kroků nazývané *řetězec proložení*. Řetězec začíná analytickým
popisem každého zákrytového páru (krok 1), volitelně ukotví publikované elementy k naměřeným minimům
a končí jedním N-částicovým modelem celé soustavy proloženým všemi daty najednou (krok 2), s robustním
zacházením s chybami, volitelnými nejistotami z metody Monte Carlo s Markovovými řetězci (MCMC)
a předpověďmi budoucích minim.

![**Obrázek 19.** Řetězec proložení (nový obrázek pro tento svazek). Rámečky jsou kroky, šipky tok
informací. Horní řada: data (okamžiky minim podle zdroje, radiální rychlosti, úseky světelných křivek,
sklony nebo meze zákrytů v daných epochách) vstupují do kroku 1, analytického modelu pro každý pár,
jehož vnější dráha může spolu s publikovanými elementy a parametry hvězd poskytnout výchozí
konfiguraci. Vpravo: volitelné ukotvení periody a fáze každého zákrytového páru. Uprostřed: krok 2,
proložení jednoho N-částicového modelu všemi daty metodou nejmenších čtverců, se střídá
s „opakovanými proloženími“ (vlevo), která přeškálují chyby pro každý zdroj, vyřadí odlehlé body
a odhadnou faktor červeného šumu, dokud se nic nemění. Dole: nejistoty z kovarianční matice a volitelně
z emcee a výsledky – proložená konfigurace použitelná všemi ostatními příkazy a předpovězená minima
s jejich pásy. Modré rámečky jsou analytické nebo lineární kroky; oranžové rámečky volají N-částicový
model mnohokrát.](../../figures/methods/fit_chain.png)

Obrázek 19 je mapou této kapitoly. Každý rámeček odpovídá jednomu z následujících oddílů a smyčka mezi
oběma oranžovými rámečky je robustní částí proložení: model chyb není pevně dán předem, ale odhaduje se
z reziduí předchozího proložení.

## Data {#fit-data}

Proložení používá až čtyři druhy dat, všechny na jedné časové ose (BJD v systému TDB; tabulky v HJD
nebo JD v UTC se převádějí pomocí nebeských souřadnic objektu):

- **Okamžiky minim** (tabulka `fit.minima`): pár, typ (primární nebo sekundární), čas, chyba, zdroj
  (přístroj, pozorovatel nebo článek) a metoda. Na *zdroji* záleží: chyby a odlehlé body se zpracovávají
  pro každý zdroj zvlášť.
- **Radiální rychlosti** (`fit.rv`): hvězda, čas, rychlost, chyba, zdroj; modelová rychlost se počítá
  v okamžiku vyslání a volitelně jako rychlost středu světla
  ([Pozorovatelné veličiny](#radial-velocities)).
- **Světelné křivky** (`fit.lc`): úseky normovaného toku s chybami, například okna kolem zákrytů nebo
  binované sektory, v daném fotometrickém pásmu a s danou dobou expozice.
- **Sklony v závislosti na čase** (`fit.inclination`): sklony zákrytového páru naměřené v několika
  epochách, nebo jednostranné podmínky ([Data o sklonu](#inclination-data)).

*Pár* je obvykle dvojhvězda ze dvou hvězd, může to však být libovolná dráha: konjunkce širší dráhy –
těžiště jejích dvou složek seřazená podél zorného paprsku – jsou středy zákrytů třetím tělesem, takže
zákryt třetí hvězdou lze prokládat stejně jako minimum dvojhvězdy.

## Krok 1: analytický model O−C pro každý pár {#stage-1}

Krok 1 je analytický model O−C, který skupina SQUADRA prokládá programem OCFit (Gajdoš & Parimucha
2019), po jednom páru, nebo několik párů s jednou vnější dráhou:

$$T(n, j) = T_0 + P\,(n + j/2) + Q\,n^2 + \Delta_\mathrm{LTTE}(t) + D_j(\omega_n),$$

s číslem cyklu $n$, typem minima $j$ (0 primární, 1 sekundární), referenční epochou $T_0$, periodou
$P$, kvadratickým koeficientem $Q$ (konstantní změna periody, $Q = \tfrac12 P\dot P$), členem
světelného času vnější dráhy ([Pozorovatelné veličiny](#ltte)) a apsidálním členem $D_j$. Členy se
volí pro každý pár (`fit.stage1.terms`, `pairs`).

**Apsidální člen.** Apsidální stáčení excentrického páru, $\omega_n = \omega_0 + \dot\omega\,n$,
vstupuje přes okamžik nejmenší promítnuté vzdálenosti. Při $i = 90°$ je to konjunkce
$\nu + \omega = \pi/2$ (primární) nebo $3\pi/2$ (sekundární) a

$$ D_j = \frac{P}{2\pi}\left[M(\nu_j) + \omega - \left(\frac{\pi}{2} + j\pi\right)\right], $$

je doba od střední ke skutečné konjunkci ($M(\nu)$ je střední anomálie příslušná k pravé anomálii
$\nu$). To je přesný tvar řady Giméneze & Bastera (1995), kterou používá OCFit; pro $i < 90°$ se minimum
promítnuté vzdálenosti hledá numericky (členy $e\cot^2 i$ podle Giméneze & Garcii-Pelaya 1983).
Apsidální perioda je $U = 2\pi P/\dot\omega$.

**Společná vnější dráha.** S `shared_outer` sdílejí páry jednu vnější dráhu (společné $P_3$, okamžik
průchodu periastrem $T_3$ a $e_3$, periastrum druhého páru je opačné k prvnímu) s oddělenými hodnotami
$a\sin i$; poměr hodnot $a\sin i$ dává poměr hmotností párů, $m_A/m_B = (a\sin i)_B/(a\sin i)_A$.

**Metoda.** Vážená metoda nejmenších čtverců startovaná z mřížky kruhových proložení LTTE přes vnější
periodu (a pro apsidální člen přes protifázi primárních a sekundárních minim); časy se prokládají vůči
referenční epoše blízko středu dat, aby si derivace zachovaly přesnost. Chyby pocházejí z kovarianční
matice, škálované faktorem $\sqrt{\chi^2_\mathrm{red}}$, když redukované $\chi^2$ převyšuje 1, nebo
z běhu MCMC (níže). Odlehlé body lze vyřadit a chyby přeškálovat pro každý zdroj jako v kroku 2. Model
souhlasí s řadou programu OCFit na $10^{-3}$ s při $e = 0{,}3$ a proložení vlastních N-částicových
okamžiků minim soustavy 2+2 spočtených SWINGERsimem vrátí keplerovské $a\sin i$ obou párů na 0,3 %
a vnější periodu na $4\times10^{-5}$ (viz *Testování, část I*).

**Meze.** Krok 1 pojímá každý pár jako keplerovskou dvojhvězdu na keplerovské vnější dráze. Nedokáže
popsat dynamické poruchy, precesi ani vazbu párů v kompaktní soustavě; jeho vnější dráha a $a\sin i$
jsou výchozím bodem pro krok 2.

## Ukotvení k současnosti {#anchoring}

Publikovaná fotodynamická řešení udávají oskulační elementy v dané epoše. Když z nich N-částicový
model SWINGERsimu vystartuje, má střední vnitřní periodu mírně odlišnou od periody naměřené
v minimech (elementy jsou oskulační, perioda viděná v minimech je střední; druhý kód se může lišit
fyzikou) a spočtená minima se v průběhu let od naměřených vzdalují. *Ukotvení* (`fit.stage2.anchor`)
tento drift odstraní před jakýmkoli srovnáním: perioda $P$ a střední anomálie $M_0$ každého
zákrytového páru se upraví tak, aby se střední lineární efemerida modelových zákrytů rovnala efemeridě
naměřených minim. Hrubý krok přeškáluje $P$ faktorem $P_\mathrm{data}/P_\mathrm{model}$ a posune $M_0$
o zpoždění modelového zákrytu nejbližšího $t = 0$; jemný krok proloží přímku rozdíly naměřený −
modelový okamžik v epochách minim každého páru (bez zkreslení nerovnoměrným vzorkováním zakřiveného
O−C) a opakuje se, dokud se žádné minimum neposune o více než 0,1 min. Vše ostatní zůstává podle
konfigurace. Krok 2 odstraňuje tentýž rozdíl ze své podstaty, protože prokládá oskulační elementy,
které minima reprodukují; ukotvení poskytuje dobrý výchozí bod a samo o sobě je užitečné, když je
potřeba jen předpověď.

## Krok 2: společný N-částicový model {#stage-2}

Krok 2 (`fit.stage2`) prokládá jeden N-částicový model celé soustavy všemi daty najednou. Model je
stejný jako při běžném běhu (newtonovská gravitace, relativita 1PN, slapy, předepsané změny period
podle konfigurace), takže efekt světelného času, apsidální stáčení, dynamické poruchy i precese plynou
z jednoho modelu a proloženou konfiguraci lze přímo spustit (stabilita, vývoj dopředu a zpět, světelné
křivky, předpovědi).

**Modelové okamžiky.** Pro každé naměřené minimum se soustava integruje do okamžiku krátce před
očekávanou konjunkcí; dvoučásticový pohyb páru z tohoto stavu dá nejbližší konjunkci správného typu
(nejmenší promítnutá vzdálenost; primární, když je první složka ta vzdálenější); N-částicová integrace
pokračuje až k ní a zpřesní ji. Pozorovaný okamžik přičítá dobu letu světla od těžiště páru,
$t_\mathrm{obs} = t + z_\mathrm{cm}/c$. Minima před $t = 0$ pocházejí z druhé integrace zpět v čase.
Všechny požadavky všech druhů dat (minima, rychlosti v jejich okamžicích vyslání, sklony v jejich
epochách) jdou do jedné fronty seřazené podle času, takže jedno vyhodnocení modelu stojí jednu
integraci přes rozsah dat v každém směru. Okamžiky souhlasí s okamžiky syntetických minim při 400
vzorcích na periodu na 0,06 s.

### Volné parametry {#free-parameters}

Parametry, které lze prokládat, uvádí tabulka 13.

**Tabulka 13.** Parametry, které může krok 2 prokládat (`fit.stage2.free`). *Název* je tvar
v konfiguraci; *Význam* veličina; *Výchozí meze* oblast, v níž optimalizátor hledá, pokud `bounds`
neudává jinou. Vše, co není uvedeno jako volné, zůstává podle konfigurace.

| Název | Význam | Výchozí meze |
|---|---|---|
| `orbits.<o>.P` | perioda [d] | 0,5–2 × počáteční |
| `orbits.<o>.M0`, `omega`, `Omega` | střední anomálie, argument periastra, uzel [deg] | žádné |
| `orbits.<o>.e` | excentricita | 0–0,95 |
| `orbits.<o>.ecosw`, `esinw` | $e\cos\omega$, $e\sin\omega$ (nahrazují $e$, $\omega$ pro malé $e$) | ±0,95, $e \le 0{,}95$ |
| `orbits.<o>.i` | sklon [deg] | 0–180 |
| `orbits.<o>.Pdot` | předepsaná změna periody [d yr⁻¹] | žádné |
| `masses.<m>` | hmotnost [M☉] | 0,05–20 × počáteční |
| `stars.<m>.R`, `Teff` | poloměr [R☉], efektivní teplota [K] | 0,2–5 ×, 0,5–2 × počáteční |
| `stars.<m>.apsidal_k2`, `rg2` | strukturní konstanty | $10^{-5}$–0,74, 0,005–0,4 |
| `rv.gamma`, `rv.offset.<source>` | systémová rychlost, nulový bod zdroje [km s⁻¹] | žádné |

Prokládání $k_2$ (nebo poloměrů) mění apsidální stáčení v podmínku pro stavbu hvězd; prokládání `Pdot`
umožňuje proložit parabolu v O−C společně s LTTE ([Přenos hmoty a vývoj](#semidetached-contact)).
Otočení celé soustavy kolem zorného paprsku nezmění žádná data, takže při volném uzlu každé dráhy jsou
měřeny jen jejich rozdíly; SWINGERsim pak drží uzel nejvnější dráhy na jeho nakonfigurované hodnotě
a oznámí to.

**Co data omezují.** Samotná minima určují periody, fáze, excentricity a argumenty periastra všech drah
a prostřednictvím amplitudy světelného času a dynamických poruch i poměry hmotností. Sklony
zákrytových párů a absolutní hmotnosti potřebují světelné křivky a radiální rychlosti; bez nich je
třeba je zafixovat na hodnotách z řešení světelné křivky nebo spektroskopického řešení, nebo je omezit
apriorními rozděleními. Vzájemné sklony a uzly jsou omezeny jen v kompaktních soustavách s měřitelným
dynamickým (nikoli světelně-časovým) signálem.

### Věrohodnost {#likelihood}

Proložení minimalizuje

$$ \chi^2 = \sum_k \left(\frac{d_k - m_k(\mathbf{p})}{\sigma_k}\right)^2 + \sum_j \left(\frac{p_j - \mu_j}{s_j}\right)^2 , $$

kde $d_k$ jsou data (okamžiky, rychlosti, toky, sklony), $m_k(\mathbf{p})$ model při parametrech
$\mathbf{p}$, $\sigma_k$ chyby po níže popsaných opravách a druhý součet gaussovská apriorní rozdělení
parametrů $p_j$ se středními hodnotami $\mu_j$ a šířkami $s_j$. Až na konstantu je to $-2\ln$ gaussovské
věrohodnosti násobené gaussovskými apriorními rozděleními. Úhly v apriorních rozděleních se převádějí
do rozsahu ±180°. Datový bod, pro nějž model neexistuje (žádná konjunkce daného typu, hvězda za koncem
své vývojové dráhy), se počítá jako $10^3\sigma$, takže se optimalizátor od takových oblastí vzdálí.

**Světelné křivky** vstupují s kalibrací pro každý úsek, která se řeší lineárně při každém vyhodnocení
a není součástí volných parametrů:

$$ F_\mathrm{obs}(t) = a\,F_\mathrm{model}(t) + b + c\,(t - t_\mathrm{mid}) + d\,(t - t_\mathrm{mid})^2 , $$

kde $a$ škáluje hloubky zákrytů (zředění třetím světlem nebo jiné normování), $b$ je nulový bod a $c$,
$d$ popisují pomalé trendy (stupeň polynomu 0–2); se `scale: false` musejí hloubky plynout z hvězd tak,
jak jsou zadány. Počet stupňů volnosti proložení se o tyto koeficienty sníží. Modelová světelná křivka
používá dobu expozice dat. Tvary hvězd se volí pro každou světelnou křivku: se `shapes: auto` dostane
hvězda kruhové dvojhvězdy ze dvou hvězd svou Rocheovu plochu, když její objemový poloměr přesáhne podíl
`roche_threshold` (výchozí 0,5) jejího Rocheova laloku – koule se tam mýlí asi o 2 ppt
([Pozorovatelné veličiny](#distorted-stars)) –, jinak zůstane koulí.

**Deformované hvězdy v proložení.** Sítě jsou výpočetně drahé. S `fidelity: mapped` (výchozí) běží
proložení na hrubé síti (24 × 48 směrů) plus rozdíl *přesná − hrubá* vyhodnocený v aktuálním řešení
s přesnou sítí (64 × 128); rozdíl se přepočítává v každém novém řešení, dokud se nezmění o méně než
0,3σ (nejvýše pět průchodů). Toto *mapování výstupního prostoru* konverguje k řešení přesného modelu
s jedním přesným vyhodnocením na průchod; $\chi^2$ a světelné křivky uváděné v řešení jsou ty z přesné
sítě.

### Data o sklonu {#inclination-data}

Tabulka sklonů zákrytového páru naměřených v několika epochách (řešení světelných křivek oddělených
souborů dat, jak je pro V907 Sco publikovali Zasche et al. 2023a) se porovnává se sklonem dráhy páru
v N-částicovém modelu v týchž okamžicích. Zákryty měří jen $|\cos i|$, takže modelový sklon se sklápí
do 0°–90°, stejně jako naměřené hodnoty nad 90° ([Hierarchie a dráhy](#orbital-elements)). Řádek může
být i jednostranný:

- `upper` / `lower`: $i_\mathrm{model} \le i$ nebo $\ge i$ (například „zhruba 75° nebo méně“ ze
  světelné křivky bez zákrytů);
- `no_eclipse`: pár v dané epoše nezakrýval, tj. $i_\mathrm{model}$ leží pod mezí zákrytů modelu
  $\arccos[(R_1+R_2)/r_\mathrm{conj}]$, s modelovými poloměry a vzdáleností $r_\mathrm{conj}$ v bližší
  konjunkci;
- `eclipse`: pár zakrýval, $i_\mathrm{model}$ je nad touto mezí.

Jednostranný řádek přispívá hodnotou $(\mathrm{bound} - i_\mathrm{model})/\sigma$ jen tehdy, když jej
model porušuje (*kloub*, hinge), a jinak ničím; $\sigma$ určuje, jak ostře. Splněné jednostranné řádky
se nezapočítávají jako data do počtu stupňů volnosti. Dolní mez chyb `min_sigma_deg` vyjadřuje
systematické chyby řešení světelných křivek.

### Optimalizátor a kovariance {#optimiser}

$\chi^2$ se minimalizuje reflektivní metodou nejmenších čtverců s oblastí důvěry a s mezemi (SciPy
`least_squares`), volitelně s robustní ztrátovou funkcí (`soft_l1`, `huber`, `cauchy`, které snižují
vliv reziduí za hranicí `f_scale` σ). Jacobiho matice pochází z **centrálních diferencí**, dvou
N-částicových integrací na volný parametr a iteraci, počítaných paralelně; dopředné diference jsou
zkresleny křivostí reziduí (například kumulativním účinkem změny periody). Na mezi nebo tam, kde model
na jedné straně selhává (hvězda za koncem své vývojové dráhy, nestabilní dráha), se použije
jednostranná diference a parametr se označí jako ležící na okraji platnosti modelu, kde jeho chyba
neplatí. Optimalizátor se zastaví, když se $\chi^2$ za iteraci změní o méně než $10^{-5}$ relativně
nebo 0,05 absolutně (hluboko pod změnou 1, která posune parametr o jedno σ), nebo po `max_iterations`
(výchozí 30).

Kovarianční matice se počítá ze singulárního rozkladu Jacobiho matice se škálovanými sloupci, nikoli
inverzí $J^\mathsf{T}J$, která číslo podmíněnosti umocní na druhou: s parametry velmi rozdílných škál
(perioda známá na $10^{-5}$ d, úhly na 0,1°) a korelacemi blízkými 0,999 pseudoinverze $J^\mathsf{T}J$
tiše vypustí nejslabší směr a dá chyby příliš malé. Směry se singulární hodnotou pod $10^{-12}$
největší se označí jako neomezené. Chyby datových řádků se škálují faktorem
$\sqrt{\max(\chi^2_\mathrm{red}, 1)}$; apriorní rozdělení si ponechávají své šířky. Silně korelované
parametry ($|r| > 0{,}95$) se uvádějí ve zprávě – typicky $M_0$, $\omega$ a $P$ vnější dráhy pokryté
méně než asi dvěma periodami –, protože jejich lineární chyby platí jen tak dlouho, dokud je údolí
$\chi^2$ přímé.

Test uzavření na TIC 120362137 (syntetická minima ve 147 epochách skutečných minim z TESS,
12 porušených volných parametrů) konverguje z $\chi^2 = 1{,}9\times10^8$ na 107
($\chi^2_\mathrm{red} = 0{,}79$) za 331 vyhodnocení modelu a vrátí všech 12 parametrů v rámci 1,5σ
lineárních chyb (viz *Testování, část I*).

## Chyby podle zdroje, jitter a odlehlé body {#errors-outliers}

Okamžiky minim z různých zdrojů mají chyby různé kvality. Uváděné chyby starých fotografických,
vizuálních nebo přehlídkových minim bývají hluboko pod jejich skutečným rozptylem a i formální chyby
minim z TESS jsou pod tím, co kterýkoli model reprodukuje (redukované $\chi^2$ 86–1432 v testu na
datech SQUADRA). Ponecháno, jak je, několik set starých minim s optimistickými chybami převáží ta
přesná. SWINGERsim proto odhaduje model chyb z reziduí v **opakovaných proloženích**: proložit, opravit
chyby a vyřadit odlehlé body, proložit znovu z řešení, dokud se ani vyřazené body, ani opravy nezmění
o více než 2 % (nejvýše `refits` proložení, výchozí 4). Minima, rychlosti (skupiny podle zdroje)
a světelné křivky (skupiny podle souboru a úseku) se zpracovávají stejně.

**Oprava chyb podle zdroje** (`error_scale: per_source`; `global` používá jednu opravu pro všechna
minima). Zdroje s méně než `min_source_minima` (výchozí 5) body sdílejí jednu skupinu. Každá skupina $s$
dostane buď

- **faktor** (`error_model: scale`), $\sigma_k \to f_s\,\sigma_k$ s $f_s = \sqrt{\chi^2_s/\nu_s} \ge 1$,
  nebo
- **jitter** přičtený kvadraticky (`error_model: jitter`), $\sigma_k \to \sqrt{\sigma_k^2 + j_s^2}$,
  kde $j_s \ge 0$ je kořen rovnice $\sum_{k \in s} r_k^2/(\sigma_k^2 + j_s^2) = \nu_s$,

kde $\chi^2_s$ je příspěvek skupiny, $r_k$ její rezidua a $\nu_s$ její počet stupňů volnosti: volné
parametry se mezi skupiny rozdělují úměrně počtu jejich bodů (úsek světelné křivky navíc odečítá své
kalibrační koeficienty). Faktor zachovává relativní váhy uvnitř zdroje; jitter působí jako dolní mez
chyb – systematické chyby modelu a časování se k chybě měření přičítají kvadraticky –, což se hodí pro
zdroje velmi přesných minim. Pro každý zdroj se chyby nikdy nesníží pod uváděné.

**Vyřazování odlehlých bodů** (`clip_sigma`, např. 5). Bod se vynechá, když jeho normované reziduum
překročí `clip_sigma` robustních směrodatných odchylek jeho skupiny, $1{,}4826\,\mathrm{median}|r_k/\sigma_k|$
(alespoň 1, aby zdroj, který se prokládá lépe než jeho chyby, neztratil žádné dobré body). Nejhorší
body jdou první, nejvýše `max_clip_fraction` (výchozí 10 %) všech; dříve vynechaný bod se vrátí, když
jej současný model proloží v rámci 0,8 × prahu. Vynechané body zůstávají ve výstupech.

**Červený šum** (`fit.lc.red_noise`). Rezidua světelných křivek bývají v čase korelovaná (přístrojové
trendy, oscilace, skvrny) a chyby bílého šumu pak podhodnocují chyby parametrů. *Faktor červeného šumu*
$\beta$ metody časového středování (Pont, Zucker & Queloz 2006; Winn et al. 2008) porovnává rozptyl
reziduí středovaných v binech po $N$ bodech s očekáváním pro bílý šum,

$$ \beta_N = \frac{\sigma_N}{\sigma_1 N^{-1/2}\sqrt{M/(M-1)}}, $$

kde $\sigma_1$ je rozptyl nebinovaných reziduí, $\sigma_N$ rozptyl středních hodnot binů a $M$ počet
binů; $\beta$ je medián přes šířky binů 20, 40, 60 a 120 min, nejméně 1. Chyby každého úseku se
v opakovaných proloženích násobí $\beta$.

## Nejistoty z MCMC {#mcmc}

Lineární chyby spoléhají na lokálně lineární model. V dlouhém zakřiveném údolí $\chi^2$ – typickém pro
vnější dráhu pokrytou méně než asi dvěma periodami – vycházejí příliš malé. Se `sampler: emcee`
vzorkuje SWINGERsim aposteriorní rozdělení kolem řešení metodou nejmenších čtverců afinně invariantním
ansámblovým vzorkovačem Goodmana & Weareho (2010) v implementaci emcee (Foreman-Mackey et al. 2013).

- **Aposteriorní rozdělení.**
  $\ln p = -\tfrac12[\chi^2_\mathrm{data}/\max(\chi^2_\mathrm{red}, 1) + \chi^2_\mathrm{prior}]$
  uvnitř mezí, takže chyby dat jsou zvětšeny stejně jako u lineárních chyb, zatímco apriorní
  rozdělení si ponechávají své. Pro dráhy prokládané pomocí $e\cos\omega$ a $e\sin\omega$ se přidá
  jakobián $1/e$: rovnoměrné apriorní rozdělení v $(e\cos\omega, e\sin\omega)$ je apriorním
  rozdělením $\propto e$, které téměř kruhové dráhy odtlačuje od $e = 0$ (Lucy & Sweeney 1971; Ford
  2006; Eastman, Gaudi & Agol 2013), a oprava dělá apriorní rozdělení rovnoměrným v $e$ a $\omega$.
- **Walkerů** – paralelních řetězců ansámblového vzorkovače – je nejméně $2n + 2$ pro $n$ volných
  parametrů a sudý počet (emcee posouvá vždy polovinu z nich). Počet nezávisí na počítači, takže stejné
  semínko dá stejný řetězec při libovolném počtu jader. Startují v kouli o velikosti poloviny kovarianční matice kolem řešení, tj. rozložení podél
  údolí.
- **Kroky.** Každý walker udělá `mcmc_steps` kroků (výchozí 300); první polovina se zahodí jako
  zahřívací fáze (burn-in). Chyby jsou polovinou rozsahu 16–84 % vzorků.
- **Autokorelační čas.** Integrovaný autokorelační čas $\tau$ je počet kroků, po nichž řetězec
  zapomene svůj stav; řetězec o $N$ krocích obsahuje asi $N/\tau$ nezávislých vzorků na walkera.
  SWINGERsim uvádí $\tau$ každého parametru, jeho maximum a poměr kroky/$\tau$ a varuje, když je řetězec
  po zahřívací fázi kratší než 10 $\tau$.

**Praktické pravidlo konvergence.** Řetězec dlouhý alespoň 50 $\tau$ je doporučením autorů emcee
(Foreman-Mackey et al. 2013) pro spolehlivý odhad $\tau$ a kvantilů aposteriorního rozdělení.
S `mcmc_tau_target: 50` se řetězec prodlužuje vždy o čtvrtinu své délky, dokud nedosáhne 50 $\tau$
nebo `mcmc_max_steps`. Každé vyhodnocení modelu je N-částicová integrace, takže 50 $\tau$ může být
drahé. **Je-li řetězec krátký**, je třeba to tak uvést: udat počet walkerů, kroků, zahozenou zahřívací
fázi, $\tau_\mathrm{max}$, kroky/$\tau$ a podíl přijatých kroků; intervaly z MCMC uvádět jako orientační
(bývají příliš úzké, protože walkeři dosud neprozkoumali okraje rozdělení) a vedle nich uvést lineární
chyby a seznam silně korelovaných parametrů.

## Předpovědi {#predictions}

S `fit.output.predict` předpovídá proložený model zákryty zvolených párů na daný počet let dopředu
(od současnosti, nebo od posledního naměřeného minima, je-li pozdější). Každé předpovězené minimum nese
pás nejistoty: sady parametrů se losují ze vzorků emcee (pokud proložení emcee použilo) nebo
z kovarianční matice, každá sada se jednou integruje (paralelně) a pás je rozsah 16–84 % okamžiků
téhož zákrytu napříč sadami. Uvádějí se jen zákryty, které v proloženém modelu nastanou (zákrytový
parametr pod 1). Pásy z kovarianční matice dědí její lineární aproximaci; v zakřivených údolích dávají
vzorky MCMC poctivější pás.

## Meze {#fitting-limits}

- Proložení najde lokální minimum blízko výchozího bodu; špatný start (chybná vnější perioda) může
  skončit ve špatném minimu. Start poskytuje krok 1, ukotvení a publikovaná řešení.
- Lineární chyby jsou spolehlivé jen v přímých údolích $\chi^2$; chyby z MCMC jen pro zkonvergované
  řetězce.
- Systematické chyby modelu (chybějící fyzika, pevné parametry hvězd) vstupují do oprav chyb podle
  zdroje, nikoli do parametrů; jitter velký ve srovnání s uváděnými chybami je známkou chybějící
  fyziky.
- Model světelných křivek dědí meze z kapitoly [Pozorovatelné veličiny](#observables): žádné skvrny,
  Rocheovy tvary jen pro kruhové synchronizované dvojhvězdy.

## Implementace {#implementation-fitting}

**Externí knihovny.** `least_squares` ze SciPy (reflektivní metoda s oblastí důvěry) jako
optimalizátor; emcee (Foreman-Mackey et al. 2013) jako vzorkovač MCMC a pro autokorelační časy;
REBOUND/REBOUNDx v každém vyhodnocení modelu; Astropy pro převody časových systémů. **Vlastní kód.**
Model kroku 1 (LTTE, přesný apsidální člen, společná vnější dráha, startovací mřížky); ukotvení; model
kroku 2 s frontou událostí, která obsluhuje minima, rychlosti a sklony z jedné integrace v každém
směru; zobrazení parametrů, meze, zacházení s $e\cos\omega$/$e\sin\omega$ a pevný uzel nejvnější
dráhy; Jacobiho matice z centrálních diferencí s jednostrannými kroky na mezích a na okraji platnosti
modelu; kovariance ze singulárního rozkladu; opakovaná proložení s faktory a jittery podle zdroje,
robustní vyřazování s opětovným zařazením a faktor červeného šumu; kalibrace světelných křivek
a věrnost sítí s mapováním výstupního prostoru; kloubové řádky jednostranných dat o sklonu;
prodlužování řetězců emcee na cílový počet $\tau$ a oprava apriorního rozdělení $1/e$; a předpovědi
s jejich pásy.
