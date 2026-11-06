# Řešené příklady {#examples}

Příklady v této kapitole postupují od jedné jednoduché soustavy k ansámblům, proložením, rozptylu
a vyvíjejícím se hvězdám. Konfigurace jsou v `configs/examples/` (mapa chaosu v
`configs/validation/`); parametry příkladů 1–3 jsou ilustrativní, nepatří žádné skutečné soustavě.
Uváděná čísla jsou ta, která příkazy vypíšou; doby výpočtu platí pro jedno jádro stolního počítače.

## Příklad 1: precedující čtyřhvězda 2+2 {#ex-2p2}

Konfigurace `configs/examples/example_2p2.yaml` je soustavou z kapitoly
[Rychlý start](#quick-start): dvojhvězdu A (2,2 M☉, 6,74 d) vidíme téměř z boku ($i = 85°$),
dvojhvězda B (1,7 M☉, 2,68 d) je skloněna o 60° k obloze a obě kolem sebe obíhají každých 500 d po
excentrické dráze ($e = 0{,}3$). Z úhlů plynou vzájemné sklony asi 18° mezi `out` a A a 29° mezi
`out` a B.

```bash
swingersim run configs/examples/example_2p2.yaml -o output/example_2p2.h5
```

| Řádek souhrnu | Hodnota | Význam |
| --- | --- | --- |
| `out a`, `out e` | 1,940–1,942 au; 0,2995–0,3007 | perioda 500 d s 3,9 M☉ dává $a_\mathrm{out} = 1{,}94$ au (třetí Keplerův zákon); rozptyl 0,1 % je oskulační zákmit, ne skutečná změna |
| `A a`, `A e` | 0,0908 au; 0,044–0,054 | 19,5 R☉, stálá na 0,07 %; malá oscilace $e$ vynucená vnější dráhou; žádné ZLK, protože 18° je pod 39,2° |
| `i_mut out-A`, `out-B` | 16,3–18,1°; 28,5–32,9° | každá dvojhvězda preceduje kolem vnější dráhy a zachovává si svůj náklon |
| `i_mut A-B` | 25,6–46,7° | A a B precedují různou rychlostí, takže se úhel mezi jejich rovinami mění |
| max \|dE/E\| | $3{,}2\times10^{-13}$ | integrace je přesná |
| `precession period` | 283 let (A), 598 let (B), odhad | běh na 200 let je kratší než jeden cyklus; běh na 1000 let naměří 255 a 590 let |
| `stability ratio` | 3,67 (Mardling–Aarseth), 4,52 (Vynatheya+22) pro A; 6,80 a 8,46 pro B | bezpečně stabilní |
| `binaries A and B` | $P_A/P_B = 2{,}515$, nejbližší 5:2 (+0,6 %) | periody nejsou v rezonanci |
| `eclipses` | A 2 %, B 43 %, obě 0 % | v těchto 200 letech A zakrývá jen krátce na začátku, B téměř polovinu doby |
| `binary A` | perioda apsidálního stáčení 25 660 let ze samotné dvojhvězdy, naměřeno 270 let | bez slapů stáčí periastrum izolované dvojhvězdy jen relativita; ve čtyřhvězdě ho vnější dráha nechá obíhat stokrát rychleji |
| doba výpočtu | 80 s | škáluje s nejkratší periodou a délkou běhu |

**Tabulka 41.** Souhrn příkladu 1, řádek po řádku. Vlevo: řádek souhrnu; uprostřed: vypsaná
hodnota; vpravo: její význam.

![**Obrázek 26.** Vzájemné sklony tří dvojic drah příkladu 1 za 1000 let (`run.t_max_yr: 1000`).
Vodorovná osa: čas v letech; svislá osa: vzájemný sklon ve stupních. Modře: A–B; oranžově: out–B;
fialově: out–A. Úhel mezi A a B osciluje mezi 17,8° a 46,7° s periodou asi 450 let; úhel B vůči
vnější dráze osciluje mezi asi 29° a 34° a úhel A mezi asi 16° a 18°, ve stejném
rytmu.](../../figures/example_2p2_mutual.png)

Cyklus úhlu mezi A a B s periodou 450 let (obrázek 26) je záznějem obou period precese,
$1/(1/255 - 1/590) \approx 450$ let: obě roviny jsou střídavě na stejné a na opačných stranách
vnější dráhy, takže se úhel mezi nimi kývá mezi rozdílem a součtem jejich náklonů k vnější dráze.

![**Obrázek 27.** Sklon A (modře) a B (oranžově) k rovině oblohy za 1000 let. Vodorovná osa: čas
v letech; svislá osa: sklon ve stupních. Vybarvené pásy jsou rozsahy, v nichž každá dvojhvězda
zakrývá ($\lvert i - 90°\rvert < 6{,}0°$ pro A, $< 9{,}1°$ pro B, z výchozích poloměrů). A se kývá
mezi 54° a 88° a zakrývá jen blízko svých maxim; B se kývá mezi 42° a
102°.](../../figures/example_2p2_sky_inclination.png)

Obrázek 27 ukazuje, co z precese vidí pozorovatel: dvojhvězda zakrývá jen tehdy, je-li její sklon
uvnitř jejího pásu, takže A zakrývá v oknech dlouhých několik desetiletí každých 255 let a dnes
($t = 0$) se v takovém okně nachází, zatímco B zakrývá asi třetinu doby. Uvnitř každého okna se
hloubka zákrytů plynule mění. Tentýž běh s `--backend secular` trvá 2,5 s místo 80 s a dá vzájemné
sklony s přesností na několik setin stupně.

## Příklad 2: cykly ZLK {#ex-zlk}

`configs/examples/example_2p2_zlk.yaml` má všechny uzly souhlasné, takže vzájemné sklony jsou
rozdíly sklonů: B je skloněna o 65° vůči vnější dráze (200 d), uvnitř okna ZLK 39,2–140,8°; A jen
o 5°.

```bash
swingersim run configs/examples/example_2p2_zlk.yaml -o output/zlk.h5
```

```text
  B in out: i_mut0 = 65.0 deg, g1_0 = 90 deg
    quadrupole theory : e = 0.010..0.807, i_mut = 39.7..65.0 deg, e-cycle 18.9 yr (inside ZLK window)
    simulation        : e = 0.007..0.811, i_mut = 40.6..65.1 deg, 6 ZLK cycles of ~14 yr
    stability ratio Mardling-Aarseth = 1.64 (stable), Vynatheya+22 = 1.73 (stable), eps_oct = 0.002, eps_SA = 0.032  (mild hierarchy: secular estimates approximate)
```

![**Obrázek 28.** Cykly ZLK dvojhvězdy B v příkladu 2. Nahoře: excentricita B (oranžově) a A (modře)
v závislosti na čase v letech; čárkovaná čára je největší excentricita v kvadrupólové aproximaci
testovací částice, $\sqrt{1-\frac53\cos^2 65°} = 0{,}84$. Dole: vzájemný sklon B (oranžově) a A
(fialově) vůči vnější dráze ve stupních; tečkovaná čára označuje 39,2°, dolní okraj okna ZLK.
V každém cyklu vzroste excentricita B na 0,81, zatímco její sklon klesne asi na 41°, a pak se obě
veličiny vrátí.](../../figures/example_zlk.png)

Rozsah excentricity souhlasí s kvadrupólovou teorií středovanou přes oběh, cyklus je však asi
o 25 % kratší, než teorie předpovídá (14 let proti 18,9 roku; obrázek 28). Souhrn na to upozorňuje:
$\varepsilon_\mathrm{SA} = 0{,}032$ znamená mírnou hierarchii (vnější perioda je jen 20krát delší
než perioda B), pro niž středované rovnice ztrácejí přesnost (Luo et al. 2016); při dvakrát širší
vnější dráze souhlasí středovaná perioda na 1,5 %. V maximu excentricity se hvězdy B přiblíží na
$a(1-e) = 4{,}4$ R☉, tedy na trojnásobek součtu jejich poloměrů.

## Příklad 3: hvězdy s vnitřní stavbou {#ex-tides}

`configs/examples/example_2p2_tides.yaml` je kompaktní soustava 2+2 ze dvou zákrytových dvojhvězd
s periodami 1,8 a 2,9 dne, s `physics: {gr: true, tides: true, rotation: true}`: hvězdy jsou
rozlehlá tělesa deformovaná slapy a rotací, rotující synchronně se svou dvojhvězdou.

```bash
swingersim run configs/examples/example_2p2_tides.yaml -o output/example_tides.h5
```

| Řádek | Hodnota | Význam |
| --- | --- | --- |
| `tidal+rotational apsidal period` (A) | 64 let proti času ZLK 949 let | deformace stáčí periastrum A patnáctkrát rychleji, než by vnější dráha stačila vybudit excentricitu, takže cykly ZLK dvojhvězdy A jsou potlačeny |
| `stars` | $k_2 = 0{,}0144$, $P_\mathrm{rot}$ 1,8 → 1,8 d | stavba z výchozích polytrop; bez disipace se rychlosti rotace nemění |
| `Roche lobe filled up to` | 0,29–0,51 | oddělené dvojhvězdy |
| `binary A` | očekáváno 62,3 roku, naměřeno 59,4 roku | převládají vlastní hvězdy; poruchy od B a vnější dráhy přidávají několik procent |
| `binary B` | očekáváno 1116 let, naměřeno 505 let | B je širší, její vlastní deformace je slabá; druhá dvojhvězda rychlost zdvojnásobí |
| max \|dE/E\| | $2{,}2\times10^{-11}$ | zachovává se, protože slapové zpoždění je nulové |

**Tabulka 42.** Klíčové řádky souhrnu příkladu 3.

Ve čtyřhvězdě nemusí pozorovaná perioda apsidálního stáčení zákrytové dvojhvězdy souhlasit s hodnotou
očekávanou jen ze stavby hvězd (tabulka 42). Zadání slapového časového zpoždění
(`stars: {m3: {lag_s: 10}}`) zapne slapovou disipaci; souhrn pak uvede časové škály zakulacení
a synchronizace. Pro studie apsidálního stáčení hvězd o hmotnosti 2–3 M☉ používejte místo výchozích polytrop
`stellar_model.kind: mist` nebo `sse` (hodnoty $k_2$ podle Clareta 2023); polytropy dávají pro
takové hvězdy $k_2$ asi třikrát větší.

## Příklad 4: do minulosti a stabilita {#ex-past}

Tatáž soustava jako v příkladu 1, integrovaná 200 let do minulosti, potřebuje jen záporný koncový
čas, `run: {t_max_yr: -200}` nebo `--t-max-yr -200`:

```text
integrated to t = -200 yr (2000 outputs, osculating elements)
  i_mut A-B          17.82..  45.32 deg
  eclipses (observer along z): A 27 %, B 0 %, both at once 0 % of the run
```

V 200 letech před dneškem zakrývala A asi čtvrtinu doby a B vůbec; v příštích 200 letech je tomu
naopak. Kontrola stability příkladu 2:

```bash
swingersim stability configs/examples/example_2p2_zlk.yaml
```

```text
nested triples (ratio > 1 = stable):
  A in out (+B)   a_in/a_out = 0.0706, e_in = 0.010, e_out = 0.300, i_mut = 5.0 deg, q_out = 0.773: Mardling-Aarseth 2.38, Vynatheya+22 3.04
  B in out (+A)   a_in/a_out = 0.1029, e_in = 0.010, e_out = 0.300, i_mut = 65.0 deg, q_out = 1.294: Mardling-Aarseth 1.64, Vynatheya+22 1.73
N-body ghost test (100 outer orbits, Newtonian): stable after 100.0 outer orbits, max relative divergence 1.70e-03
```

Kritičtější je dvojhvězda B, skloněná o 65° a procházející cykly ZLK; obě kritéria i přímý test se
shodují, že soustava je stabilní.

## Příklad 5: ansámbl – V994 Her {#ex-ensemble}

U V994 Her určuje fotometrie periody, excentricity a sklony obou dvojhvězd (84° a 85°), nikoli však
to, jak jsou dráhy natočeny kolem zorného paprsku. `configs/examples/example_v994_ensemble.yaml`
vylosuje 200 soustav: hmotnosti v rozmezí 2 %, sklony v rozmezí 1°, dvojhvězdu B s 85° nebo 95°
(smysl oběhu neznámý), její uzel vůči A kdekoli a vnější dráhu s náhodnou orientací; každý člen se
integruje sekulárním backendem 10 000 let dopředu i dozadu.

```bash
swingersim ensemble configs/examples/example_v994_ensemble.yaml -o output/v994/ensemble -j 4
swingersim ensemble-report output/v994/ensemble --plot output/v994/overview.png --x orbits.out.i --y B_didt0_deg_yr
```

```text
forward: 200 members (fractions with 68 % Wilson intervals)
  fate intact           100.0 % (99.5-100.0) [200/200]
  i_mut A_B     now 90.3 deg (16-84 %: 27.7-150.9), range over the run: median 40.0...141.9 deg, members that cross 90 deg: 69.5 %
  A      di/dt now +0.00462 deg/yr (16-84 %: -0.0765..+0.0723); of the members eclipsing now, eclipses end within the run in 92 % (median after 440.4 yr)
  eclipses (time averaged over members): A 48.9 %, B 51.6 %, both 27.9 %
  stable by Vynatheya+22       100.0 % (99.5-100.0) [200/200]
```

Těch 400 běhů trvá s Numba asi 10 minut na jednom jádře. Soustava je bezpečná v každém členu a její
zákryty jsou dočasné: u 92 % členů zákryty A do 10 000 let skončí, typicky asi po 400 letech, a obě
dvojhvězdy zakrývají současně jen 28 % doby – dnes pozorovaný dvojitě zákrytový stav je jen fází.
Současná rychlost změny sklonu závisí na neznámé orientaci vnější dráhy, takže naměřená rychlost
tuto orientaci omezuje; přesně to dělá podmínění z kapitoly [Pracovní postupy](#wf-ensemble).
Minulost vypadá statisticky stejně jako budoucnost, jak tomu bez naměřeného směru změny musí být.

## Příklad 6: proložení okamžiků minim, krok 1 {#ex-fit1}

`configs/examples/example_fit.yaml` je soustava 2+2, jejíž minima jsou syntetická: tabulka
`example_fit_minima.csv` (50 minim A, 354 minim B za 17,7 roku; vytvořená skriptem
`tools/make_example_fit_minima.py`) vznikla s vnější dráhou 2200 d ($e_3 = 0{,}35$,
$\omega_3 = 1{,}2$ rad), apsidálním stáčením A o 0,0015 rad za cyklus a $a\sin i$ 1,1 au (A) a 0,95 au
(B), v náhodně vybraných cyklech a s chybami 2·10⁻⁴ d. Proložení používá efekt světelného času pro oba
páry, apsidální stáčení pro A a jednu společnou vnější dráhu:

```bash
swingersim fit configs/examples/example_fit.yaml
```

```text
stage 1: chi^2 = 385.8 for 391 degrees of freedom (reduced 0.99)
  out.P3        2199.8 ± 1.2       d
  out.e3        0.3508 ± 0.0058
  out.w3         1.191 ± 0.015     rad
  A.e          0.24986 ± 0.00040
  A.dw       0.0015010 ± 0.0000028 rad/cycle
  A.asini       1.0965 ± 0.0073    AU
  B.asini       0.9461 ± 0.0032    AU
```

Všechny hodnoty se obnoví v rámci 1,3 směrodatné odchylky a redukované $\chi^2$ je 0,99, jak má být
u chyb, které jsou z konstrukce správné. Cykly jsou náhodné záměrně: minima vzdálená o pevný počet $N$
cyklů nerozliší rychlost stáčení $\dot\omega$ od $\dot\omega + 2\pi k/N$ (alias). Dřívější verze této
tabulky měla minima A každých 40 cyklů a proložení vrátilo 0,1586 rad za cyklus ($0{,}0015 + 2\pi/40$)
s redukovaným $\chi^2$ rovným 14. Do stejné pasti mohou spadnout i skutečná data s pravidelným
odstupem (například jedno minimum za sezónu); varovným znamením je redukované $\chi^2$ daleko nad 1
při chybách, kterým věříme.

## Příklad 7: setkání dvou dvojhvězd {#ex-scatter}

Může čtyřhvězda 2+2 vzniknout ze dvou dvojhvězd, které se setkají v mladé hvězdokupě?
`configs/examples/example_scatter_2p2.yaml` posílá dvojhvězdu s $a = 0{,}5$ au rychlostí 1 km/s
proti dvojhvězdě s $a = 1$ au (všechny hvězdy 1 M☉ a 1 R☉, excentricity losované z tepelného
rozdělení, pericentra až 3 au):

```bash
swingersim scatter configs/examples/example_scatter_2p2.yaml -o output/scatter_2p2 --n 200 -j 0
```

```text
  outcome                               n  fraction    1-sigma interval  sigma [AU^2]
  collision                            97     0.485  [0.450, 0.520]  3.245e+04
  flyby                                56     0.280  [0.249, 0.313]  1.874e+04
  exchange + ionization                19     0.095  [0.076, 0.118]  6357
  triple + single (stable)             12     0.060  [0.045, 0.079]  4015
  single ionization                    10     0.050  [0.037, 0.068]  3346
  exchange                              2     0.010  [0.005, 0.020]  669.2
  triple + single (unstable)            2     0.010  [0.005, 0.020]  669.2
  bound 3+1 (unstable)                  1     0.005  [0.002, 0.013]  334.6
  bound 2+2 (unstable)                  1     0.005  [0.002, 0.013]  334.6
```

Sloupce uvádějí výsledek (collision – srážka, flyby – průlet, exchange – výměna, ionization –
ionizace, triple + single – trojhvězda a osamocená hvězda, bound – vázaná soustava), počet setkání,
podíl s jeho 68% intervalem a účinný průřez. Při 1 km/s jsou obě dvojhvězdy po setkání silně
vázané: téměř každé těsné setkání se změní v dlouhou rezonanční interakci a v polovině z nich se dvě
hvězdy dotknou. Nejčastějším dynamickým produktem je stabilní trojhvězda s unikající hvězdou (6 %);
vázaná soustava 2+2 nebo 3+1 se objevuje jen jako dočasný, nestabilní stav. Zachycení do stabilních
čtyřhvězd vyžaduje třetího účastníka nebo disipaci (slapy, plyn). Těch 200 setkání trvalo 22 minut
na dvou jádrech.

## Příklad 8: kde je soustava stabilní? Mapa chaosu {#ex-chaos}

`configs/validation/quarles18_circumbinary.yaml` umístí testovací částici kolem dvojhvězdy se
stejnými hmotnostmi složek ($a = 1$ au); mapa prochází její velkou poloosu a excentricitu:

```bash
swingersim chaosmap configs/validation/quarles18_circumbinary.yaml -o output/chaosmap_cb \
  --x orbits.out.a --x-range 1.8 4.0 --nx 23 --y orbits.out.e --y-range 0 0.5 --ny 11 --t-max 3000 -j 0
```

![**Obrázek 29.** Mapa MEGNO cirkumbinární dráhy po 3000 letech. Vodorovná osa: velká poloosa vnější
dráhy v au; svislá osa: její excentricita. Barva: MEGNO (stupnice vpravo), 2 pro pravidelný pohyb
(modře), nad 4 chaotický, 10 označuje únik nebo srážku (tmavě červeně). Hranice mezi pravidelnou
oblastí (vpravo dole) a úniky (vlevo nahoře) sleduje téměř stálou vzdálenost pericentra
$a(1-e)$.](../../figures/chaosmap_circumbinary.png)

Těch 253 běhů trvá 7 minut na dvou jádrech. Při $e = 0$ je dráha pro tuto jedinou počáteční fázi
pravidelná od 2,0 au dále; publikovaná mez, která uvažuje všechny fáze, je reprodukována
v dokumentu *Testování, část I*. Týž příkaz zmapuje libovolnou dvojici parametrů čtyřhvězdy, např.
`orbits.out.P` proti `orbits.B.P` blízko poměru period dvojhvězd 3:2.

## Příklad 9: co by naměřili pozorovatelé? Okamžiky zákrytů {#ex-timings}

`configs/examples/example_timings_czev4315.yaml` je dvojitě zákrytová čtyřhvězda CzeV4315 s jednou
možnou vzájemnou dráhou (120 d, blízko signálu 119 d v okamžicích zákrytů B; orientace je
ilustrativní):

```bash
swingersim timings configs/examples/example_timings_czev4315.yaml --years 5 \
  -o output/timings_czev4315.csv --plot output/timings_czev4315.png
```

```text
eclipse timings over 5 yr
  A: 45 primary and 73 secondary eclipses, P = 6.74571156 d, O-C of primaries -840.62 .. +787.62 min (semi-amplitude 814.12 min)
    eclipses missing between 0.330 and 3.158 yr
    no primary eclipses after t = 3.635 yr
  B: 1987 primary and 1988 secondary eclipses, P = 0.91885817 d, O-C of primaries -4.30 .. +3.59 min (semi-amplitude 3.94 min)
```

![**Obrázek 30.** Syntetické O−C soustavy CzeV4315 pro vzájemnou dráhu s periodou 120 d za 5 let.
Nahoře: dvojhvězda A, dole: dvojhvězda B; vodorovná osa: čas v letech od $t = 0$; svislá osa: O−C
v minutách; modré body jsou primární, oranžové sekundární zákryty. B ukazuje vlnu ±4 min s periodou
120 d. A zakrývá jen ve dvou intervalech; její primární a sekundární zákryty jsou od sebe vzdáleny
tisíce minut a driftují.](../../figures/timings_czev4315.png)

B ukazuje očekávanou oscilaci s periodou 120 dní (doba letu světla plus dynamický vliv A; obrázek
30). A, excentrická a jen 18krát rychlejší než vzájemná dráha, reaguje silně: její apsidy se otočí
během let, její sekundární zákryty se posouvají asi o ±3 dny vůči primárním a téměř tři roky
nezakrývá vůbec. CzeV4315 A v datech zakrývá, takže tato konkrétní dráha je vyloučena – to je
druh testu, který ansámbl uplatní na tisíce drah.

## Příklad 10: šest hvězd – V994 Her {#ex-sextuple}

V994 Her je trojnásobně zákrytová šestihvězda: dvojitě zákrytová čtyřhvězda A+B (vzájemná dráha
1062 d) a třetí zákrytová dvojhvězda C na široké dráze (Zasche et al. 2023b). Popisuje ji předvolba
`2+2+2`; `configs/examples/example_v994_sextuple.yaml` přebírá hmotnosti, poloměry a dráhy
z publikovaného řešení a označuje, co je předpokládáno.

```bash
swingersim stability configs/examples/example_v994_sextuple.yaml --no-nbody
swingersim run configs/examples/example_v994_sextuple.yaml -o output/v994_sextuple.h5
```

Kontrola stability použije kritéria na tři vnořené trojice (AB v `out` s C: Mardling–Aarseth 21,00;
A v AB s B: 4,43; B v AB s A: 5,65) a, jsou-li nainstalovány klasifikátory, na každou dílčí
čtyřhvězdu. Běh (40 s na 30 let) naměří periody apsidálního stáčení A, B a C odpovídající 1,42, 3,51
a 0,71 °/yr proti pozorovaným 1,86, 3,66 a 1,68 °/yr: B, jejíž apsidální konstanta byla odvozena
z těchto dat, souhlasí na 4 %; A a C se s předpokládaným teoretickým $k_2$ stáčejí příliš pomalu,
takže data by omezovala právě jejich konstanty stavby nebo poloměry.

## Příklad 11: vyvíjející se čtyřhvězda a její světelná křivka {#ex-evolution}

`configs/examples/example_evolution.yaml` je ilustrativní soustava 2+2 ve stáří 1,1 miliardy let:
dvojhvězda A s primární složkou 2 M☉ (4 d) a dvojhvězda B ze dvou hvězd slunečního typu (5 d).
S `stellar_model.kind: sse` a `physics.evolution: true` sledují poloměry, zářivé výkony a konstanty
stavby stáří hvězd a `tidal_lag: hurley2002` bere slapovou disipaci každé hvězdy z její obálky.

```bash
swingersim run configs/examples/example_evolution.yaml -o output/evolution.h5
swingersim lightcurve configs/examples/example_evolution.yaml --days 30 -o output/evolution_lc.csv
```

```text
EVENT: roche_overflow at t = 7.445e+07 yr {'age_myr': 1174.45, 'stars': ['m1']}
  max |dE/E| = 2.22e-01  (tidal dissipation: energy decreases)
  max |dL/L| = 1.51e-05
  m1  ... P_rot = 4 -> 3.509 d, Roche lobe filled up to 1.01
      evolving: age 1100 -> 1174.3 Myr, R 3.72 -> 5.96 Rsun, Teff 6761 -> 5057 K, main sequence -> giant branch
```

![**Obrázek 31.** Příklad 11. Nahoře: poloměry čtyř hvězd v R☉ v závislosti na čase běhu v Myr
(stáří 1100 Myr v $t = 0$); primární složka m1 o hmotnosti 2 M☉ (červeně) opouští hlavní posloupnost
po 64 Myr a rozpíná se na 6 R☉. Uprostřed: oběžná perioda A (modře) a rotační perioda m1 (čárkovaně
červeně) ve dnech; schody jsou úseky operátorového štěpení, pokles na konci je slapová synchronizace,
která odebírá dráze moment hybnosti. Dole: světelné křivky TESS v délce 12 dní na začátku (stáří
1100 Myr, šedě) a na konci (1174 Myr, červeně, posunuto o −0,25), relativní tok v závislosti na čase
ve dnech.](../../figures/example_evolution.png)

Primární složka A přejde Hertzsprungovu mezeru a ve stáří 1174 Myr vyplní svůj Rocheův lalok; tam se
běh zastaví (`roche_overflow`), protože tento příklad nezapíná `physics.mass_transfer` (obrázek 31).
Zatímco se hvězda rozpíná, její rotace se zpomaluje; jakmile se její obálka stane konvektivní, slapy
ji během několika set tisíc let synchronizují a odeberou moment hybnosti dráze: perioda A se zkrátí
ze 4,0 na 3,5 d. Energie se proto nezachovává (slapy ji disipují), zatímco celkový moment hybnosti
se zachovává na $1{,}5\times10^{-5}$. Běh trvá 3 minuty na 74 Myr. Velká elipsoidální proměnnost
pozdní světelné křivky je na hranici aproximace malé deformace kulových hvězd;
`stars.<m>.shape: roche` dá hvězdě její Rocheův povrch.

## Příklad 12: další druhy soustav {#ex-zoo}

`configs/examples/` obsahuje připravenou konfiguraci pro každý druh soustavy; kteroukoli z nich
spustíte příkazem `swingersim run <file> -o <name>.h5`.

| Druh soustavy | Konfigurace | Co se děje |
| --- | --- | --- |
| mladé hvězdy ve svém plynu | `example_gas_2p2.yaml` | plyn zmenší, zakulatí a srovná vnější dráhu a pak se rozptýlí (sekulární backend) |
| trojhvězda se čtvrtou hvězdou | `example_3p1.yaml` | precese vnořených drah |
| obr vedle dvojhvězdy | `example_giant_triple.yaml` | třetí hvězda na AGB ztrácí hmotu větrem; vnější dráha se rozšiřuje |
| obr předávající hmotu dvojhvězdě | `example_tertiary_transfer.yaml` | třetí hvězda přenáší hmotu na vnitřní dvojhvězdu |
| přenos hmoty (Algol) | `example_algol_triple.yaml` | dárce přetéká na svého průvodce; poměr hmotností se obrátí |
| společná obálka | `example_ce_triple.yaml` | nestabilní přenos pohltí dvojhvězdu; dráha se zmenší |
| kontaktní dvojhvězda | `example_contact_evolution.yaml` | těsná dvojhvězda se dostane do kontaktu, zmenšuje se a splyne |
| pulzující hvězdy | `example_pulsating_triple.yaml`, `example_cepheid_triple.yaml` | pulzace typu δ Scuti, dlouhoperiodické a cefeidové ve světelné křivce |
| hmotné hvězdy, Wolfova–Rayetova hvězda | `example_wr_triple.yaml` | hmotná třetí hvězda se stane Wolfovou–Rayetovou hvězdou; její vítr rozšíří vnější dráhu |
| vítr, který destabilizuje | `example_mass_loss.yaml` | ztráta hmoty dostane trojhvězdu za mez stability |

**Tabulka 43.** Další ukázkové konfigurace. Vlevo: druh soustavy; uprostřed: soubor v
`configs/examples/`; vpravo: co běh ukazuje. Fyziku, kterou každá potřebuje, popisuje dokument
*Fyzika a metody*, kapitola Hvězdy ve vícenásobných soustavách.
