\appendix

# Všech 219 srovnávacích testů {#appendix-benchmarks}

Tato příloha uvádí každý srovnávací test registrovaný ve SWINGERsimu s jeho výsledkem z konečné verze
kódu, podle úrovně důkazu. Každá úroveň má dvě tabulky se stejnými řádky ve stejném pořadí: první říká,
co každý test srovnává a jaký kód ověřuje (tabulky A1, A3, A5, A7), druhá uvádí jeho výsledek (tabulky
A2, A4, A6, A8); test se v obou najde podle svého Id. V rámci úrovně jsou testy uvedeny v pořadí, v němž
je program registruje, takže řádek popsaný jako „totéž pro …“ následuje za řádkem, na který odkazuje.
Sloupce jsou:

- **Id** – identifikátor; `swingersim validate <Id>` test spustí a vypíše jeho zdroj (článek a rovnici,
  obrázek nebo tabulku), jeho vstupy a důvod jeho tolerance;
- **Kód** – třída kódu, který test ověřuje (tabulka 1): *own* – rovnice a numerika napsané pro
  SWINGERsim; *own/tab* – interpolace publikovaných tabulek ve SWINGERsimu; *R* – REBOUND IAS15
  s efekty REBOUNDx, řízený nastavením SWINGERsimu; *ext* – externí knihovna přes obal SWINGERsimu.
  Test vlastního kódu je silnějším důkazem, protože nic mimo SWINGERsim nemůže způsobit, že projde;
- **Veličina** – co se srovnává, tak jak to test popisuje (symboly jako ve zdrojovém kódu: `e`
  excentricita, `i0` počáteční sklon, `MA01` kritérium Mardlinga–Aarsetha, `tier A` kulový model
  světelné křivky);
- **Očekáváno** – referenční hodnota; **Spočteno** – hodnota, kterou spočte SWINGERsim, se všemi
  vypsanými číslicemi;
- **Rozdíl** – SWINGERsim minus očekávaná hodnota, relativně (%) nebo absolutně, podle definice testu;
- **Tolerance** – největší rozdíl, který se počítá jako shoda, stanovený před srovnáním;
- **Verdikt** – PASS (vyhovuje), nebo **EXPLAINED** (vysvětleno: mimo toleranci z doloženého
  důvodu; tabulka 48).

Při čtení řádku srovnejte sloupec *Rozdíl* se sloupcem *Tolerance*: PASS má absolutní hodnotu rozdílu
nejvýše rovnou toleranci. Je-li *Očekáváno* 0 nebo 1, test srovnává poměr nebo rozdíl s jeho ideální
hodnotou. Rychlé testy pocházejí z posledního běhu `swingersim validate`, pomalé z úplného běhu
`swingersim validate --all` na produkčním stroji a devět ručních (`phoebe-*`, `sbi-*`, `pejcha13-*`,
`vyn23-*`) ze samostatných běhů; vše s konečnou verzí kódu. Celkem: 93 testů úrovně A (93 PASS), 60
úrovně B (55 PASS, 5 EXPLAINED), 31 úrovně C (29 PASS, 2 EXPLAINED) a 35 úrovně D (28 PASS,
7 EXPLAINED); 165 z 219 ověřuje kód tříd own a own/tab.

\footnotesize

## Úroveň A: analytické výsledky {#appendix-level-a}

93 testů: 93 PASS, 0 EXPLAINED.

**Tabulka A1.** Úroveň A: co každý test srovnává a jaký kód ověřuje (sloupce *Id*, *Kód*, *Veličina*). Výsledky následují v tabulce A2.

| Id | Kód | Veličina |
|---|---|---|
| `gr-mercury` | R | stáčení perihelia Merkuru v řádu 1PN |
| `zlk-emax-65` | R | maximální excentricita ZLK, i0 = 65 deg, testovací částice |
| `zlk-emax-80` | R | maximální excentricita ZLK, i0 = 80 deg, testovací částice |
| `zlk-imin-65` | R | vzájemný sklon v maximu excentricity, i0 = 65 deg |
| `zlk-below-critical` | R | maximální excentricita pod kritickým úhlem, i0 = 30 deg |
| `zlk-kozai-constant` | R | sekulární drift sqrt(1-e\^2) cos i během cyklů ZLK, i0 = 65 deg |
| `stability-unstable` | R | 2+2 při 0,5 × meze Mardlinga–Aarsetha změní hierarchii (1 = ano) |
| `stability-stable` | R | 2+2 při 2,5 × meze Mardlinga–Aarsetha si hierarchii uchová (0 = beze změny) |
| `quad-tp-emax` | own | sekulární kvadrupólový kód: e\_max ZLK, testovací částice, i0 = 65 deg |
| `zlk-period-nbody` | own | perioda ZLK, N-částicově vůči kvadrupólové sekulární periodě (i0 = 65 deg) |
| `megno-regular` | R | MEGNO pravidelné (kvaziperiodické) hierarchické trojhvězdy |
| `megno-chaotic` | R | MEGNO > 4 pro 2+2 uvnitř nestabilní oblasti (1 = ano) |
| `brouwer-law` | R | exponent růstu chyby energie IAS15 (Brouwerův zákon) |
| `nodal-precession` | R | perioda precese uzlů, kruhová vnitřní dráha, i = 20 deg |
| `polytrope-k2-n1.5` | own | apsidální konstanta polytropy n = 1,5 |
| `polytrope-k2-n3` | own | apsidální konstanta polytropy n = 3 |
| `roche-eggleton-q1` | own | poloměr Rocheova laloku / vzdálenost složek, q = 1 (Eggletonův vzorec proti přesné hodnotě) |
| `roche-eggleton-q0.1` | own | poloměr Rocheova laloku / vzdálenost složek, q = 0,1 (Eggletonův vzorec proti přesné hodnotě) |
| `apsidal-sterne-static` | R | apsidální perioda ze slapové deformace, e = 0,3, bez rotace |
| `apsidal-sterne-rotating` | R | apsidální perioda ze slapové + rotační deformace, e = 0,3, Omega = n |
| `tides-leconte-a` | R | slapové da/dt, N-částicově / teorie, e = 0,3, nerotující hvězdy |
| `tides-leconte-e` | R | slapové de/dt, N-částicově / teorie, e = 0,3, nerotující hvězdy |
| `tides-leconte-spin` | R | slapové roztáčení dOmega/dt, N-částicově / teorie, e = 0,6, Omega = 1,2 n |
| `tides-leconte-e-eccentric` | R | slapové de/dt, N-částicově / teorie, e = 0,6, Omega = 1,2 n |
| `secular-kozai-emax` | own | sekulární backend: maximální excentricita ZLK, i0 = 65 deg, testovací částice |
| `secular-nodal-precession` | own | sekulární backend: perioda precese uzlů, i = 20 deg \[roky\] |
| `secular-gr-apsidal` | own | sekulární backend: relativistická apsidální perioda, a = 0,05 au, e = 0,3 |
| `secular-3p1-kozai-emax` | own | sekulární backend 3+1: e\_max ZLK vnitřní dráhy, testovací částice, i0 = 65 deg |
| `gas-migration-a` | own | plynná fáze: velká poloosa po 3000 letech s tau\_a = 5000 let \[au\] |
| `gas-damping-e` | own | plynná fáze: excentricita po 3000 letech s tau\_e = 2000 let |
| `gas-alignment-i` | own | plynná fáze: sklon k disku po 3000 letech, i0 = 40 deg, tau\_i = 1500 let |
| `scatter-hyperbola-deflection` | own | nastavení rozptylu: úhel odklonu hyperbolického průletu 1 + 0,5 M☉, v\_inf = 5 km/s, b = 3 au |
| `escape-theory-energy-slope-L0` | own | statistická teorie úniku, stejné hmotnosti, L0 = 0: sklon dP/dE\_B pro /E\_B/ >> /E0/ |
| `escape-theory-energy-slope-highL` | own | statistická teorie úniku, L0 = L\_c: sklon dP/dE\_B |
| `escape-theory-mean-e-L0` | own | statistická teorie úniku, L0 = 0: střední excentricita dvojhvězdy |
| `escape-theory-mean-e-highL` | own | statistická teorie úniku, L0 = L\_c: střední excentricita dvojhvězdy |
| `timing-ltte-amplitude` | own | syntetický diagram O-C: poloviční amplituda světelného času dvojhvězdy A, kruhová vnější dráha viděná z boku, a = 10 au, stejné hmotnosti |
| `phenom-romer-delay` | own | fenomenologická minima: zpoždění sekundárního minima za polovinou periody, 2 + 1 M☉ ZAMS, P = 5 d, kruhová dráha, i = 89 deg |
| `phoebe-romer-delay` | ext | hybrid PHOEBE: zpoždění sekundárního minima za polovinou periody, 2 + 1 M☉, P = 5 d, kruhová dráha (se světelným časem minus bez něj) |
| `pulsation-time-delay-circular` | own | amplituda zpoždění daných dobou šíření světla u pulzující hvězdy v kruhové dvojhvězdě (1,7 + 0,5 M☉, P = 100 d, i = 60 deg), změřená ze syntetické světelné křivky jako data Kepleru (kadence 30 min, fáze pulzací ve dvoudenních úsecích) |
| `mt-soberman-general` | own | vzdálenost složek po přenosu s konstantními podíly (dárce 3 -> 1,5 M☉, příjemce 2 M☉; 30 % vítr dárce, 40 % izotropní reemise, 30 % akreováno): numerická odezva dráhy proti uzavřenému tvaru |
| `mt-conservative` | own | vzdálenost složek po konzervativním přenosu (dárce 3 -> 1 M☉, příjemce 2 -> 4) |
| `mt-isotropic-reemission` | own | vzdálenost složek po přenosu s úplnou izotropní reemisí (dárce 3 -> 1 M☉, příjemce 2 M☉) |
| `ce-impulsive-outer-a` | R | vnější velká poloosa po společné obálce uvnitř trojhvězdy (N-částicově) vůči hodnotě předtím, dělená předpovědí impulzní ztráty hmoty pro kruhovou dráhu |
| `ce-impulsive-outer-e` | R | vnější excentricita po společné obálce téže trojhvězdy minus předpověď impulzní ztráty hmoty |
| `matt15-unsaturated-powerlaw` | own | Omega(4,55 Gyr) / Omega(1,14 Gyr) bržděného modelu Slunce |
| `galactic-tide-emax-60` | own | největší excentricita dráhy 2·10⁴ au skloněné o 60 deg ke galaktickému disku (sekulární backend, physics.galactic\_tide) |
| `flyby-first-order-kick` | own | relativní rozdíl změny vektoru excentricity dvojhvězdy při parabolickém průletu s Q = 20 a: vzorec prvního řádu proti přímé N-částicové integraci |
| `sbi-linear-gaussian-std` | ext | SBI (neuronový odhad aposteriorního rozdělení): šířka aposteriorního rozdělení lineárního gaussovského modelu |
| `quint-out-zlk-emax` | R | 2+2+1: e\_max ZLK dráhy obou dvojhvězd buzené pátou hvězdou, i0 = 65 deg |
| `mist-grid-tams-age-1` | own/tab | mřížka MIST: stáří 1 M☉ na konci hlavní posloupnosti (EEP 454), \[Fe/H\] = 0, v/vcrit = 0 |
| `mist-interpolation` | own/tab | MIST: poloměr 1,14 M☉ v TAMS interpolovaný z 1,12 a 1,16 M☉ / tabelovaný |
| `claret2023-k2-2msun` | own/tab | mřížka Claret (2023): log k2 pro 2,0 M☉, \[Fe/H\] = 0, log g = 4,069 |
| `claret2023-beta-1msun` | own/tab | mřížka Claret (2023): poloměr setrvačnosti beta pro 1,0 M☉, \[Fe/H\] = 0, log g = 4,216 |
| `sse-zams-radius-1` | own/tab | mřížka vývojových drah SSE: poloměr ZAMS pro 1 M☉, Z = 0,02 |
| `sse-zams-radius-2` | own/tab | mřížka vývojových drah SSE: poloměr ZAMS pro 2,5 M☉, Z = 0,02 |
| `ml-adiabatic-invariant` | R | max abs(a M / (a M)\_0 - 1), zatímco m1 pomalu ztrácí 0,3 M☉ (Psi = 5,6·10⁻⁶), N-částicově |
| `ml-osculating-e-03` | R | max abs(e - e\_theory) / Psi během dvou oběhů, e0 = 0,3, Psi = 5,6·10⁻⁵, N-částicově |
| `ml-osculating-e-07` | R | max abs(e - e\_theory) / Psi během dvou oběhů, e0 = 0,7, Psi = 5,6·10⁻⁵, N-částicově |
| `ml-circular-emax` | R | největší e / Psi počátečně kruhové dvojhvězdy při pomalé ztrátě hmoty |
| `ml-impulsive-a` | R | a'/a poté, co kruhová dvojhvězda ztratí 15 % hmotnosti během 10⁻³ periody |
| `ml-impulsive-e` | R | excentricita poté, co kruhová dvojhvězda ztratí 15 % hmotnosti během 10⁻³ P |
| `ml-zlk-nbody-secular` | own | okamžik 4. vzestupu ZLK (e > 0,6), zatímco rušící těleso ztrácí polovinu hmotnosti, N-částicově / sekulární kvadrupól |
| `mist-interpolation-massive-r` | own/tab | MIST: poloměr 60 M☉ uprostřed hlavní posloupnosti (EEP 330) interpolovaný z 55 a 65 M☉ / tabelovaný |
| `mist-interpolation-massive-m` | own/tab | MIST: současná hmotnost 100 M☉ v TAMS interpolovaná z 95 a 105 M☉ / tabelovaná |
| `mist-vink2001` | own/tab | vítr MIST pro 40 M☉ ve 30 % hlavní posloupnosti / Vink et al. (2001) rov. 24 |
| `claret2004-log-k2` | own/tab | Claret (2004): log k2 pro 39,8 M☉ při log g = 3,416 (3,65 Myr) |
| `claret2004-beta` | own/tab | Claret (2004): poloměr setrvačnosti beta pro 39,8 M☉ při log g = 3,416 |
| `ld-tlusty-row` | own/tab | okrajové ztemnění TLUSTY: a1 modelu OG40000g400v10 v Bessellově V |
| `ocfit-gimenez-bastero` | own | apsidální člen okamžiků minim, e = 0,3: největší rozdíl mezi přesnými okamžiky konjunkcí a řadou Giméneze & Bastera (1995) použitou v OCFit |
| `ocfit-lite-amplitude` | own | polovina rozkmitu efektu světelného času během jednoho oběhu / K podle Irwina (1952), e3 = 0,55 |
| `ocfit-nbody-asini-A` | own | proložení N-částicových okamžiků zákrytů čtyřhvězdy 2+2 v kroku 1 (vnější P = 1500 d, e = 0,3, i = 60 deg): proložené / keplerovské a sin i páru A |
| `ocfit-nbody-asini-B` | own | totéž pro pár B (v protifázi, společná vnější dráha) |
| `ocfit-nbody-P3` | own | totéž proložení: proložená / nastavená vnější perioda |
| `wind-transfer-jeans` | own | dvojhvězda MIST 25 + 18 M☉ s větrem a přenosem hmoty: největší relativní změna a (M1 + M2) před začátkem přenosu |
| `secular-nested-tides-hut81` | own | 3+1 ve vnořeném sekulárním backendu: slapové de/dt vnitřní dvojhvězdy (P = 4 d, e = 0,05, synchronní) / Hut (1981) rov. 10 |
| `galactic-tide-emax-60-nbody` | R | největší excentricita dvojhvězdy 2·10⁴ au skloněné o 60 deg k disku, N-částicový backend se slapy Galaxie (od e = 0,01) |
| `gw-secular-circular-decay` | own | kruhová dvojhvězda 1,4 + 1,4 M☉ v trojhvězdě, sekulární backend s physics.gw: a / a0 (1 - t/tau)\^(1/4) v 90 % doby do splynutí |
| `gw-nbody-decay` | R | tatáž dvojhvězda s e = 0,3, N-částicový backend: změna a během 10 let / Petersova rychlost × 10 let |
| `disc-strip-uniform` | own | neprůhledná deska viděná z boku (poloviční tloušťka 0,3 R) přes hvězdu s rovnoměrnou jasností: propuštěný podíl |
| `disc-face-on-tau` | own | disk viděný shora s optickou hloubkou 0,7 před hvězdou: propuštěný podíl |
| `disc-inclined-tau` | own | disk s optickou hloubkou 0,7 skloněný o 60 deg před hvězdou: propuštěný podíl |
| `roche-lobe-volume` | own | poloměr ekvivalentního objemu Rocheova laloku (q = M2/M1 = 0,3) z povrchové sítě Rocheova modulu / Eggletonův vzorec |
| `rotating-star-critical` | own | rovníkový / polární poloměr hvězdy v Rocheově modelu při kritické rotaci |
| `rotating-star-equator` | own | R\_eq/R\_pol - (1 + w\^2/2) s w = Omega / Omega\_K(R\_eq), Rocheův model, Omega / Omega\_crit = 0,6 |
| `mesh-sphere-eclipse` | own | největší rozdíl světelných křivek dvou koulí během obou zákrytů: povrchové sítě (modul sítí) proti přesným zákrytům tier A |
| `disc-ring-precession` | own | rychlost precese úzkého tuhého prstence ve vzdálenosti 4 a\_b kolem kruhové dvojhvězdy / rychlost testovací částice |
| `disc-polar-critical` | own | kritický sklon \[deg\] pro polární libraci disku kolem dvojhvězdy s e = 0,5, z integrací rovnice tuhého disku ve SWINGERsimu |
| `rejuvenation-hurley` | own | stáří příjemce 2 M☉ na hlavní posloupnosti poté, co získá 1 % své hmotnosti / (t'\_MS / t\_MS)(M / M') t |
| `k2-rotation` | own | rotační oprava k2 (2 M☉, 2 R☉, V = 200 km/s) / 10\^(-0,87 × 2 V\^2 / 3 g R) (Claret 1999) |
| `contact-jdot` | own | doba ztráty momentu hybnosti \[Gyr\] kontaktní dvojhvězdy 1,0 + 0,5 M☉, P = 0,35 d (jednotky SWINGERsimu proti cgs) |
| `pdot-parabola` | R | kvadratický koeficient okamžiků zákrytů dvojhvězdy s dP/dt = 4·10⁻⁶ d/rok během 60 let (N-částicově s impulzy periody) / analyticky (1/2) P dP/dt |
| `dynamic-onset-circularisation` | own | vyplnění Rocheova laloku dárce po zakulacení dráhy Kozaiovy–Lidovovy trojhvězdy v okamžiku dotyku (e\_c = 0,89), krát (1 + e\_c) |

**Tabulka A2.** Úroveň A: výsledky (sloupce *Očekáváno*, *Spočteno*, *Rozdíl*, *Tolerance*, *Verdikt*) ve stejném pořadí jako v tabulce A1.

| Id | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|
| `gr-mercury` | 42,98 ″/století | 42,9803 ″/století | +0,00 % | 0,5 % | PASS |
| `zlk-emax-65` | 0,838047 | 0,84023 | +0,26 % | 1 % | PASS |
| `zlk-emax-80` | 0,974548 | 0,976155 | +0,16 % | 1 % | PASS |
| `zlk-imin-65` | 39,2315° | 39,3441° | +0,113° | 0,5° | PASS |
| `zlk-below-critical` | 0 | 0,00173037 | +0,00173 | 0,01 | PASS |
| `zlk-kozai-constant` | 0 | 5,16811·10⁻⁵ | +5,17·10⁻⁵ | 0,001 | PASS |
| `stability-unstable` | 1 | 1 | +0 | 0 | PASS |
| `stability-stable` | 0 | 0 | +0 | 0 | PASS |
| `quad-tp-emax` | 0,838047 | 0,838082 | +0,00 % | 0,1 % | PASS |
| `zlk-period-nbody` | 0 | 0,00660493 | +0,0066 | 0,02 | PASS |
| `megno-regular` | 2 | 1,99959 | −0,00041 | 0,05 | PASS |
| `megno-chaotic` | 1 | 1 | +0 | 0 | PASS |
| `brouwer-law` | 0,5 | 0,59301 | +0,093 | 0,15 | PASS |
| `nodal-precession` | 11351,4 roku | 11339,8 roku | −0,10 % | 1 % | PASS |
| `polytrope-k2-n1.5` | 0,1433 | 0,143279 | −0,01 % | 0,2 % | PASS |
| `polytrope-k2-n3` | 0,0144 | 0,014443 | +0,30 % | 0,5 % | PASS |
| `roche-eggleton-q1` | 0,379861 | 0,378921 | −0,25 % | 1 % | PASS |
| `roche-eggleton-q0.1` | 0,2054 | 0,206773 | +0,67 % | 1 % | PASS |
| `apsidal-sterne-static` | 318,227 roku | 318,201 roku | −0,01 % | 0,5 % | PASS |
| `apsidal-sterne-rotating` | 291,213 roku | 291,185 roku | −0,01 % | 0,5 % | PASS |
| `tides-leconte-a` | 1 | 0,999285 | −0,07 % | 1 % | PASS |
| `tides-leconte-e` | 1 | 0,999356 | −0,06 % | 1 % | PASS |
| `tides-leconte-spin` | 1 | 1,00113 | +0,11 % | 1 % | PASS |
| `tides-leconte-e-eccentric` | 1 | 0,996923 | −0,31 % | 1 % | PASS |
| `secular-kozai-emax` | 0,838047 | 0,838082 | +0,00 % | 0,5 % | PASS |
| `secular-nodal-precession` | 11351,4 roku | 11351,4 roku | −0,00 % | 0,5 % | PASS |
| `secular-gr-apsidal` | 7113,77 roku | 7113,77 roku | −0,00 % | 0,5 % | PASS |
| `secular-3p1-kozai-emax` | 0,838047 | 0,838082 | +0,00 % | 0,5 % | PASS |
| `gas-migration-a` | 5,48812 au | 5,48812 au | −0,00 % | 0,0001 % | PASS |
| `gas-damping-e` | 0,111565 | 0,111565 | +0,00 % | 0,0001 % | PASS |
| `gas-alignment-i` | 5,63999° | 5,63943° | −0,000561° | 0,01° | PASS |
| `scatter-hyperbola-deflection` | 173,548° | 173,548° | +1,36·10⁻¹²° | 1·10⁻⁶° | PASS |
| `escape-theory-energy-slope-L0` | −3,5 | −3,5026 | −0,0026 | 0,05 | PASS |
| `escape-theory-energy-slope-highL` | −4 | −4,02479 | −0,0248 | 0,1 | PASS |
| `escape-theory-mean-e-L0` | 0,8133 | 0,812817 | −0,000483 | 0,01 | PASS |
| `escape-theory-mean-e-highL` | 0,7 | 0,710776 | +0,0108 | 0,02 | PASS |
| `timing-ltte-amplitude` | 41,5837 min | 41,5863 min | +0,01 % | 0,1 % | PASS |
| `phenom-romer-delay` | 27,4517 s | 27,3578 s | −0,34 % | 1,5 % | PASS |
| `phoebe-romer-delay` | 27,452 s | 27,287 s | −0,60 % | 2 % | PASS |
| `pulsation-time-delay-circular` | 53,859 s | 53,859 s | +0,00 % | 0,01 % | PASS |
| `mt-soberman-general` | 1,29573 a/a0 | 1,29573 a/a0 | −0,00 % | 1·10⁻⁵ % | PASS |
| `mt-conservative` | 2,25 a/a0 | 2,25 a/a0 | −0,00 % | 1·10⁻⁵ % | PASS |
| `mt-isotropic-reemission` | 2,03003 a/a0 | 2,03003 a/a0 | −0,00 % | 1·10⁻⁵ % | PASS |
| `ce-impulsive-outer-a` | 1 | 1,01274 | +1,27 % | 3 % | PASS |
| `ce-impulsive-outer-e` | 0 | 0,00442155 | +0,00442 | 0,02 | PASS |
| `matt15-unsaturated-powerlaw` | 0,5 | 0,491949 | −0,00805 | 0,02 | PASS |
| `galactic-tide-emax-60` | 0,829156 | 0,829181 | +2,5·10⁻⁵ | 0,01 | PASS |
| `flyby-first-order-kick` | 0 | 0,00594917 | +0,00595 | 0,03 | PASS |
| `sbi-linear-gaussian-std` | 0,3092 | 0,291 | −5,8 % | 10 % | PASS |
| `quint-out-zlk-emax` | 0,81272 | 0,8152 | +0,31 % | 1,0 % | PASS |
| `mist-grid-tams-age-1` | 9919,42 Myr | 9919,42 Myr | +0,00 % | 0,001 % | PASS |
| `mist-interpolation` | 1 | 1,01068 | +1,07 % | 2 % | PASS |
| `claret2023-k2-2msun` | −2,478 | −2,47794 | +6·10⁻⁵ | 0,002 | PASS |
| `claret2023-beta-1msun` | 0,2265 | 0,226499 | −6,34·10⁻⁷ | 0,0005 | PASS |
| `sse-zams-radius-1` | 0,88825 R☉ | 0,888312 R☉ | +0,01 % | 0,1 % | PASS |
| `sse-zams-radius-2` | 1,79434 R☉ | 1,79435 R☉ | +0,00 % | 0,1 % | PASS |
| `ml-adiabatic-invariant` | 0 | 3,56804·10⁻⁹ | +3,57·10⁻⁹ | 1·10⁻⁶ | PASS |
| `ml-osculating-e-03` | 0 | 0,0126471 | +0,0126 | 0,05 | PASS |
| `ml-osculating-e-07` | 0 | 0,0123538 | +0,0124 | 0,05 | PASS |
| `ml-circular-emax` | 2 | 2,00265 | +0,13 % | 2 % | PASS |
| `ml-impulsive-a` | 1,21429 | 1,21491 | +0,05 % | 1 % | PASS |
| `ml-impulsive-e` | 0,176471 | 0,176892 | +0,24 % | 1 % | PASS |
| `ml-zlk-nbody-secular` | 1 | 1,0021 | +0,21 % | 1,0 % | PASS |
| `mist-interpolation-massive-r` | 1 | 0,999431 | −0,06 % | 1 % | PASS |
| `mist-interpolation-massive-m` | 1 | 1,00107 | +0,11 % | 3 % | PASS |
| `mist-vink2001` | 1 | 1,05511 | +5,51 % | 10 % | PASS |
| `claret2004-log-k2` | −2,81523 | −2,81569 | −0,000463 | 0,002 | PASS |
| `claret2004-beta` | 0,161706 | 0,161699 | −6,61·10⁻⁶ | 0,001 | PASS |
| `ld-tlusty-row` | 1,18262 | 1,18262 | +0 | 0,0001 | PASS |
| `ocfit-gimenez-bastero` | 0 s | 0,000629837 s | +0,00063 s | 0,01 s | PASS |
| `ocfit-lite-amplitude` | 1 | 1 | −0,00 % | 0,0001 % | PASS |
| `ocfit-nbody-asini-A` | 1 | 0,99864 | −0,14 % | 1,0 % | PASS |
| `ocfit-nbody-asini-B` | 1 | 1,0024 | +0,24 % | 1,0 % | PASS |
| `ocfit-nbody-P3` | 1 | 0,99996 | −0,00 % | 0,5 % | PASS |
| `wind-transfer-jeans` | 0 | 1,68718·10⁻¹¹ | +1,69·10⁻¹¹ | 1·10⁻⁸ | PASS |
| `secular-nested-tides-hut81` | 1 | 1 | −0,00 % | 0,0001 % | PASS |
| `galactic-tide-emax-60-nbody` | 0,82916 | 0,82908 | −0,01 % | 1,0 % | PASS |
| `gw-secular-circular-decay` | 1 | 1 | −0,00 % | 0,01 % | PASS |
| `gw-nbody-decay` | 1 | 1,0003 | +0,03 % | 1,0 % | PASS |
| `disc-strip-uniform` | 0,623838 | 0,62475 | +0,000912 | 0,005 | PASS |
| `disc-face-on-tau` | 0,496585 | 0,496585 | +0 | 1·10⁻⁹ | PASS |
| `disc-inclined-tau` | 0,246597 | 0,246597 | −4,16·10⁻¹⁶ | 1·10⁻⁹ | PASS |
| `roche-lobe-volume` | 1 | 1,00509 | +0,51 % | 1 % | PASS |
| `rotating-star-critical` | 1,5 | 1,5 | −0,00 % | 0,0001 % | PASS |
| `rotating-star-equator` | 0 | −2,22045·10⁻¹⁶ | −2,22·10⁻¹⁶ | 1·10⁻⁹ | PASS |
| `mesh-sphere-eclipse` | 0 | 0,000259643 | +0,00026 | 0,001 | PASS |
| `disc-ring-precession` | 1 | 1,00018 | +0,02 % | 0,1 % | PASS |
| `disc-polar-critical` | 37,7612° | 37,7612° | +6,29·10⁻⁹° | 0,2° | PASS |
| `rejuvenation-hurley` | 1 | 0,999174 | −0,08 % | 0,2 % | PASS |
| `k2-rotation` | 1 | 1 | +0,00 % | 1·10⁻⁷ % | PASS |
| `contact-jdot` | 2,18893 Gyr | 2,18898 Gyr | +0,00 % | 0,1 % | PASS |
| `pdot-parabola` | 1 | 1,00005 | +0,00 % | 0,1 % | PASS |
| `dynamic-onset-circularisation` | 1 | 1 | −0,00 % | 0,01 % | PASS |

## Úroveň B: publikované numerické výsledky {#appendix-level-b}

60 testů: 55 PASS, 5 EXPLAINED.

**Tabulka A3.** Úroveň B: co každý test srovnává a jaký kód ověřuje (sloupce *Id*, *Kód*, *Veličina*). Výsledky následují v tabulce A4.

| Id | Kód | Veličina |
|---|---|---|
| `naoz13-fig4-imin` | R | Naoz et al. 2013, obr. 4 (uspořádání PSR B1620-26): minimální vzájemný sklon |
| `naoz13-fig4-imax` | R | Naoz et al. 2013, obr. 4: maximální vzájemný sklon (překlopení dráhy na retrográdní) |
| `naoz13-fig3-flip-time` | R | Naoz et al. 2013, obr. 3 (excentrické ZLK): okamžik prvního překlopení přes 90 deg |
| `naoz13-fig3-imax` | R | Naoz et al. 2013, obr. 3: maximální sklon po překlopení |
| `quad-naoz13-fig4-imin` | own | sekulární kvadrupólový kód: Naoz et al. 2013, obr. 4, minimální sklon |
| `quad-naoz13-fig4-imax` | own | sekulární kvadrupólový kód: Naoz et al. 2013, obr. 4, maximální sklon |
| `vok16-v994-JA` | own | Vokrouhlický 2016, uspořádání V994 Her: počáteční sklon J\_A dvojhvězdy A k vnější dráze |
| `vok16-v994-JB` | own | Vokrouhlický 2016, uspořádání V994 Her: počáteční sklon J\_B dvojhvězdy B k vnější dráze |
| `vok16-jab-max-analytic` | own | Vokrouhlický 2016, analytická mez J\_AB\^max pro V994 Her (uspořádání obr. 2) |
| `vok16-fig2-jab-min` | own | sekulární backend: Vokrouhlický 2016, obr. 2, minimální J\_AB během 100 kyr |
| `vok16-fig2-jab-max` | own | sekulární backend: Vokrouhlický 2016, obr. 2, maximální J\_AB během 100 kyr |
| `vok16-fig2-jab-1myr` | own | sekulární backend: uspořádání obr. 2, maximální J\_AB během 1 Myr |
| `vok16-fig3-jab-max` | own | sekulární backend: Vokrouhlický 2016, obr. 3, maximální J\_AB během 20 kyr |
| `vok16-eclipse-prob-A` | own | V994 Her: pravděpodobnost, že dvojhvězda A zakrývá, izotropní orientace |
| `vok16-eclipse-prob-B` | own | V994 Her: pravděpodobnost, že dvojhvězda B zakrývá, izotropní orientace |
| `secular-naoz13-fig3-flip-time` | own | sekulární backend 3+1: Naoz et al. 2013, obr. 3, okamžik prvního překlopení \[roky\] |
| `secular-naoz13-fig3-imax` | own | sekulární backend 3+1: Naoz et al. 2013, obr. 3, maximální sklon po překlopení |
| `at16-binbin-single-ionization` | R | rozptyl dvojhvězda–dvojhvězda, stejné hmotnosti, a1 = a2, e = 0, v\_inf = v\_crit: účinný průřez „prosté ionizace“ / (pi a\^2) |
| `at16-binbin-exchange-ionization` | R | rozptyl dvojhvězda–dvojhvězda, stejné hmotnosti, a1 = a2, e = 0, v\_inf = v\_crit: účinný průřez „výměny + ionizace“ / (pi a\^2) |
| `at16-binbin-exchange` | R | rozptyl dvojhvězda–dvojhvězda, stejné hmotnosti, a1 = a2, e = 0, v\_inf = v\_crit: účinný průřez „výměny“ / (pi a\^2) |
| `chaosmap-circumbinary-limit` | R | mapa MEGNO: vnitřní mez stability cirkumbinární dráhy, kruhová dvojhvězda se stejnými hmotnostmi, všech 8 fází pravidelných \[a\_bin\] |
| `qzcar-broz-a1` | own | QZ Car: velká poloosa Ac1+Ac2 (P = 6,0 d) z elementů a hmotností tabulky 2 |
| `qzcar-broz-a2` | own | QZ Car: velká poloosa Aa1+Aa2 (P = 20,7 d) z elementů a hmotností tabulky 2 |
| `qzcar-broz-a3` | own | QZ Car: velká poloosa Ac+Aa (P = 14 728 d) z elementů a hmotností tabulky 2 |
| `mt-zeta-lobe-q2` | own | odezva zeta\_L = dln R\_L / dln M\_d Rocheova laloku dárce (Eggleton 1983) na konzervativní přenos při q = 2 |
| `ce-tout97-period` | own | oběžná perioda po společné obálce předchůdce kataklyzmické proměnné podle Touta et al. (1997): hvězda AGB 5,36 M☉ (jádro 1,014 M☉) s průvodcem 1 M☉, alpha\_CE = 1 |
| `ce-bse-algol-second-period` | own | oběžná perioda po druhé společné obálce v příkladu Algolu z BSE (obr 2,85 M☉ s jádrem 0,416 M☉, bílý trpaslík 0,417 M☉, P = 18 d, alpha\_CE = 3) |
| `matt15-sun-period` | own | perioda rotace modelu Slunce ve stáří 4,55 Gyr s momentem síly větru podle Matta et al. (2015), začátek při P = 1 d |
| `applegate-huaqr` | own | energie potřebná pro Applegateův mechanismus vůči energii vyzářené za jednu periodu modulace, HU Aqr (model tenké slupky) |
| `applegate-nnser` | own | energie potřebná pro Applegateův mechanismus vůči energii vyzářené za jednu periodu modulace, NN Ser (model tenké slupky) |
| `applegate-v471tau` | own | energie potřebná pro Applegateův mechanismus vůči energii vyzářené za jednu periodu modulace, V471 Tau (model tenké slupky) |
| `vyn23-2p2-ma01-score` | own | celkové skóre kritéria MA01 proti N-částicovému testu duchem na trénovací množině 2+2 SWINGERsimu (vzorkováno jako ve Vynatheya et al. 2023, oddíl 3.1) |
| `vyn23-2p2-mlp-score` | ext | skóre emulátoru stability MLP SWINGERsimu na odložených členech (2+2; 11 parametrů, 4 × 50 neuronů) proti N-částicovému testu duchem |
| `vyn23-3p1-ma01-score` | own | celkové skóre kritéria MA01 proti N-částicovému testu duchem na trénovací množině 3+1 SWINGERsimu (vzorkováno jako ve Vynatheya et al. 2023, oddíl 3.1) |
| `vyn23-3p1-mlp-score` | ext | skóre emulátoru stability MLP SWINGERsimu na odložených členech (3+1; 11 parametrů, 4 × 50 neuronů) proti N-částicovému testu duchem |
| `escape-theory-vs-nbody-light-escaper` | own | rezonanční setkání dvojhvězdy a osamocené hvězdy (dvojhvězda 1 + 0,6 M☉, hvězda 0,8 M☉): N-částicový podíl, v němž uniká hvězda 0,6 M☉, proti statistické teorii (očekávaná hodnota: vlastní vyčíslení teorie Stonea & Leigha 2019 ve SWINGERsimu, nikoli publikované číslo) |
| `czev343-apsidal-tidal` | own | CzeV343, dvojhvězda A: slapová rychlost apsidálního stáčení \[rad/den\] |
| `czev343-apsidal-rotational` | own | CzeV343, dvojhvězda A: rotační rychlost apsidálního stáčení, pseudosynchronní rotace \[rad/den\] |
| `czev343-apsidal-gr` | own | CzeV343, dvojhvězda A: relativistická rychlost apsidálního stáčení \[rad/den\] |
| `vok16-pejcha-fA` | own | populace 2+2 podle Pejchy et al., sekulárně: dvojhvězda A dosáhne 1 - e < 10⁻³ |
| `vok16-pejcha-fB` | own | populace 2+2 podle Pejchy et al., sekulárně: dvojhvězda B dosáhne 1 - e < 10⁻³ |
| `pejcha13-F6-fA` | R | populace F6 Pejchy et al., N-částicově: dvojhvězda A dosáhne 1 - e < 10⁻³ |
| `pejcha13-F6-fB` | R | populace F6 Pejchy et al., N-částicově: dvojhvězda B dosáhne 1 - e < 10⁻³ |
| `mist-ms-lifetime-1` | own/tab | doba života 1 M☉ na hlavní posloupnosti (ZAMS až TAMS), SSE proti MIST |
| `mist-ms-lifetime-2` | own/tab | doba života 2 M☉ na hlavní posloupnosti (ZAMS až TAMS), SSE proti MIST |
| `mist-tams-radius-2` | own/tab | poloměr 2 M☉ na konci hlavní posloupnosti, SSE proti MIST |
| `iw-two-layer` | own | disipace inerciálních vln dvouvrstvé hvězdy (alpha 0,6, beta 0,95): impulzní formalismus / uzavřený vzorec |
| `tides-circularization-pms` | own | perioda, pod níž jsou dvojhvězdy dvou hvězd 1 M☉ (e0 = 0,3) kruhové při dosažení hlavní posloupnosti (MIST před hlavní posloupností od 1 Myr do 50 Myr; rovnovážný slap a inerciální vlny) |
| `ml-tedi-onset` | R | počátek vnitřní nestability (e\_in > 0,1) / překročení meze MA01, vnitřní primární složka ztrácí 3 z 5 M☉ za 8000 let |
| `mist-geneva-tauh-15-0` | own/tab | doba hoření vodíku 15 M☉, bez rotace, MIST proti ženevským modelům |
| `mist-geneva-tauh-15-4` | own/tab | doba hoření vodíku 15 M☉, s rotací, MIST proti ženevským modelům |
| `mist-geneva-tauh-60-0` | own/tab | doba hoření vodíku 60 M☉, bez rotace, MIST proti ženevským modelům |
| `mist-geneva-tauh-60-4` | own/tab | doba hoření vodíku 60 M☉, s rotací, MIST proti ženevským modelům |
| `mist-geneva-tauh-120-0` | own/tab | doba hoření vodíku 120 M☉, bez rotace, MIST proti ženevským modelům |
| `mist-geneva-tauh-120-4` | own/tab | doba hoření vodíku 120 M☉, s rotací, MIST proti ženevským modelům |
| `flyby-nbody-backend` | R | parabolický průlet s Q = 20 a integrovaný N-částicovým backendem (physics.encounters): relativní rozdíl impulzu excentricity od výsledku prvního řádu podle Hamerse & Samsinga (2019) |
| `tidal-capture-radius-n3` | own | největší vzdálenost periastra pro slapový záchyt / R, dvě hvězdy n = 3 (0,6 M☉, 0,65 R☉) při 10 km/s |
| `tidal-capture-radius-n15` | own | totéž pro dvě hvězdy n = 1,5 |
| `tidal-capture-radius-n15-mcmillan` | own | největší vzdálenost periastra pro slapový záchyt / R, dvě stejné hvězdy n = 1,5 při 10 km/s |
| `mt-bse-algol-giant-phase` | own | perioda na konci přenosu v příkladu Algolu z BSE, spočtená z jeho stavu na konci Hertzsprungovy mezery (0,53 + 2,74 M☉, 17 d; konzervativně) |

**Tabulka A4.** Úroveň B: výsledky (sloupce *Očekáváno*, *Spočteno*, *Rozdíl*, *Tolerance*, *Verdikt*) ve stejném pořadí jako v tabulce A3.

| Id | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|
| `naoz13-fig4-imin` | 57,5° | 57,2391° | −0,261° | 0,5° | PASS |
| `naoz13-fig4-imax` | 106,7° | 106,93° | +0,23° | 0,5° | PASS |
| `naoz13-fig3-flip-time` | 3,8·10⁶ let | 3,9868·10⁶ let | +4,92 % | 10,0 % | PASS |
| `naoz13-fig3-imax` | 142,6° | 143,01° | +0,408° | 2° | PASS |
| `quad-naoz13-fig4-imin` | 57,5° | 57,5376° | +0,0376° | 0,1° | PASS |
| `quad-naoz13-fig4-imax` | 106,7° | 106,721° | +0,0207° | 0,1° | PASS |
| `vok16-v994-JA` | 11,6° | 11,6462° | +0,0462° | 0,1° | PASS |
| `vok16-v994-JB` | 15,8° | 15,7932° | −0,00678° | 0,1° | PASS |
| `vok16-jab-max-analytic` | 5,2° | 5,20074° | +0,000736° | 0,05° | PASS |
| `vok16-fig2-jab-min` | 3,5° | 3,47152° | −0,0285° | 0,15° | PASS |
| `vok16-fig2-jab-max` | 5,2° | 5,20069° | +0,000687° | 0,15° | PASS |
| `vok16-fig2-jab-1myr` | 5,2° | 5,2006° | +0,000633° | 0,15° | PASS |
| `vok16-fig3-jab-max` | 6° | 6,20874° | +0,209° | 0,5° | PASS |
| `vok16-eclipse-prob-A` | 0,366 | 0,369265 | +0,89 % | 2 % | PASS |
| `vok16-eclipse-prob-B` | 0,397 | 0,399075 | +0,52 % | 2 % | PASS |
| `secular-naoz13-fig3-flip-time` | 3,8·10⁶ let | 4,04981·10⁶ let | +6,57 % | 10 % | PASS |
| `secular-naoz13-fig3-imax` | 142,6° | 142,556° | −0,0435° | 2° | PASS |
| `at16-binbin-single-ionization` | 2,645 | 2,6873 | +1,60 % | 10,0 % | PASS |
| `at16-binbin-exchange-ionization` | 3,259 | 3,2305 | −0,88 % | 10,0 % | PASS |
| `at16-binbin-exchange` | 0,302 | 0,37165 | +23,06 % | 30,0 % | PASS |
| `chaosmap-circumbinary-limit` | 2,27 | 2,3 | +0,03 | 0,06 | PASS |
| `qzcar-broz-a1` | 54,032 R☉ | 53,9946 R☉ | −0,07 % | 0,2 % | PASS |
| `qzcar-broz-a2` | 136,132 R☉ | 136,081 R☉ | −0,04 % | 0,2 % | PASS |
| `qzcar-broz-a3` | 13053 R☉ | 13046 R☉ | −0,05 % | 0,2 % | PASS |
| `mt-zeta-lobe-q2` | 2,59 | 2,60707 | +0,0171 | 0,05 | PASS |
| `ce-tout97-period` | 35 h | 34,84 h | −0,46 % | 10 % | PASS |
| `ce-bse-algol-second-period` | 0,04 d | 0,0365638 d | −8,59 % | 15 % | PASS |
| `matt15-sun-period` | 27,97 d | 27,7635 d | −0,74 % | 2 % | PASS |
| `applegate-huaqr` | 0,1 | 0,105659 | +5,66 % | 10 % | PASS |
| `applegate-nnser` | 3,3 | 3,30338 | +0,10 % | 10 % | PASS |
| `applegate-v471tau` | 0,014 | 0,0143688 | +2,63 % | 10 % | PASS |
| `vyn23-2p2-ma01-score` | 0,83 | 0,869 | +0,039 | 0,05 | PASS |
| `vyn23-2p2-mlp-score` | 0,94 | 0,9375 | −0,0025 | 0,03 | PASS |
| `vyn23-3p1-ma01-score` | 0,56 | 0,7896 | +0,23 | 0,05 | **EXPLAINED** |
| `vyn23-3p1-mlp-score` | 0,93 | 0,9195 | −0,0105 | 0,03 | PASS |
| `escape-theory-vs-nbody-light-escaper` | 0,769 | 0,75369 | −0,0153 | 0,07 | PASS |
| `czev343-apsidal-tidal` | 0,000414 rad/d | 0,000422272 rad/d | +2,00 % | 3 % | PASS |
| `czev343-apsidal-rotational` | 6,4·10⁻⁵ rad/d | 6,32533·10⁻⁵ rad/d | −1,17 % | 5 % | PASS |
| `czev343-apsidal-gr` | 1,47·10⁻⁵ rad/d | 1,47874·10⁻⁵ rad/d | +0,59 % | 2 % | PASS |
| `vok16-pejcha-fA` | 0,31 | 0,21 | −0,1 | 0,07 | **EXPLAINED** |
| `vok16-pejcha-fB` | 0,481 | 0,4025 | −0,0785 | 0,07 | **EXPLAINED** |
| `pejcha13-F6-fA` | 0,358 | 0,24 | −0,118 | 0,14 | PASS |
| `pejcha13-F6-fB` | 0,530 | 0,41 | −0,12 | 0,15 | PASS |
| `mist-ms-lifetime-1` | 9878 Myr | 11003,5 Myr | +11,39 % | 15 % | PASS |
| `mist-ms-lifetime-2` | 1054 Myr | 1164,17 Myr | +10,45 % | 15 % | PASS |
| `mist-tams-radius-2` | 3,913 R☉ | 3,56583 R☉ | −8,87 % | 15 % | PASS |
| `iw-two-layer` | 1 | 1,00031 | +0,03 % | 0,5 % | PASS |
| `tides-circularization-pms` | 7,85 d | 7,9307 d | +1,03 % | 10 % | PASS |
| `ml-tedi-onset` | 1 | 1,06475 | +6,47 % | 15 % | PASS |
| `mist-geneva-tauh-15-0` | 11,015 Myr | 12,3538 Myr | +12,15 % | 15 % | PASS |
| `mist-geneva-tauh-15-4` | 13,447 Myr | 12,6825 Myr | −5,69 % | 15 % | PASS |
| `mist-geneva-tauh-60-0` | 3,53 Myr | 3,71791 Myr | +5,32 % | 15 % | PASS |
| `mist-geneva-tauh-60-4` | 4,465 Myr | 3,94657 Myr | −11,61 % | 15 % | PASS |
| `mist-geneva-tauh-120-0` | 2,672 Myr | 2,78283 Myr | +4,15 % | 15 % | PASS |
| `mist-geneva-tauh-120-4` | 3,178 Myr | 3,33128 Myr | +4,82 % | 15 % | PASS |
| `flyby-nbody-backend` | 0 | 0,00738403 | +0,00738 | 0,05 | PASS |
| `tidal-capture-radius-n3` | 2,72 | 2,36599 | −13,01 % | 5 % | **EXPLAINED** |
| `tidal-capture-radius-n15` | 3,84 | 3,33454 | −13,16 % | 5 % | **EXPLAINED** |
| `tidal-capture-radius-n15-mcmillan` | 3,4 | 3,33454 | −1,93 % | 5 % | PASS |
| `mt-bse-algol-giant-phase` | 30 d | 30,3561 d | +1,19 % | 3 % | PASS |

## Úroveň C: jiné programy {#appendix-level-c}

31 testů: 29 PASS, 2 EXPLAINED.

**Tabulka A5.** Úroveň C: co každý test srovnává a jaký kód ověřuje (sloupce *Id*, *Kód*, *Veličina*). Výsledky následují v tabulce A6.

| Id | Kód | Veličina |
|---|---|---|
| `sextuple-secular-emax` | own | šestihvězda 2+2+2: maximální excentricita A, sekulární backend proti N-částicovému |
| `sextuple-secular-zlk` | own | šestihvězda 2+2+2: perioda ZLK dvojhvězdy A, sekulární backend proti N-částicovému |
| `sextuple-nbody-emax` | R | šestihvězda 2+2+2: maximální excentricita A, N-částicově (reprodukovatelnost) |
| `mse-3p1-inner-emax` | own | 3+1 (a = 1, 15, 300 au, i\_mid = 85 deg): maximální e vnitřní dvojhvězdy během 10⁵ let |
| `mse-3p1-inner-t-emax` | own | 3+1: okamžik maximální e vnitřní dvojhvězdy |
| `mse-3p1-mid-emax` | own | 3+1: maximální e prostřední dráhy |
| `mse-2p2-A-emax` | own | 2+2 (a = 1, 1,5, 40 au): maximální e dvojhvězdy A během 8·10⁴ let |
| `mse-2p2-B-emax` | own | 2+2: maximální e dvojhvězdy B během 8·10⁴ let |
| `mse-3p1-rlof-time` | own | 3+1 (a = 0,3, 5, 200 au, i = 89 deg): počátek přetečení Rocheova laloku vnitřní dvojhvězdy vyvolaného cykly ZLK |
| `phoebe-sphere-vs-tier-a` | own | největší rozdíl normovaných světelných křivek, tier A proti PHOEBE 2.5 se stejnou fyzikou (koule, černá tělesa, Claretovo okrajové ztemnění), e = 0,1 |
| `qzcar-broz-lv-ac1` | own | QZ Car: podíl Ac1 na světle v pásmu V z jejího poloměru a Teff |
| `qzcar-broz-lv-ac2` | own | QZ Car: podíl Ac2 na světle v pásmu V z jejího poloměru a Teff |
| `qzcar-broz-lv-aa1` | own | QZ Car: podíl Aa1 na světle v pásmu V z jejího poloměru a Teff |
| `qzcar-broz-lv-aa2` | own | QZ Car: podíl Aa2 na světle v pásmu V z jejího poloměru a Teff |
| `cepheid-period-radius-6` | own/tab | log R cefeidy MIST 6 M☉ (modrá smyčka, střed pásu nestability, perioda podle Bona et al. 2000) minus vztah perioda–poloměr rotujících ženevských modelů cefeid |
| `cepheid-period-radius-7` | own/tab | jako cepheid-period-radius-6 pro cefeidu 7 M☉ |
| `mt-bse-algol-donor` | own | hmotnost dárce na konci prvního přenosu hmoty v příkladu Algolu z BSE (2,9 + 0,9 M☉, Z = 0,02, začátek při P = 8 d, e = 0,7; od počátku přetečení, kruhová dráha, P = 2,7 d) |
| `mt-bse-algol-accretor` | own | hmotnost příjemce na konci prvního přenosu hmoty v příkladu Algolu z BSE |
| `mt-bse-algol-period` | own | oběžná perioda na konci prvního přenosu hmoty v příkladu Algolu z BSE |
| `mt-cosmic-algol-donor` | own | hmotnost dárce na konci prvního přenosu hmoty v příkladu Algolu z BSE proti kódu BSE v COSMIC 4.2.1 (2,9 + 0,9 M☉, kruhová dráha 2,7 d, Z = 0,02, acc\_lim = -1) |
| `lc-occultation-small` | own | tok hvězdy s okrajovým ztemněním za tmavým diskem, rp = 0,1, centrálně |
| `lc-occultation-limb` | own | tok hvězdy s okrajovým ztemněním za tmavým diskem, rp = 0,3 u okraje (z = 0,9) |
| `lc-occultation-large` | own | tok hvězdy s okrajovým ztemněním za větším tmavým diskem, rp = 1,2, z = 0,5 |
| `lc-detached-primary` | own | hloubka primárního zákrytu oddělené dvojhvězdy z koulí |
| `lc-detached-secondary` | own | hloubka sekundárního zákrytu oddělené dvojhvězdy z koulí |
| `lc-ellipsoidal` | own | poloviční amplituda elipsoidální proměnnosti, R/a = 0,15, q = 0,8, i = 60 deg |
| `ld-tlusty-vs-atlas` | own/tab | intenzita středovaná přes disk <I>/I(1) ve V při 25 000 K, log g 4: TLUSTY bez LTE / ATLAS LTE |
| `sse-wind-rgb-1msun` | own/tab | hmotnost ztracená hvězdou 1 M☉ (Z = 0,02) na větvi rudých obrů s Reimersovým větrem, eta = 0,1: vývojová dráha SSE s větrem SSE proti modelu MIST v1.2 (MESA) na vrcholu větve obrů |
| `contact-vs-phoebe` | own | amplituda ve V kontaktní dvojhvězdy (q = 0,42, přeplnění 0,3, i = 80 deg): síť SWINGERsimu minus PHOEBE 2.5 se stejnou fyzikou |
| `rv-light-vs-phoebe` | own | rychlost středu světla minus rychlost těžiště sekundární složky kontaktní dvojhvězdy (q = 0,42, přeplnění 0,3, i = 80 deg): největší rozdíl SWINGERsim - PHOEBE / její amplituda (52 km/s) |
| `dynamic-onset-nbody-vs-secular` | own | okamžik dotyku Rocheova laloku vyvolaného Kozaiovými–Lidovovými cykly (první maximum, \~61 let): N-částicově / sekulárně (dvojitě středováno, oktupól) |

**Tabulka A6.** Úroveň C: výsledky (sloupce *Očekáváno*, *Spočteno*, *Rozdíl*, *Tolerance*, *Verdikt*) ve stejném pořadí jako v tabulce A5.

| Id | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|
| `sextuple-secular-emax` | 0,75207 | 0,751821 | −0,03 % | 2 % | PASS |
| `sextuple-secular-zlk` | 70,2 roku | 71,04 roku | +1,20 % | 5 % | PASS |
| `sextuple-nbody-emax` | 0,75207 | 0,75207 | +0,00 % | 0,2 % | PASS |
| `mse-3p1-inner-emax` | 0,9707 | 0,970714 | +1,37·10⁻⁵ | 0,002 | PASS |
| `mse-3p1-inner-t-emax` | 49325 let | 49322,8 roku | −0,00 % | 0,2 % | PASS |
| `mse-3p1-mid-emax` | 0,2245 | 0,224461 | −3,86·10⁻⁵ | 0,002 | PASS |
| `mse-2p2-A-emax` | 0,9494 | 0,949492 | +9,2·10⁻⁵ | 0,002 | PASS |
| `mse-2p2-B-emax` | 0,7758 | 0,776227 | +0,000427 | 0,002 | PASS |
| `mse-3p1-rlof-time` | 380,03 roku | 380,1 roku | +0,02 % | 0,2 % | PASS |
| `phoebe-sphere-vs-tier-a` | 0 | 0,00037 | +0,00037 | 0,001 | PASS |
| `qzcar-broz-lv-ac1` | 0,22 | 0,220775 | +0,000775 | 0,005 | PASS |
| `qzcar-broz-lv-ac2` | 0,133 | 0,133856 | +0,000856 | 0,005 | PASS |
| `qzcar-broz-lv-aa1` | 0,635 | 0,635024 | +2,36·10⁻⁵ | 0,005 | PASS |
| `qzcar-broz-lv-aa2` | 0,01 | 0,0103458 | +0,000346 | 0,005 | PASS |
| `cepheid-period-radius-6` | 0 dex | −0,00765757 dex | −0,00766 dex | 0,05 dex | PASS |
| `cepheid-period-radius-7` | 0 dex | 0,00392438 dex | +0,00392 dex | 0,05 dex | PASS |
| `mt-bse-algol-donor` | 0,42 M☉ | 0,404336 M☉ | −0,0157 M☉ | 0,05 M☉ | PASS |
| `mt-bse-algol-accretor` | 2,85 M☉ | 2,51105 M☉ | −0,339 M☉ | 0,15 M☉ | **EXPLAINED** |
| `mt-bse-algol-period` | 30 d | 21,9301 d | −26,90 % | 15 % | **EXPLAINED** |
| `mt-cosmic-algol-donor` | 0,401 M☉ | 0,404336 M☉ | +0,00334 M☉ | 0,03 M☉ | PASS |
| `lc-occultation-small` | 0,987891 | 0,987891 | −1,92·10⁻¹² | 1·10⁻⁷ | PASS |
| `lc-occultation-limb` | 0,945236 | 0,945236 | +3,66·10⁻⁹ | 1·10⁻⁷ | PASS |
| `lc-occultation-large` | 0,128979 | 0,128979 | −5,52·10⁻⁹ | 1·10⁻⁷ | PASS |
| `lc-detached-primary` | 0,380821 | 0,380868 | +4,72·10⁻⁵ | 0,0001 | PASS |
| `lc-detached-secondary` | 0,144658 | 0,144654 | −3,85·10⁻⁶ | 0,0001 | PASS |
| `lc-ellipsoidal` | 0,0027949 | 0,00260619 | −6,75 % | 10 % | PASS |
| `ld-tlusty-vs-atlas` | 1 | 0,997587 | −0,24 % | 2 % | PASS |
| `sse-wind-rgb-1msun` | 0,04158 M☉ | 0,0421789 M☉ | +1,44 % | 10 % | PASS |
| `contact-vs-phoebe` | 0 mag | 0,0030377 mag | +0,00304 mag | 0,03 mag | PASS |
| `rv-light-vs-phoebe` | 0 | 0,022028 | +0,022 | 0,05 | PASS |
| `dynamic-onset-nbody-vs-secular` | 1 | 1,0173 | +1,73 % | 3,0 % | PASS |

## Úroveň D: pozorované soustavy {#appendix-level-d}

35 testů: 28 PASS, 7 EXPLAINED.

**Tabulka A7.** Úroveň D: co každý test srovnává a jaký kód ověřuje (sloupce *Id*, *Kód*, *Veličina*). Výsledky následují v tabulce A8.

| Id | Kód | Veličina |
|---|---|---|
| `v994s-dyn-A` | R | šestihvězda V994 Her: podíl apsidálního stáčení A způsobený B (a C) |
| `v994s-dyn-B` | R | šestihvězda V994 Her: podíl apsidálního stáčení B způsobený A (a C) |
| `v994s-dyn-C` | R | šestihvězda V994 Her: podíl apsidálního stáčení C způsobený A a B |
| `v994s-gr-A` | R | V994 Her: relativistický podíl apsidálního stáčení A |
| `v994s-gr-B` | R | V994 Her: relativistický podíl apsidálního stáčení B |
| `v994s-gr-C` | R | V994 Her: relativistický podíl apsidálního stáčení C |
| `qzcar-evolutionary-mass-ac1` | own/tab | QZ Car: současná hmotnost modelu osamocené hvězdy MIST s poloměrem a Teff složky Ac1 proti její dynamické hmotnosti |
| `qzcar-evolutionary-mass-ac2` | own/tab | totéž pro Ac2 |
| `qzcar-evolutionary-mass-aa1` | own/tab | totéž pro Aa1 |
| `agb-final-mass-1p5` | own/tab | hmotnost bílého trpaslíka (jádro na konci AGB) hvězdy MIST 1,5 M☉ (sluneční \[Fe/H\]) proti vztahu mezi počáteční a konečnou hmotností bílých trpaslíků |
| `agb-final-mass-2p5` | own/tab | hmotnost bílého trpaslíka (jádro na konci AGB) hvězdy MIST 2,5 M☉ (sluneční \[Fe/H\]) proti vztahu mezi počáteční a konečnou hmotností bílých trpaslíků |
| `agb-final-mass-3p0` | own/tab | hmotnost bílého trpaslíka (jádro na konci AGB) hvězdy MIST 3,0 M☉ (sluneční \[Fe/H\]) proti vztahu mezi počáteční a konečnou hmotností bílých trpaslíků |
| `agb-final-mass-4p0` | own/tab | hmotnost bílého trpaslíka (jádro na konci AGB) hvězdy MIST 4,0 M☉ (sluneční \[Fe/H\]) proti vztahu mezi počáteční a konečnou hmotností bílých trpaslíků |
| `murphy2014-kic9651065-a1` | own | první orbitální harmonická A1 časových zpoždění KIC 9651065 (delta Sct na excentrické dráze s periodou 272,7 d) z její dráhy (a1 sin i = 0,37 au, e = 0,468, varpi = 2,01) |
| `murphy2014-kic9651065-a2-a1` | own | poměr druhé a první orbitální harmonické časových zpoždění KIC 9651065 (určuje excentricitu; e = 0,468, varpi = 2,01 rad) |
| `delta-scuti-pl-1p8` | own/tab | M\_V hvězdy MIST 1,8 M☉ při Teff = 7500 K na hlavní posloupnosti minus M\_V z pozorovaného vztahu P-L hvězd delta Scuti v základním módu při její periodě z Q = 0,0329 d |
| `delta-scuti-pl-2p0` | own/tab | jako delta-scuti-pl-1p8 pro hvězdu 2,0 M☉ |
| `cmdra-a-radius-excess-solar` | own/tab | pozorovaný poloměr CM Dra A vůči poloměru MIST pro její hmotnost (4 Gyr, \[Fe/H\] = 0), minus 1 |
| `cmdra-a-radius-excess-lowz` | own/tab | totéž pro \[Fe/H\] = -0,3 (metalicita podle Terriena et al. 2012) |
| `czev343-apsidal-period` | R | CzeV343: pozorovaná apsidální perioda dvojhvězdy A, úplná N-částicová čtyřhvězda \[roky\] |
| `sse-sun-radius` | own/tab | Slunce: poloměr hvězdy SSE 1 M☉, Z = 0,02 ve stáří 4,57 Gyr |
| `sse-sun-luminosity` | own/tab | Slunce: svítivost hvězdy SSE 1 M☉, Z = 0,02 ve stáří 4,57 Gyr |
| `sse-sun-convection-zone` | own/tab | Slunce: spodní hranice konvektivní obálky r\_bcz/R (model obálky BSE) |
| `iw-q-solar` | own | log10 Q' inerciálních vln, hvězda hlavní posloupnosti 1 M☉, P\_rot = 10 d |
| `iw-q-m-dwarf` | own | log10 Q' inerciálních vln, plně konvektivní hvězda 0,3 M☉, P\_rot = 10 d |
| `tides-circularization-m67` | own | perioda zakulacení dvojhvězd slunečního typu ve stáří 4 Gyr (M67), slapy podle Hurleyho et al. (2002) (pozorováno: 12,1 +1,0/−1,5 d, Meibom & Mathieu 2005, tabulka 3) |
| `tides-circularization-m67-iw` | own | perioda zakulacení ve stáří 4 Gyr (M67) s inerciálními vlnami, jen hlavní posloupnost (jen rovnovážný slap: 3,8 d; pozorováno jako výše) |
| `gw-hulse-taylor-pdot` | own | zkracování dráhy dP/dt pulzaru PSR B1913+16 gravitačními vlnami (hmotnosti 1,438 a 1,390 M☉, P = 0,323 d, e = 0,617) |
| `epsaur-mid-eclipse` | own | epsilon Aur: střed zákrytu 2009–2011 ze spektroskopické dráhy |
| `epsaur-duration` | own | epsilon Aur: trvání zákrytu, kontakty I–IV |
| `epsaur-depth` | own | epsilon Aur: hloubka zákrytu ve V (neprůhledný disk, tloušťka 0,6 au) |
| `hd98800-polar-disc` | own | HD 98800 B: sklon polárního disku k rovině oblohy (podél vektoru excentricity dráhy B) |
| `eecep-disc-precession` | own | perioda precese \[roky\] disku EE Cep z momentu síly hvězdy B na tuhý disk (ilustrativní disk z configs/catalog/ee\_cep.yaml, r\_out 0,8 au) |
| `tides-circularization-pms-population` | own | perioda zakulacení (konečná perioda, e < 0,01) dvojhvězd dvou hvězd 1 M☉ ve stáří 3 Myr (MIST před hlavní posloupností od 1 Myr, rovnovážný slap s lineárním tlumením a inerciální vlny; pozorované dvojhvězdy před hlavní posloupností: 7,1 ± 1,2 d, Meibom & Mathieu 2005, tabulka 3) |
| `tides-circularization-m67-final` | own | perioda zakulacení (konečná perioda) dvojhvězd dvou hvězd 1 M☉ ve stáří 4 Gyr (MIST od 1 Myr, rovnovážný slap s lineárním tlumením a inerciální vlny; pozorováno v M67: 12,1 +1,0/−1,5 d) |

**Tabulka A8.** Úroveň D: výsledky (sloupce *Očekáváno*, *Spočteno*, *Rozdíl*, *Tolerance*, *Verdikt*) ve stejném pořadí jako v tabulce A7.

| Id | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|
| `v994s-dyn-A` | 0,12 | 0,0999365 | −16,72 % | 25 % | PASS |
| `v994s-dyn-B` | 0,05 | 0,0519361 | +3,87 % | 25 % | PASS |
| `v994s-dyn-C` | 0 | 1,23366·10⁻⁶ | +1,23·10⁻⁶ | 0,01 | PASS |
| `v994s-gr-A` | 0,1 | 0,0938886 | −6,11 % | 10 % | PASS |
| `v994s-gr-B` | 0,08 | 0,0748582 | −6,43 % | 10 % | PASS |
| `v994s-gr-C` | 0,06 | 0,0811201 | +35,20 % | 10 % | **EXPLAINED** |
| `qzcar-evolutionary-mass-ac1` | 25,5 M☉ | 27,0258 M☉ | +5,98 % | 15 % | PASS |
| `qzcar-evolutionary-mass-ac2` | 33,2 M☉ | 26,0243 M☉ | −21,61 % | 15 % | **EXPLAINED** |
| `qzcar-evolutionary-mass-aa1` | 69,8 M☉ | 42,5172 M☉ | −39,09 % | 15 % | **EXPLAINED** |
| `agb-final-mass-1p5` | 0,609 M☉ | 0,565426 M☉ | −0,0436 M☉ | 0,1 M☉ | PASS |
| `agb-final-mass-2p5` | 0,689 M☉ | 0,610396 M☉ | −0,0786 M☉ | 0,1 M☉ | PASS |
| `agb-final-mass-3p0` | 0,745 M☉ | 0,653775 M☉ | −0,0912 M☉ | 0,1 M☉ | PASS |
| `agb-final-mass-4p0` | 0,899 M☉ | 0,835342 M☉ | −0,0637 M☉ | 0,1 M☉ | PASS |
| `murphy2014-kic9651065-a1` | 165,57 s | 167,767 s | +1,33 % | 6 % | PASS |
| `murphy2014-kic9651065-a2-a1` | 0,231986 | 0,219651 | −0,0123 | 0,025 | PASS |
| `delta-scuti-pl-1p8` | 0 mag | −0,0251584 mag | −0,0252 mag | 0,25 mag | PASS |
| `delta-scuti-pl-2p0` | 0 mag | −0,115268 mag | −0,115 mag | 0,25 mag | PASS |
| `cmdra-a-radius-excess-solar` | 0,035 | 0,0335989 | −0,0014 | 0,015 | PASS |
| `cmdra-a-radius-excess-lowz` | 0,06 | 0,0641756 | +0,00418 | 0,015 | PASS |
| `czev343-apsidal-period` | 33,4 roku | 34,491 roku | +3,27 % | 6,0 % | PASS |
| `sse-sun-radius` | 1 R☉ | 0,984843 R☉ | −1,52 % | 5 % | PASS |
| `sse-sun-luminosity` | 1 L☉ | 0,95786 L☉ | −4,21 % | 5 % | PASS |
| `sse-sun-convection-zone` | 0,713 | 0,733903 | +0,0209 | 0,03 | PASS |
| `iw-q-solar` | 7 | 7,0122 | +0,0122 | 0,3 | PASS |
| `iw-q-m-dwarf` | 7 | 6,94697 | −0,053 | 0,3 | PASS |
| `tides-circularization-m67` | 12,1 d | 3,7847 d | −68,72 % | 30,0 % | **EXPLAINED** |
| `tides-circularization-m67-iw` | 12,1 d | 7,4572 d | −38,37 % | 30 % | **EXPLAINED** |
| `gw-hulse-taylor-pdot` | −2,40263·10⁻¹² | −2,4021·10⁻¹² | +0,02 % | 0,1 % | PASS |
| `epsaur-mid-eclipse` | 2,45543·10⁶ JD | 2,45541·10⁶ JD | −15 JD | 30 JD | PASS |
| `epsaur-duration` | 740 d | 638 d | −13,78 % | 15 % | PASS |
| `epsaur-depth` | 0,75 mag | 0,83134 mag | +0,0813 mag | 0,15 mag | PASS |
| `hd98800-polar-disc` | 26° | 28,4735° | +2,47° | 3° | PASS |
| `eecep-disc-precession` | 65 let | 183,163 roku | +181,79 % | 15 % | **EXPLAINED** |
| `tides-circularization-pms-population` | 7,1 d | 6,5831 d | −0,517 d | 1,2 d | PASS |
| `tides-circularization-m67-final` | 12,1 d | 7,5337 d | −4,57 d | 1,2 d | **EXPLAINED** |

\normalsize
