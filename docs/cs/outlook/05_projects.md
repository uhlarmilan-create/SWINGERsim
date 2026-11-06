# Výzkumné a studentské projekty {#projects}

SWINGERsim je určen k modelování konkrétních soustav z jejich dat. Níže uvedené projekty využívají
to, co program již umí. U každého je uvedena otázka, nástroje, které by použil, a data, která
potřebuje.

## Projekty s daty skupiny SQUADRA {#projects-squadra}

**Tabulka 3.** Projekty na soustavách pozorovaných skupinou SQUADRA. Sloupce: projekt a jeho
otázka, co poskytuje SWINGERsim a co ještě chybí (obvykle data).

| projekt | otázka | SWINGERsim poskytuje | ještě chybí |
|:--|:--|:--|:--|
| A. Soustavy s proměnným sklonem | které vnější dráhy vysvětlují naměřené změny sklonu (CzeV4315, ASASSN-V J075203 a podobné) a kdy se jejich zákryty objeví nebo zmizí? | ansámbly s vahami pozorování a evidencemi pro alternativní hierarchie; dynamický fit minim, radiálních rychlostí a světelných křivek | okamžiky minim obou párů; sklony ze světelných křivek v několika epochách |
| B. Předpovědi zákrytů pro seznam SQUADRA | pro každý cíl pravděpodobnost, že zákryty kteréhokoli z párů zmizí během 10–50 let: které soustavy pozorovat nejdříve | konfigurace z publikovaných řešení, ansámbly, worker pro dávkové výpočty | úplná tabulka period, hmotností a vnějších drah |
| C. CzeV343 a rezonance 3:2 | je soustava v rezonanci, nebo jí prochází, a jak se změní její O−C? | apsidální stáčení, mapy chaosu (MEGNO), dlouhé integrace, řetězec fitu | diagnostika rezonančních úhlů (tabulka 2) |
| D. Kompaktní soustavy 2+2 s vnější periodou pod jeden rok | sekulární teorie selhává; co dává úplný N-částicový model jejich O−C? | dynamický fit (fáze 2) se všemi páry | minima cílů skupiny |

## Populační a teoretické otázky {#projects-theory}

**Tabulka 4.** Širší otázky. Sloupce: projekt, otázka, co poskytuje SWINGERsim.

| projekt | otázka | SWINGERsim poskytuje |
|:--|:--|:--|
| E. Orientace drah 2+2 | jsou dráhy soustav 2+2 většinou srovnané do jedné roviny (vznik v disku), nebo izotropní (dynamický záchyt), soudě podle statistiky dvojitě zákrytových soustav? | pravděpodobnosti zákrytů populací (ověřené proti Vokrouhlickému 2016); chybí: model výběrových efektů přehlídky |
| F. Přebytek poměru period 3:2 | proč soustavy 2+2 vykazují přebytek blízko 3:2? Z 781 kandidátů leží 16 do 1 % od 3:2 proti 11 očekávaným pro rovnoměrné periody (27 proti 21 do 2 %; Kolář et al. 2025) a 8 proti 12 ve vzdálenosti 2–3 % (Kolář 2025, disertace): přebytek v rezonanci a nedostatek vedle ní, jak by vytvořil záchyt, ale v malých počtech. Záchyt do rezonance během migrace drah, jak je znám u planet? Dnešní CzeV343 v rezonanci není (část IV). | mapy chaosu, parametrická plynná fáze, slapy |
| G. Splynutí a kontaktní dvojhvězdy z ZLK | kolik soustav 2+2 skončí jako splynulé hvězdy, kontaktní dvojhvězdy nebo dvojhvězdy s přenosem hmoty, a kdy? | ZLK se slapy, vývoj hvězd, přenos hmoty, společná obálka, kontaktní dvojhvězdy; část IV je jedním příkladem |
| H. Nestabilita po přenosu hmoty | jak často přenos hmoty destabilizuje čtyřhvězdu a co následuje: splynutí působením ZLK jako u CzeV343, nebo vymrštění? | dvoufázová metoda z části IV |
| I. Vznik rozptylem dvojhvězda–dvojhvězda | účinné průřezy pro vznik stabilních soustav 2+2 a 3+1 v mladých hvězdokupách | rozptylové experimenty, teorie úniku, slapový záchyt |
| J. Plynná fáze a orientace drah | jaké rozdělení vzájemných sklonů po sobě plynná fáze zanechá? | parametrická plynná fáze |
| K. Pulzující hvězdy jako hodiny | dávají zpoždění světelného času pulzací hvězd typu δ Scuti, β Cephei nebo γ Dor spolu s okamžiky zákrytů hmotnosti vnější dráhy bez spektroskopie? | pulzace vyhodnocené v čase vyzáření; minima ze syntetických světelných křivek |
| L. Přenos hmoty v trojhvězdách | mění třetí hvězda výsledek přenosu hmoty (algolidy, obnažené hvězdy) a mění přenos vnější dráhu a zákryty? | přenos hmoty ve vnitřních dvojhvězdách i z třetí hvězdy, s jeho vlivem na ostatní dráhy |
| M. Společná obálka ve vícenásobných soustavách | jak často společná obálka uvolní nebo destabilizuje vnější dráhu? | energetický předpis s odezvou vnějších drah |
| N. Mladé hmotné vícenásobné soustavy | která hvězda mladé čtyřhvězdy s hvězdami typu O začne interagovat první (soustavy podobné QZ Car)? | MIST až do vysokých hmotností, hvězdné větry, efektivní stáří (část I, QZ Car) |
| O. Skryté vícenásobné soustavy | kolik zákrytových dvojhvězd v TESS, ZTF a ASAS-SN začne nebo přestane zakrývat během časové základny přehlídky? | statistika objevování zákrytů v populaci |

## Studentské projekty {#projects-students}

**Tabulka 5.** Studentské projekty. Sloupce: úroveň, téma (s písmenem projektu z tabulek 3–4),
úkol, obtížnost.

| úroveň | téma | úkol | obtížnost |
|:--|:--|:--|:--|
| bakalářská | předpovědi zákrytů pro deset soustav SQUADRA (B) | konfigurace z publikovaných řešení, ansámbly, shrnutí, kdy zákryty zmizí | nízká |
| bakalářská | rozptyl dvojhvězda–dvojhvězda (I) | prohledání rychlostí a poměrů drah, srovnání s Antognini & Thompson (2016) | nízká až střední |
| bakalářská | časová zpoždění hvězdy typu δ Scuti v zákrytové dvojhvězdě z TESS (K) | změření zpoždění v datech TESS, srovnání s modelem SWINGERsimu | nízká až střední |
| magisterská | dynamický model CzeV4315 (A) | fit minim obou párů modelem 2+(2+1); článek společně s pozorovateli | střední |
| magisterská | QX Cas: poloměry, třetí světlo a staré zákryty (otevřená otázka z oddílu 3) | společný fit hloubek, „heartbeat“ z TESS a okamžiků minim; desky DASCH | střední |
| magisterská | mapy chaosu a rezonance 3:2 v CzeV343 (C, F) | mapy MEGNO, srovnání s rezonanční sekulární teorií Breitera & Vokrouhlického (2018) | střední až vysoká |
| magisterská | algolidy s třetími hvězdami (L) | ansámbly dvojhvězd s přenosem hmoty a třetími hvězdami, srovnání se známými soustavami | střední |
| doktorská | orientace drah 2+2 ze zákrytové populace (E) | model výběrových efektů přehlídky, bayesovská inverze | vysoká |
| doktorská | původ přebytku 3:2 (F) | N-částicové výpočty s migrací a slapy, srovnání s pozorovanými soustavami | vysoká |
| doktorská | osudy soustav 2+2 (G, H) | ZLK s vývojem hvězd a přenosem hmoty pro populaci, srovnání s kódy populační syntézy | vysoká |

Ansámbly o $10^4$ členech proběhnou na dvanáctijádrovém počítači během hodin, takže se student může
soustředit na fyziku.
