---
title: "Testování SWINGERsimu"
subtitle: "Část I: Ověření proti literatuře"
---

```{=latex}
\makeatletter
\addtocontents{toc}{\protect\renewcommand*{\protect\numberline}{\protect\hb@xt@ 3.4em}}
\makeatother
```

**Abstrakt.** SWINGERsim simuluje hierarchické vícenásobné hvězdy – trojhvězdy, čtyřhvězdy a větší
vnořené soustavy – spolu s vývojem jejich hvězd, slapy, přenosem hmoty, zákryty, světelnými
křivkami, okamžiky zákrytů, prokládáním pozorování a ansámbly možných soustav. Tato první ze čtyř
zpráv o testování se ptá, zda každý fyzikální modul SWINGERsimu reprodukuje nezávislé publikované
výsledky. Důkazy mají čtyři úrovně: analytické výsledky (A), publikované numerické výsledky (B),
jiné programy spuštěné se stejným vstupem (C) a pozorované soustavy (D). SWINGERsim obsahuje 219
registrovaných srovnávacích testů (benchmarků), každý se svým zdrojem, očekávanou hodnotou
a tolerancí stanovenou před srovnáním. S konečnou verzí kódu proběhlo všech 219: 205 leží v rámci
své tolerance a 14 mimo ni z doloženého důvodu, který není chybou SWINGERsimu (přehlédnutí nebo
nekonzistence v článku, odlišná definice, délka publikovaných běhů nebo omezení společné všem
publikovaným modelům). Žádný test neselhává. Tři čtvrtiny testů (165 z 219) ověřují kód napsaný
přímo pro SWINGERsim – sekulární backendy, vazbu hvězdného vývoje, slapů a přenosu hmoty na dráhy,
světelné křivky na koulích a povrchových sítích, okamžiky zákrytů, prokládání a ansámbly – a nejen
externí N-částicový integrátor, který program používá. Kromě srovnávacích testů SWINGERsim
reprodukuje publikovaný model hmotné čtyřhvězdy QZ Car, dynamiku kandidátů skupiny SQUADRA CzeV4315
a ASASSN-V J075203, 23 publikovaných vícenásobných hvězd napříč Hertzsprungovým–Russellovým
diagramem a skóre stability publikovaného neuronového klasifikátoru. Zpráva končí problémy
zjištěnými během testování, popsanými jako vlastnosti programu, které musí uživatel znát,
a přehledem toho, co je ověřeno, s jakou přesností, a co ověřeno není.

# Úvod {#introduction}

## Cíl této zprávy {#aim}

Simulační program je astronomům užitečný jen tehdy, je-li jeho fyzika správná. To, že program
běží bez chyb nebo že souhlasí sám se sebou, dokazuje málo: chybná síla může dokonale zachovávat
energii. Cíl této zprávy je proto úzký a přísný. Pro každý fyzikální jev ve SWINGERsimu se ptáme:
**reprodukuje SWINGERsim nezávislý publikovaný výsledek pro tento jev, a s jakou přesností?** Kde
tomu tak není, ptáme se proč a rozdíl sledujeme až k jeho příčině.

Část I prochází fyzikální moduly jeden po druhém ([kapitoly 2–4](#dynamics)) a poté celé
publikované soustavy: hmotnou čtyřhvězdu QZ Car ([kapitola 5](#qz-car)), kandidáty skupiny SQUADRA
([kapitola 6](#squadra-systems)), katalog 23 publikovaných vícenásobných hvězd
([kapitola 7](#catalogue)), příklady, které sledují vícenásobnou hvězdu během jedné fáze jejího
života ([kapitola 8](#life-cycles)), a nástroje pro stabilitu ([kapitola 9](#stability)).
[Kapitola 10](#problems-found) uvádí problémy, které testy odhalily, a to, jak s nimi SWINGERsim
zachází, [kapitola 11](#conclusions) shrnuje, co je ověřeno a co ne, a [příloha A](#appendix-benchmarks)
uvádí všech 219 srovnávacích testů s jejich výsledky.

## Co je SWINGERsim a které jeho části jsou vlastním kódem {#own-code}

SWINGERsim sleduje hierarchické vícenásobné hvězdy: soustavy tří až asi deseti hvězd, v nichž
kolem sebe obíhají těsné páry a tyto páry obíhají kolem společných těžišť po širších drahách.
Dráhy integruje jedním ze tří **backendů** (výpočetních jader, která posouvají soustavu v čase):

- **N-částicový backend** přímo integruje pohybové rovnice všech hvězd. Používá integrátor IAS15
  knihovny REBOUND (Rein & Liu 2012; Rein & Spiegel 2015) a tři efekty knihovny REBOUNDx (Tamayo
  et al. 2020): `gr_full` (relativita prvního postnewtonovského řádu pro všechna tělesa),
  `tides_spin` (rovnovážný slap s rotacemi, Eggleton et al. 1998; Lu et al. 2023) a `modify_mass`
  (hmotnosti proměnné v čase);
- **sekulární backendy** středují pohyb přes oběžné periody a sledují jen pomalé změny drah
  (precesi, cykly ZLK, definované v [kapitole 2](#dynamics)). Jsou napsány pro SWINGERsim:
  kvadrupólově-oktupólový kód pro trojhvězdy, kód 2+2 podle Vokrouhlického (2016) a vnořený kód
  pro libovolnou hierarchii podle Hamerse et al. (2015) a Hamerse & Portegiese Zwarta (2016);
- **vývojový backend** spojuje sekulární backend s vyvíjejícími se hvězdami: běh rozdělí na úseky,
  v nichž se žádný poloměr hvězdy nezmění o více než 1 %, a mezi nimi uplatní hvězdný vítr, přenos
  hmoty, společné obálky, slapy a magnetické brzdění. Je napsán pro SWINGERsim.

Otázka „jaký kód srovnávací test ověřuje?“ je důležitá pro váhu důkazu. Test N-částicového
backendu, který reprodukuje například relativistické stáčení perihelia Merkuru, potvrzuje hlavně
to, že je správně použit REBOUNDx – samotný silový zákon pochází z REBOUNDx. Test sekulárního
backendu, vazby přenosu hmoty nebo modelu světelné křivky ověřuje rovnice a numeriku, které
existují jen ve SWINGERsimu; je silnějším důkazem, protože nic mimo SWINGERsim nemohlo způsobit,
že projde. Prošli jsme zdrojový kód (`src/swingersim/`) a každý srovnávací test zařadili do jedné
ze čtyř tříd (tabulka 1). Třída je uvedena u každého testu v tabulkách kapitol 2–9
a v [příloze A](#appendix-benchmarks).

**Tabulka 1.** Třídy srovnávacích testů podle kódu, který ověřují. Sloupec 1 je označení
používané ve všech dalších tabulkách; sloupec 2 říká, jaký kód počítá ověřované číslo; sloupec 3
uvádí hlavní moduly; poslední sloupec udává, kolik z 219 testů do třídy patří.

| Třída | Co počítá ověřované číslo | Moduly (příklady) | Testy |
|---|---|---|---|
| own | rovnice a numerika napsané pro SWINGERsim | sekulární backendy, vývojový backend, přenos hmoty a společná obálka, slapy a vzorce apsidálního stáčení, inerciální vlny, magnetické brzdění, Rocheova geometrie a povrchové sítě, světelné křivky, prachové disky, okamžiky zákrytů, časování pulzací, prokládání minim, ansámbly, kritéria stability, teorie úniku, gravitační vlny, slapy Galaxie a procházející hvězdy | 126 |
| own/tab | interpolace a použití publikovaných tabulek ve SWINGERsimu | vývojové dráhy MIST a SSE, strukturní konstanty Clareta (2004, 2023), okrajové ztemnění Clareta & Bloemenové (2011) a TLUSTY | 39 |
| R | REBOUND IAS15 s efekty REBOUNDx, řízený nastavením SWINGERsimu | hierarchická inicializace, převod elementů, detekce srážek a změn hierarchie, časový plán ztráty hmoty, konvence časového zpoždění, operátorové štěpení vlastních sil SWINGERsimu | 50 |
| ext | externí knihovna přes obal (wrapper) SWINGERsimu | hybridní světelné křivky PHOEBE 2.5, balík `sbi`, neuronová síť scikit-learn natrénovaná na vlastních trénovacích souborech SWINGERsimu | 4 |

Ani test třídy „R“ neověřuje jen REBOUND: SWINGERsim sestavuje počáteční polohy a rychlosti
z hierarchických dráhových elementů konfigurace, převádí zintegrovaný stav zpět na elementy,
detekuje srážky a změny hierarchie a nastavuje parametry efektů REBOUNDx (například faktor mezi
časovým zpožděním SWINGERsimu a parametrem `tau` efektu `tides_spin`). Test třídy R proto ověřuje
toto nastavení a externí integrátor společně. Kde test srovnává backend SWINGERsimu
s N-částicovým během téže soustavy (například sekulární backend proti N-částicovému), třída
označuje ověřovanou stranu, kterou je vlastní kód SWINGERsimu.

Z 219 srovnávacích testů jich 165 (75 %) ověřuje vlastní kód SWINGERsimu (třídy own a own/tab),
50 N-částicové nastavení kolem REBOUND a 4 externí knihovnu přes obal.

## Pravidlo ověření {#rule}

Každý fyzikální jev ve SWINGERsimu splňuje pět pravidel.

1. **Nezávislá reference.** Jev se ověřuje proti výsledku, který SWINGERsim nevytvořil: proti
   rovnici z literatury, obrázku nebo tabulce článku, jinému programu nebo pozorované soustavě.
   Vnitřní kontroly (zpětné převody souřadnic, zachování energie) jsou rovněž součástí
   automatických testů, nepočítají se však jako ověření.
2. **Zaznamenaný zdroj.** Každý srovnávací test uchovává článek a rovnici, obrázek nebo tabulku,
   kterou reprodukuje, z nich převzaté vstupní hodnoty, očekávanou hodnotu a toleranci. Příkaz
   `swingersim validate` vypíše zdroj vedle každého výsledku.
3. **Fyzikální tolerance.** Tolerance vyplývá z přiblížení, na kterém reference stojí (například
   z vynechaného řádu rozvoje, přesnosti aproximačního vzorce, chyby měření nebo statistické chyby
   ansámblu). Stanoví se před srovnáním a nikdy se nerozšiřuje, aby test prošel.
4. **Každý rozdíl je vysvětlen, jinak je chybou.** Výsledek mimo toleranci je buď chybou
   SWINGERsimu, nebo má doložený důvod. Důvod je uložen u testu a vypíše se s výsledkem.
5. **Výsledky jsou reprodukovatelné.** Každý srovnávací test se spouští jedním příkazem a jeho
   číslo v této zprávě je číslo, které tento příkaz vypíše s konečnou verzí kódu.

## Úrovně důkazu a verdikty {#levels}

Reference spadají do čtyř **úrovní důkazu** (tabulka 2). Úrovně odpovídají na různé otázky.
Úroveň A ověřuje fyziku v limitách, v nichž je přesná odpověď známa. Úrovně B a C ověřují režimy,
pro které vzorec neexistuje – chaotickou dynamiku, provázaný vývoj hvězd a drah. Úroveň D ověřuje,
zda model jako celek popisuje skutečné hvězdy, včetně všeho, co model vynechává.

**Tabulka 2.** Úrovně důkazu. Sloupec 2 říká, čím je reference, sloupec 3 uvádí typické
příklady, sloupec 4 počet srovnávacích testů dané úrovně.

| Úroveň | Reference | Příklady | Testy |
|---|---|---|---|
| A | analytický výsledek z literatury | maximální excentricita ZLK, relativistické stáčení perihelia, Sterneův vzorec apsidálního stáčení, poloměr Rocheova laloku, uzavřené tvary přenosu hmoty | 93 |
| B | publikovaný numerický výsledek (obrázek nebo tabulka článku) | překlápění drah podle Naoz et al. (2013), vzájemné sklony podle Vokrouhlického (2016), účinné průřezy rozptylu podle Antogniniho & Thompsona (2016) | 60 |
| C | jiný program se stejným vstupem | MSE, BSE v COSMIC, PHOEBE 2.5, batman, ellc, model Xitau pro QZ Car, ženevské modely cefeid | 31 |
| D | pozorovaná soustava | V994 Her, CzeV343, QZ Car, Slunce, M67, bílí trpaslíci, PSR B1913+16, ε Aur, HD 98800 | 35 |

Každý srovnávací test končí jedním ze čtyř **verdiktů**. Ponecháváme jim anglická označení, pod
nimiž je vypisuje příkaz `swingersim validate`:

- **PASS** (vyhovuje) – výsledek leží v rámci tolerance;
- **EXPLAINED** (vysvětleno) – výsledek leží mimo toleranci z doloženého důvodu, který není chybou
  SWINGERsimu: vstup nebo definice použitá v článku, publikovaná čísla, která si navzájem odporují,
  délka nebo vzorkování publikovaných běhů, nebo uvedené omezení fyzikálního modelu SWINGERsimu,
  které je otevřenou otázkou, a nikoli chybou v kódu. Tolerance se nemění, aby rozdíl zakryla;
- **FAIL** (selhává) – mimo toleranci bez takového důvodu; šlo by o chybu SWINGERsimu;
- **SKIP** (přeskočeno) – není nainstalován volitelný balík, který test potřebuje.

![**Obrázek 1.** 219 srovnávacích testů SWINGERsimu podle úrovně důkazu (řádky A–D, viz tabulka 2) s konečnou verzí kódu. Délka sloupce je počet testů dané úrovně; modrá část počítá verdikty PASS, oranžová EXPLAINED. Čísla na konci každého sloupce opakují počty. Nevyskytuje se žádný FAIL ani SKIP.](../../figures/part1/benchmark_tiers.png)

Obrázek 1 shrnuje výsledek. Všech 93 analytických testů (úroveň A) vyhovuje. 14 vysvětlených
rozdílů leží v úrovních B (5), C (2) a D (7) – tedy tam, kde je referencí publikovaná simulace,
jiný program nebo skutečná hvězda. Právě to se očekává, je-li fyzika SWINGERsimu správná:
u analytické limity není žádný vstup, o který by se dalo přít, kdežto článek může uvádět číslo podle jiné
definice a skutečná hvězda obsahuje fyziku, kterou nemá žádný model. Každý ze 14 případů je
rozebrán tam, kde se probírá jeho fyzika, a souhrnně jsou uvedeny v [kapitole 11](#conclusions).

## Odkud čísla pocházejí {#numbers}

Srovnávací testy jsou součástí programu. `swingersim validate` spouští rychlé testy (každý
sekundy), `swingersim validate --all` přidává pomalé (každý minuty, například ansámbly stovek
rozptylových experimentů) a `swingersim validate <id>` spustí jeden test podle jména. Devět
*ručních* testů běží jen tehdy, jsou-li uvedeny jménem: potřebují hodiny výpočtu (N-částicový
ansámbl Pejchy et al. 2013), trénovací soubory 10⁴ N-částicových běhů (klasifikátor stability
Vynatheyi et al. 2023) nebo volitelný balík (PHOEBE, `sbi`).

```bash
swingersim validate                  # 184 rychlých testů
swingersim validate --all            # všech 210 automatických testů
swingersim validate gr-mercury       # jeden test podle jména
swingersim validate pejcha13-F6-fA   # ruční test (hodiny), běží jen při uvedení jména
```

Každé číslo v této zprávě bylo spočteno s konečnou verzí kódu SWINGERsimu: rychlé testy
v posledním běhu sady, pomalé testy v úplném běhu na produkčním stroji (12jádrová pracovní
stanice) a ruční testy v samostatných bězích téhož kódu. Příkaz zapíše také zprávu HTML
s tabulkou všech výsledků a s důvody verdiktů EXPLAINED.

Studie celých soustav ([kapitoly 5–9](#qz-car)) jsou simulacemi s konfiguračními soubory
distribuce SWINGERsimu; jejich čísla pocházejí z výstupů těchto běhů. Je-li studie ansámblem,
uvádí se s výsledkem i jeho velikost.

## Jak číst tuto zprávu {#how-to-read}

Každý oddíl kapitol 2–4 má stejné pořadí: fyzikální jev a jeho reference, vstup, výsledek
(tabulka, často obrázek) a shoda vyjádřená číslem. Tabulky srovnávacích testů mají všude stejné
sloupce:

- **Test** – identifikátor, pod nímž jej spouští `swingersim validate`;
- **Třída** – jaký kód ověřuje (tabulka 1);
- **Očekáváno** – referenční hodnota (analytická, publikovaná nebo z jiného programu);
- **Spočteno** – hodnota, kterou spočte SWINGERsim;
- **Rozdíl** – SWINGERsim minus očekávaná hodnota, buď relativně (v %), nebo absolutně
  (v jednotce veličiny), podle definice testu;
- **Tolerance** – největší rozdíl, který se počítá jako shoda, ve stejném tvaru;
- **Verdikt** – PASS nebo EXPLAINED.

Pojmy jsou definovány při prvním použití. Značení odpovídá ostatním dokumentům SWINGERsimu: $P$ je
perioda, $a$ velká poloosa, $e$ excentricita, $i$ sklon k rovině oblohy, $i_\mathrm{mut}$ vzájemný
sklon dvou drah, $\omega$ argument periastra a $\Omega$ délka uzlu. Dráhy čtyřhvězdy 2+2 se
nazývají A a B (oba páry) a „out“ (jejich vzájemná dráha); u čtyřhvězdy 3+1 jsou to „inner“,
„mid“ a „out“. Rovnice každého modulu jsou odvozeny v dokumentu *SWINGERsim – Fyzika a metody*;
tato zpráva uvádí jen to, co je potřeba k posouzení srovnání.

## Části II–IV {#parts}

Část I se ptá, zda fyzika SWINGERsimu reprodukuje literaturu. Ostatní tři části se ptají, zda
SWINGERsim, používaný tak, jak by jej používal astronom, reprodukuje skutečná pozorování soustav
studovaných skupinou SQUADRA:

- **Část II: Tři publikované vícenásobné hvězdy** – čtyřhvězda 2+2 BU CMi, čtyřhvězda 3+1
  TIC 120362137 a šestihvězda V994 Her: jejich publikovaná řešení integrovaná dopředu i dozadu,
  syntetické světelné křivky a okamžiky zákrytů proti datům z TESS a ze země a společné proložení
  všech dat;
- **Část III: Zákrytové dvojhvězdy s proměnným sklonem** – V907 Sco, QX Cas a AH Cep, jejichž
  zákryty se objevují, prohlubují nebo mizí, protože třetí těleso naklání jejich dráhy;
- **Část IV: Zvláštní případy** – čtyřhvězdy 2+2 CzeV343 a CzeV4315, které testují meze modelů
  (rezonance, nestabilní konfigurace, chaotické členy ansámblů).
