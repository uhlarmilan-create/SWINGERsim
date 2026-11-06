---
title: "SWINGERsim"
subtitle: "Uživatelská příručka"
---

# Úvod {#introduction}

SWINGERsim je program pro modelování **konkrétních pozorovaných vícenásobných hvězd**: těsné
dvojhvězdy se třetí hvězdou, dvou dvojhvězd obíhajících kolem sebe (čtyřhvězdy 2+2), trojhvězdy se
vzdálenou čtvrtou hvězdou (3+1) a větších hierarchií až zhruba do deseti hvězd, jako je trojnásobně
zákrytová šestihvězda V994 Her. Popíšete jednu soustavu – hmotnosti hvězd, jak jsou spárovány, po
jakých drahách se pohybují a jaká fyzika je podstatná – a SWINGERsim spočítá, jak se soustava vyvíjí
a co by naměřil pozorovatel: zákryty, které se objevují a mizí, okamžiky minim a jejich diagram O−C,
světelné křivky, radiální rychlosti a také to, zda soustava přežije, nebo zda spolu dvě její hvězdy
splynou.

Tato Uživatelská příručka vysvětluje, k čemu software slouží, jak do sebe zapadají jeho části, jak
ho nainstalovat a ovládat a co do každého postupu vstupuje a co z něj vychází. Fyziku neodvozuje;
to je předmětem doprovodného dokumentu *Fyzika a metody*.

## K čemu SWINGERsim slouží {#what-for}

**Hierarchická vícenásobná hvězda** je soustava, jejíž dráhy mají velmi rozdílné velikosti, takže
každá dráha je blízká keplerovské elipse, kterou ostatní dráhy pomalu ruší. Každá dráha takové
soustavy spojuje dvě složky, z nichž každá je hvězda nebo menší podsoustava; celá soustava je stromem
dvojic (obrázek 1). Právě pomalé poruchy činí tyto soustavy zajímavými pro pozorovatele: roviny drah
precedují, takže se mění sklon zákrytové dvojhvězdy k zornému paprsku a její zákryty se prohlubují,
mělčí nebo mizí; excentricity oscilují; okamžiky minim se posouvají, protože se dvojhvězda pohybuje po
své vnější dráze (efekt světelného času, LTTE) a protože ji přitahují ostatní hvězdy.

![**Obrázek 1.** Hierarchie, které SWINGERsim popisuje, nakreslené jako stromy. Oranžové kroužky jsou
hvězdy (pojmenované `m1`, `m2`, …), zelené čtverečky jsou dráhy (pojmenované `A`, `B`, `out`, …);
každá dráha spojuje dvě složky nakreslené pod ní. Panely ukazují dvojhvězdu, trojhvězdu (2+1), oba
druhy čtyřhvězd (2+2 a 3+1), tři pětihvězdy a šestihvězdy, které mají v konfiguraci připravené názvy
(2+2+1, 2+(2+1), 2+2+2, 3+3), a desetihvězdný strom zapsaný ručně (`custom`). Strom se čte zdola: ve
3+1 tvoří `m1` a `m2` dráhu `inner`, `inner` a `m3` tvoří `mid` a `mid` a `m4` tvoří nejvzdálenější
dráhu `out`.](../../figures/concepts/hierarchies.png)

Obrázek 1 je klíčem ke každému konfiguračnímu souboru a ke každé výstupní tabulce: názvy drah ve
stromu (`A`, `B`, `out`, `inner`, `mid`, …) jsou názvy, pod nimiž se elementy drah zadávají na vstupu
a uvádějí na výstupu. Předvolby v horních dvou řadách a obě šestihvězdy mají pevné názvy; vlastní
strom používá názvy, které napíšete.

SWINGERsim pracuje vždy s jednou soustavou a nabízí tři způsoby, jak se na ni ptát.

| Způsob práce | Na jakou otázku odpovídá | Typický výsledek |
| --- | --- | --- |
| dopředný model | „Co dělá soustava s těmito parametry a co bychom pozorovali?“ | vývoj drah za roky až miliony let, intervaly zákrytů, diagram O−C, světelná křivka |
| proložení naměřených dat | „Které parametry této soustavy reprodukují naměřená minima, radiální rychlosti, světelné křivky a sklony?“ | proložené elementy s chybami, rezidua, předpovědi budoucích minim, proložená konfigurace, kterou lze spustit jako kteroukoli jinou |
| ansámbl | „Některé parametry neznáme – co je možné a co dovolují měření?“ | rozdělení výsledků přes tisíce vylosovaných soustav, vážená měřeními |

**Tabulka 1.** Tři způsoby práce se SWINGERsimem. Každý řádek uvádí otázku a pravý sloupec to, co
dostanete. Všechny tři spolu souvisejí: proložení vytvoří konfiguraci, která se pak spustí dopředu,
a ansámbl je sada dopředných modelů, jejichž členy jsou váženy stejnými druhy měření, jaké používá
proložení.

## Na jaké otázky odpovídá {#questions}

Následující otázky jsou typické pro soustavy, pro něž je SWINGERsim určen. Na každou odpovídá jeden
z postupů z kapitoly [Pracovní postupy](#workflows).

- Kdy začaly zákryty této dvojhvězdy a kdy skončí?
- Vysvětluje pozorovanou změnu hloubky zákrytů precese vyvolaná známým třetím tělesem, nebo je
  potřeba další, neviděný průvodce?
- Které vnější dráhy a hmotnosti neviděného průvodce reprodukují naměřený sklon a rychlost jeho
  změny?
- Reprodukuje N-částicový model celé soustavy naměřené okamžiky minim všech zákrytových párů
  najednou, a s jakými elementy?
- Je soustava dynamicky stabilní, a jak dlouho? Je blízko rezonance nebo v chaotické oblasti?
- Jak soustava vypadá ve světelné křivce z TESS a jak přesně by se změřila minima?
- Jak se hvězdy a dráhy mění za miliony let: kdy hvězda vyplní svůj Rocheův lalok a co následuje
  (přenos hmoty, kontaktní fáze, společná obálka, splynutí)?
- Může být vlna v diagramu O−C projevem magnetické aktivity místo třetího tělesa?

## Co nedělá {#not-for}

SWINGERsim je záměrně úzce zaměřený. Není to kód pro populační syntézu: nelosuje celé galaxie
vícenásobných hvězd a nemodeluje výběrové efekty přehlídek (ansámbl může losovat ze širokých
rozdělení, ale interpretace jako populace je na uživateli). Nezpracovává fotometrii ani spektra
a neměří okamžiky minim na skutečných datech; vychází z tabulek naměřených veličin. Jeho světelné
křivky se počítají z poloh hvězd v dynamickém modelu (fotodynamický model) s kulovými, rotujícími nebo
Rocheovými hvězdami; nenahrazuje specializovaný kód pro světelné křivky tam, kde jde o podrobnou fyziku
povrchu jedné dvojhvězdy (skvrny, podrobné atmosféry). Supernovy a vývoj neutronových hvězd a černých
děr jsou mimo jeho rozsah a nehierarchické (demokratické) soustavy se objevují jen jako výsledek
rozptylových experimentů. Meze jednotlivých fyzikálních složek uvádí dokument *Fyzika a metody*,
kapitola Omezení.

## Soubor dokumentace {#documentation-set}

Dokumentaci SWINGERsimu tvoří sedm dokumentů.

| Dokument | Obsah | Čtěte, když |
| --- | --- | --- |
| *Uživatelská příručka* (tento dokument) | co software dělá, instalace, webová aplikace, příkazový řádek, pracovní postupy, konfigurace, výstupy | chcete software používat |
| *Fyzika a metody* | rovnice a algoritmy: hierarchie a dráhy, hvězdy, co vidí pozorovatel, N-částicová a sekulární dynamika, hvězdný vývoj, slapy, přenos hmoty, stabilita, pozorovatelné veličiny, prokládání, ansámbly, omezení | potřebujete vědět, co se počítá a s jakou přesností |
| *Testování, část I: Ověření proti literatuře* | srovnávací testy: publikované výsledky reprodukované SWINGERsimem | chcete vědět, nakolik je fyzika ověřená |
| *Testování, část II: Tři publikované vícenásobné hvězdy* | BU CMi, TIC 120362137 a V994 Her modelované z jejich dat | chcete vidět řešené aplikace na skutečné soustavy |
| *Testování, část III: Zákrytové dvojhvězdy s proměnným sklonem* | V907 Sco, QX Cas a AH Cep | studujete dvojhvězdy, jejichž zákryty přicházejí a odcházejí |
| *Testování, část IV: Zvláštní případy* | otázky bez publikované odpovědi: budoucnost CzeV343 (přenos hmoty, nestabilita, splynutí vlivem ZLK), hierarchie CzeV4315 (2+2 proti 2+(2+1) podle bayesovské evidence) a konvergence proložení MCMC | kombinujete vývoj, dynamiku a ansámbly nebo srovnáváte hypotézy |
| *Otevřené otázky a další práce* | co není vyřešeno a co by mohlo následovat | plánujete další práci |

**Tabulka 2.** Sedm dokumentů dokumentace SWINGERsimu. Prostřední sloupec říká, co který obsahuje,
pravý sloupec, kdy ho číst. Odkazy mezi dokumenty se v textu uvádějí názvem dokumentu a jeho
kapitoly.

## Jak číst tuto příručku {#how-to-read}

Příručka je psána pro astronomy a studenty astrofyziky; kromě úpravy textového souboru a napsání
příkazu není potřeba žádné programování. Každý odborný pojem je vysvětlen při prvním použití nebo ve
[Slovníku pojmů](#glossary).

| Chcete-li … | čtěte |
| --- | --- |
| vědět, jak je program postaven a která část co dělá | [Jak SWINGERsim funguje](#how-it-works) |
| program nainstalovat | [Instalace](#installation) |
| vidět výsledky za pár minut | [Rychlý start](#quick-start) |
| pracovat s formuláři ve webovém prohlížeči | [Webová aplikace](#gui) |
| pracovat z terminálu nebo ze skriptu | [Příkazový řádek](#command-line) |
| splnit konkrétní úkol (proložit minima, spustit ansámbl, sledovat vývoj) | [Pracovní postupy](#workflows) a [Řešené příklady](#examples) |
| vyhledat klíč konfiguračního souboru | [Konfigurační soubor](#configuration) |
| porozumět výstupnímu souboru nebo sloupci | [Výstupy](#outputs) |
| vyřešit problém | [Řešení problémů a časté otázky](#troubleshooting) |
| vyhledat pojem | [Slovník pojmů](#glossary) |

**Tabulka 3.** Kde co v této příručce najdete. Levý sloupec je úkol, pravý sloupec kapitola, která
ho pokrývá; po prvních dvou kapitolách lze ostatní číst v libovolném pořadí.

Příkazy jsou zapsány tak, jak se píší, např. `swingersim run config.yaml`; konfigurační klíče se
zapisují s tečkami pro vnoření, např. `run.t_max_yr` je klíč `t_max_yr` uvnitř oddílu `run`
a `orbits.<name>.P` je perioda libovolné dráhy. Jednotky odpovídají zvyklostem oboru: hmotnosti
v M☉, poloměry v R☉, periody ve dnech (d), velké poloosy v au, úhly ve stupních a časy v letech (yr),
není-li uvedeno jinak.
