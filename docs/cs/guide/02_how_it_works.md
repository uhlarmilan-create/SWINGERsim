# Jak SWINGERsim funguje {#how-it-works}

Tato kapitola popisuje, co se děje mezi konfiguračním souborem a hotovými výsledky: zpracovatelský
řetězec, části programu, které provádějí jednotlivé kroky, jak jsou propojeny, a dva způsoby výpočtu
pohybu (N-částicový a sekulární backend). Končí mapou všech modulů a tabulkou, která odděluje to, co
SWINGERsim počítá vlastním kódem, od toho, co přebírá z externích knihoven a publikovaných dat.

## Zpracovatelský řetězec {#chain}

Každý výpočet prochází stejným řetězcem (obrázek 2). **Konfigurace** – textový soubor ve formátu
YAML – popisuje jednu soustavu. SWINGERsim ji zkontroluje, sestaví počáteční stav soustavy, integruje
pohyb jedním ze dvou **backendů** (výpočetních jader), analyzuje výsledek sadou **diagnostik** a zapíše
**produkty** (tabulky, obrázky, zprávu HTML) a **pozorovatelné veličiny** (okamžiky zákrytů, světelné
křivky, radiální rychlosti). Ansámbly, proložení a mapy jsou vrstvy nad tímto řetězcem: volají ho
mnohokrát s různými parametry. Příkazový řádek a webová aplikace výpočty spouštějí a zobrazují jejich
výsledky.

![**Obrázek 2.** Zpracovatelský řetězec SWINGERsimu. Bloky 1–5 (nahoře) jsou cesta jedné
konfigurace: vstupní soubor (šedý) se zkontroluje a převede na počáteční stav (modrý, blok 2),
integruje se N-částicovým backendem (3a) nebo sekulárním backendem (3b), analyzuje se (fialový,
blok 4) a převede se na produkty běhu a pozorovatelné veličiny (oranžové, 5a a 5b). Zelené bloky dole
(6a–6c) jsou metody, které tento řetězec spouštějí mnohokrát – ansámbly, proložení naměřených dat,
mapy a experimenty; šipky z nich míří do řetězce, který volají. Žlutý blok (7) je řídicí vrstva:
příkazový řádek, webová aplikace a zpracovatel úloh, který spouští soubory úloh. Každý blok uvádí
hlavní úkoly daného kroku; moduly, které je provádějí, jsou v tabulce
5.](../../figures/guide/pipeline.png)

Obrázek 2 se pro jeden běh čte zleva doprava. Šipky z bloku 2 se rozdělují, protože každý běh používá
právě jeden backend, zvolený v konfiguraci; oba backendy předávají diagnostikám stejný druh výsledku,
a proto vše napravo od nich funguje s kterýmkoli z nich. Spodní řada ukazuje, že ansámbly, proložení
a mapy nejsou samostatné programy: každý z nich vytváří mnoho konfigurací (nebo mnoho variant jedné)
a posílá je bloky 2–5.

Jednotlivé kroky podrobně:

1. **Vstup.** Konfigurace určuje hierarchii (`topology`), hmotnosti, keplerovské elementy každé
   dráhy na počátku ($t = 0$), volitelné vlastnosti hvězd, zahrnutou fyziku a délku běhu. Volitelné
   oddíly popisují ansámbl (`ensemble`), proložení naměřených dat (`fit`) a zakrývající prachové
   disky (`discs`). Naměřená data (okamžiky minim, radiální rychlosti, světelné křivky, sklony) jsou
   tabulky CSV, na které konfigurace odkazuje. Úplný přehled je v kapitole
   [Konfigurační soubor](#configuration).
2. **Kontrola a sestavení.** Soubor se před jakýmkoli výpočtem ověří podle schématu: každý klíč musí
   být známý, každá hodnota musí mít správný typ a rozsah, dráhy musí odpovídat hierarchii
   a nesmyslné kombinace jsou odmítnuty se zprávou, která jmenuje příslušný klíč. Z hierarchie se
   stane binární strom; elementy každé dráhy se převedou na polohy a rychlosti v Jacobiho
   souřadnicích (každá dráha je relativním pohybem těžišť jejích dvou složek). Každá hvězda dostane
   poloměr, efektivní teplotu, apsidální konstantu $k_2$, moment setrvačnosti a na požádání i slapové
   časové zpoždění, a to z modelu hvězdy ve stáří soustavy nebo z hodnot zadaných v souboru. Je-li
   zapnut přenos hmoty, spočítá se v tomto okamžiku z vývojových drah hvězd historie každé
   interagující dvojhvězdy (viz *Fyzika a metody*, kapitola Přenos hmoty, společná obálka
   a splynutí).
3. **Backend.** N-částicový backend integruje polohy a rychlosti všech hvězd; sekulární backend
   integruje elementy středované přes oběh. Oba jsou popsány [níže](#backends). Vyvíjejí-li se
   hvězdy, rozdělí se běh na úseky, v nichž jsou hvězdy pevné, a hvězdy se aktualizují mezi úseky
   (operátorové štěpení).
4. **Diagnostiky.** Časové řady elementů se převedou na čísla, která odpovídají na otázky: které
   hvězdy jsou v každém výstupu vázány na které (a zda se hierarchie změnila), rozsahy excentricity
   a vzájemného sklonu ve srovnání s analytickou teorií cyklů von Zeipela–Lidova–Kozaie (ZLK),
   periody precese a apsidální stáčení, podíl času, kdy každá dvojhvězda zakrývá, největší zaplnění
   Rocheova laloku a poměry stability každé vnořené trojice.
5. **Produkty a pozorovatelné veličiny.** Běh zapíše svůj souhrn, soubor HDF5 se všemi časovými
   řadami, tabulku CSV, obrázky PNG, interaktivní 3D prohlížeč a zprávu HTML. Pozorovatelné veličiny
   se na požádání počítají z téže dynamiky: okamžiky zákrytů a diagram O−C, syntetické světelné
   křivky, okamžiky minim na nich změřené, radiální rychlosti a periody pulzací.
6. **Vrstvy.** Ansámbly losují neznámé parametry a spouštějí každého člena; proložení upravují volné
   parametry, dokud model nereprodukuje naměřená data; mapy chaosu a rozptylové experimenty spouštějí
   mřížky nebo náhodné sady konfigurací.
7. **Řízení.** Příkaz `swingersim` s jedním podpříkazem pro každý úkol, webová aplikace
   (`swingersim gui`) a zpracovatel úloh (worker), který spouští soubory úloh na počítači, jenž pro
   dlouhé výpočty zůstává zapnutý.

## Dva backendy {#backends}

Pohyb hierarchické soustavy má dvě velmi odlišné časové škály: oběhy po drahách (dny až roky) a jejich
pomalé sekulární změny (precese, cykly ZLK; desetiletí až miliony let). Proto SWINGERsim nabízí dva
backendy.

| | N-částicový (`run.backend: nbody`) | sekulární (`run.backend: secular`) |
| --- | --- | --- |
| co se integruje | polohy a rychlosti všech hvězd | elementy středované přes oběh (vektory momentu hybnosti a excentricity každé dráhy) |
| integrátor | IAS15 z REBOUND, adaptivní, chyba na úrovni strojové přesnosti | DOP853 ze SciPy; rovnice 2+2 přeložené pomocí Numba, je-li nainstalována |
| relativita | 1PN pro všechny dvojice těles (REBOUNDx `gr_full`) | stáčení periastra každé dráhy |
| hvězdy | slapy a rotace jako síly na rozlehlá tělesa (REBOUNDx `tides_spin`), slapová disipace, změna hmotnosti | deformace jako apsidální a uzlová precese, slapová disipace, vyvíjející se hvězdy, ztráta hmoty a přenos hmoty |
| události | srážky, změny hierarchie, úniky, přetečení Rocheova laloku, mez chyby energie, mez skutečného času výpočtu | přetečení Rocheova laloku a události přenosu hmoty; hierarchie se změnit nemůže |
| výstupy | oskulační elementy | střední elementy |
| pozorovatelné veličiny (okamžiky zákrytů, světelné křivky, radiální rychlosti) a proložení | ano | ne |
| náročnost | sekundy až hodiny, podle počtu oběhů nejkratší dráhy | asi sekunda na běh, nezávisle na vnitřních periodách |
| použijte pro | desetiletí až tisíce let; mírné hierarchie; těsná přiblížení; nestabilitu; vše, co měří pozorovatel | miliony let; ansámbly s tisíci členy; silně hierarchické soustavy |

**Tabulka 4.** Dva backendy. Každý řádek srovnává jednu vlastnost; poslední řádek říká, kdy použít
který. „Oskulační“ elementy popisují okamžitou keplerovskou elipsu spočtenou z poloh a rychlostí;
„střední“ elementy jsou středované přes oběžný pohyb a chybějí jim krátkoperiodické zákmity.

Oba backendy vracejí stejnou strukturu – časové řady elementů každé dráhy, vzájemné sklony, chyby
energie a momentu hybnosti, události – takže každá diagnostika, výstupní soubor i stránka výsledků
funguje s kterýmkoli z nich a výsledek lze vždy zkontrolovat během s druhým backendem. Obrázek 3
ukazuje rozdíl na jedné soustavě.

![**Obrázek 3.** Tatáž hierarchická soustava spočtená oběma backendy. Vodorovná osa: čas v letech;
svislá osa: excentricita vnitřní dvojhvězdy. Tenká šedá čára je N-částicový backend, silná modrá čára
sekulární backend. Obě ukazují dva cykly ZLK, v nichž excentricita roste z 0,05 na 0,77 a vrací se.
N-částicová křivka nese drobné rychlé zákmity – oskulační elementy reagují na každý oběh vnější dráhy
–, které středované rovnice vynechávají; pomalý cyklus, podstatný pro dlouhodobý vývoj, je v obou
stejný.](../../figures/concepts/nbody_vs_secular.png)

Obrázek 3 ilustruje praktické pravidlo: jsou-li poměry period velké a soustava zůstává hierarchická,
středované rovnice věrně sledují pomalý vývoj za nepatrný zlomek výpočetní ceny. Přesnost ztrácejí,
je-li hierarchie mírná (vnější perioda jen několik desítek vnitřních period), blízko rezonancí
středních pohybů a v chaotických soustavách; souhrn běhu pak upozorní, že „sekulární odhady“ jsou
přibližné (parametr $\varepsilon_\mathrm{SA}$, viz [Výstupy](#outputs)). Sekulární rovnice čtyřhvězdy
2+2 vycházejí z práce Vokrouhlického (2016); všechny ostatní hierarchie používají obecnou formulaci pro
vnořené dvojhvězdy (Hamers et al. 2015), jejímž zvláštním případem jsou rovnice 2+2. Rovnice,
jejich platnost a testy proti N-částicovému backendu jsou v dokumentu *Fyzika a metody*, kapitola
Dynamika: N-částicová integrace, relativita a sekulární backend.

**Který backend zvolit.** Pro vše, co měří pozorovatel, a pro běhy do několika tisíc let začněte
s N-částicovým backendem. Sekulární backend použijte pro vývoj za miliony let, pro dlouhé ansámbly
a tehdy, když jsou vnitřní dvojhvězdy tak krátkoperiodické, že by N-částicová integrace potřebovala
miliardy kroků. V případě pochybností spusťte oba na krátkou dobu a porovnejte souhrny.

## Jak jsou části propojeny {#connections}

Části spolu komunikují jen prostřednictvím souborů a jedné vnitřní struktury výsledku, což je
udržuje nezávislé:

- **Konfigurace → vše.** Každý příkaz čte stejný soubor YAML a ověřuje ho stejným schématem.
  Konfigurace napsaná ve webové aplikaci, ručně nebo proložením (`stage2_fitted.yaml`) je soubor
  téhož druhu.
- **Backend → diagnostiky → soubory.** Oba backendy vracejí stejnou strukturu výsledku; diagnostiky
  čtou jen tuto strukturu; soubor HDF5 ukládá výsledek spolu s úplnou konfigurací, takže soubor běhu
  vždy říká, jak vznikl.
- **Pozorovatelné veličiny z dynamiky.** Okamžiky zákrytů, světelné křivky, radiální rychlosti
  i modelové pozorovatelné veličiny proložení integrují tytéž N-částicové rovnice jako běh; světelná
  křivka je proto v souladu s diagramem O−C téže konfigurace.
- **Vrstvy volají řetězec.** Ansámbl zapíše jeden souhrnný řádek na člena a směr času; proložení
  volá N-částicový model jednou pro každé vyhodnocení svých volných parametrů; mapa chaosu spouští
  jednu integraci v každém bodě mřížky.
- **Úlohy.** Webová aplikace sama nic nepočítá: zapisuje soubory úloh (malé soubory YAML, které
  uvádějí příkaz, konfiguraci a volby) do složky a zpracovatel úloh je provádí týmž příkazovým
  řádkem, jaký by napsal uživatel. Výsledky se načítají zpět z výstupních souborů.

## Mapa modulů {#module-map}

Program je balík Pythonu (`swingersim`, zdrojový kód v `src/swingersim/`). Tabulka 5 uvádí všechny
moduly seskupené podle role, s jedním řádkem o tom, co dělají, a zda jde o vlastní kód SWINGERsimu,
nebo o obal (wrapper) externí knihovny. „Vlastní“ znamená, že metoda je implementována ve SWINGERsimu
(pro pole, hledání kořenů a integrátory může používat NumPy a SciPy); „obal“ znamená, že podstatný
výpočet provádí uvedený externí balík a modul připravuje jeho vstup a čte jeho výstup; „vlastní
+ externí“ znamená, že vlastní algoritmus SWINGERsimu používá externí balík pro jednu přesně
vymezenou část.

| Modul | Co dělá | Vlastní kód nebo obal |
| --- | --- | --- |
| **Konfigurace a sestavení** | | |
| `config.py` | schéma konfigurace, jednotky, předvolby hierarchií, všechny kontroly konzistence | vlastní (schéma napsané v pydantic) |
| `units.py` | vnitřní jednotky (M☉, au, yr) a převody | vlastní |
| `kepler.py` | Keplerova rovnice, elementy ↔ poloha a rychlost | vlastní |
| `hierarchy.py` | strom hierarchie, inicializace v Jacobiho souřadnicích, elementy každého uzlu | vlastní |
| `tables.py`, `rootpath.py` | jednotné pravidlo pro čtení tabulek CSV a pro cesty úloh | vlastní |
| **Hvězdy** | | |
| `stellar.py` | poloměr, $T_\mathrm{eff}$, $k_2$, moment setrvačnosti, slapové zpoždění každé hvězdy; vztahy pro ZAMS, polytropy, Claretovy tabulky, uživatelské tabulky, poloměr Rocheova laloku | vlastní; data Tout et al. (1996), Claret (2004, 2023) |
| `evolution.py` | vývojové dráhy SSE (poloměr, zářivý výkon, obálka) pro libovolnou hmotnost, stáří a metalicitu interpolací | vlastní interpolace; mřížka spočtená publikovanými rutinami SSE (Hurley et al. 2000) |
| `mist.py` | vývojové dráhy MIST od fáze před hlavní posloupností do konce AGB, interpolované v ekvivalentních vývojových bodech | vlastní interpolace; data Choi et al. (2016) |
| `gridclip.py` | jednotný postup pro hodnoty mimo mřížky tabelovaných modelů | vlastní |
| `massloss.py` | hvězdný vítr v dynamice: hmotnosti z vývojových drah nebo rychlostí, index ztráty hmoty Ψ | vlastní |
| `masstransfer.py` | přetečení Rocheova laloku, stabilita, akrece, odezva dráhy, omlazení, kontaktní dvojhvězdy, přenos z třetí hvězdy | vlastní (předpisy Hurley et al. 2002 a dalších autorů) |
| `commonenv.py` | výsledek společné obálky (energetický formalismus) a odezva vnějších drah | vlastní |
| `merger.py` | pokračování běhu po splynutí dvou hvězd | vlastní |
| `tides.py`, `apsidal.py` | analytické apsidální stáčení (slapové, rotační, relativistické) a slapové časové škály | vlastní |
| `inertial_waves.py` | disipace inerciálních vln v konvektivních obálkách | vlastní |
| `magnetic.py` | magnetické brzdění a Applegateův energetický test | vlastní |
| `pulsation.py` | periody pulzací z modelů hvězd a modulace toku | vlastní |
| `mesa.py` | vlastní mřížky modelů hvězd: zapíše vstup pro MESA, spustí MESA, zredukuje profily na $k_2$, obálky a disipaci | obal MESA (externí, instaluje se zvlášť) + vlastní redukce |
| `gw.py`, `environment.py`, `tidalcapture.py`, `pdot.py` | rozpad drah gravitačními vlnami, slapy Galaxie a míjející hvězdy, slapový záchyt při setkáních, vnucená změna periody | vlastní |
| **Dynamika** | | |
| `backends/nbody.py` | N-částicový backend: sestavení simulace, události, extrémy mezi výstupy, rotace | vlastní + externí (REBOUND, REBOUNDx) |
| `backends/secular.py`, `backends/secular_nested.py`, `backends/secular_common.py` | sekulární backendy pro 2+2 a pro libovolnou vnořenou hierarchii | vlastní (integrátor SciPy) |
| `backends/evolving.py`, `backends/chunked.py` | vyvíjející se hvězdy a ztráta hmoty operátorovým štěpením; integrace po dávkách, které se vejdou do paměti | vlastní |
| `secular/twotwo.py`, `secular/kernel.py`, `secular/nested.py`, `secular/quadrupole.py` | sekulární rovnice (2+2, libovolný vnořený strom), jejich přeložené jádro, kvadrupólové předpovědi | vlastní (jádro přeložené pomocí Numba, je-li nainstalována) |
| **Diagnostiky a stabilita** | | |
| `diagnostics/hierarchy.py`, `zlk.py`, `binary.py`, `extremes.py`, `orbits.py` | rozpoznání hierarchie, ZLK a precese, apsidální stáčení a podíly zákrytů, extrémy, vzájemné sklony | vlastní |
| `stability.py` | kritéria Mardlinga–Aarsetha a Vynatheyi et al. (2022), N-částicový test duchem, volitelné neuronové klasifikátory (Vynatheya et al. 2023) | vlastní; klasifikátory jsou publikované sítě autorů (volitelné) |
| `chaosmap.py` | mapy MEGNO v jednom nebo dvou parametrech | vlastní řízení; MEGNO počítá REBOUND |
| `scatter.py`, `escape_theory.py` | rozptylové experimenty několika těles; statistická teorie rozpadu soustavy tří těles (Stone & Leigh 2019) | vlastní (rozptyl integrovaný pomocí REBOUND) |
| **Pozorovatelné veličiny** | | |
| `timing.py` | okamžiky zákrytů z dynamiky a jejich O−C | vlastní |
| `photometry.py` | fotometrická pásma, okrajové ztemnění, zákryt disků s okrajovým ztemněním, syntetické světelné křivky | vlastní; tabulky okrajového ztemnění (Claret & Bloemen 2011; Reeve & Howarth 2016) |
| `roche.py`, `mesh.py`, `shapes.py` | Rocheova geometrie, světlo deformovaných hvězd na povrchových sítích, volba tvarů hvězd | vlastní |
| `discs.py` | prachové disky a prstence jako zakrývající tělesa, jejich precese, tepelné a rozptýlené světlo | vlastní |
| `phenom.py` | fenomenologické profily zákrytů (Mikulášek 2015) a okamžiky minim změřené na syntetických světelných křivkách | vlastní |
| `radvel.py`, `rvlight.py` | radiální rychlosti těžišť a středů světla | vlastní |
| `phoebe_lc.py` | světlo vnitřních dvojhvězd spočtené programem PHOEBE 2 v hybridu s N-částicovými polohami | obal PHOEBE 2 (volitelný) |
| `systemview.py` | modelová soustava jako data pro obrázky | vlastní |
| **Prokládání naměřených dat** | | |
| `ocfit.py` | krok 1: analytický model O−C každého páru, metoda nejmenších čtverců, MCMC, předpovědi | vlastní + externí (nejmenší čtverce ze SciPy, emcee, převody času z astropy) |
| `dynfit.py` | krok 2: společný N-částicový model všech párů, škálování chyb, odlehlé body, apriorní rozdělení, ukotvení, MCMC, předpovědi | vlastní + externí (nejmenší čtverce ze SciPy, emcee) |
| `obsmodel.py`, `lcfit.py`, `incfit.py` | modelové pozorovatelné veličiny jedné integrace; světelné křivky a sklony ve společném proložení | vlastní |
| **Ansámbly a strojové učení** | | |
| `ensemble.py`, `parallel.py` | losování členů, jejich paralelní běh, souhrnné řádky, váhy podle měření, zprávy | vlastní |
| `emulator.py` | neuronové emulátory a klasifikátory natrénované na ansámblu | vlastní příznaky a vyhodnocení; sítě ze scikit-learn |
| `sbi_fit.py` | neuronový odhad aposteriorního rozdělení na ansámblu | obal `sbi` a PyTorch |
| **Produkty, řízení a zobrazení** | | |
| `io.py`, `runproducts.py` | soubor HDF5 běhu; souhrn, CSV, obrázky PNG vedle něj | vlastní (h5py, matplotlib) |
| `viewer.py`, `report.py` | interaktivní 3D prohlížeč běhu; samostatná zpráva HTML každého příkazu | vlastní (grafy v plotly) |
| `cli.py`, `__main__.py` | příkazový řádek | vlastní |
| `worker.py` | zpracovatel úloh | vlastní |
| `gui/server.py`, `gui/digest.py`, `gui/assist.py`, `gui/modelviews.py`, `gui/app.html`, `gui/help.yaml` | webová aplikace: server, stránky výsledků, nápovědy a rychlé kontroly, modelové křivky, samotná stránka | vlastní (webový server ze standardní knihovny Pythonu; plotly.js v prohlížeči) |
| `cnoise.py` | filtr hlášení, která vypisují přeložené knihovny | vlastní |
| **Ověřování** | | |
| `validation.py` a přidružené moduly, `_benchmark.py` | srovnávací testy spouštěné příkazem `swingersim validate` | vlastní |

**Tabulka 5.** Mapa modulů SWINGERsimu. Tučné řádky jsou skupiny. Prostřední sloupec uvádí úkol
každého modulu, pravý sloupec, zda je výpočet vlastním kódem SWINGERsimu, obalem, nebo vlastním kódem
s jednou externí částí (uvedenou v závorce). Cesty jsou relativní vůči `src/swingersim/`.

Balík obsahuje asi 45 000 řádků Pythonu a k tomu webovou stránku. Vlastní kód pokrývá fyziku, která
propojuje dynamiku, hvězdy a pozorovatelné veličiny; externí knihovny slouží pro přesně vymezené
numerické úkoly.

## Co SWINGERsim implementuje a co používá {#own-vs-external}

Tabulka 6 obojí výslovně odděluje. Odpovídá na otázku, kterou si oponent položí jako první: které
výsledky závisejí na kódu napsaném pro SWINGERsim a které na zavedeném externím softwaru
a publikovaných datech.

| Úkol | Implementováno ve SWINGERsimu | Použitý externí software nebo data |
| --- | --- | --- |
| hierarchie | sestavení libovolného stromu dvojic, inicializace v Jacobiho souřadnicích, rozpoznání hierarchie v každém výstupu a jejích změn | – |
| N-částicová integrace | sestavení, události (srážky, přetečení Rocheova laloku, úniky, skutečný čas výpočtu), vzorkování extrémů mezi výstupy, vazba na vyvíjející se hvězdy a změnu hmotnosti, rotace | REBOUND s IAS15 (Rein & Liu 2012; Rein & Spiegel 2015); síly REBOUNDx `gr_full`, `tides_spin` a operátor `modify_mass` (Tamayo et al. 2020; Lu et al. 2023) |
| sekulární integrace | rovnice pro 2+2 (Vokrouhlický 2016) a pro libovolný vnořený strom (Hamers et al. 2015), slapy, tlumení plynem, slapy Galaxie a míjející hvězdy, integrace po dávkách | integrátor SciPy DOP853; překladač Numba (volitelný) |
| hvězdy | interpolace vývojových drah, konstanty stavby, obálky, slapová zpoždění, omlazení, periody pulzací | vývojové dráhy: MIST (Choi et al. 2016), mřížka SSE z rutin Hurleyho et al. (2000); stavba: Claret (2004, 2023); ZAMS: Tout et al. (1996); volitelné vlastní mřížky spočtené programem MESA (Paxton et al. 2011; Jermyn et al. 2023) |
| interagující dvojhvězdy | vazba přenosu hmoty, kontaktních fází, společných obálek, přenosu z třetí hvězdy a splynutí na dynamiku | předpisy z literatury (Hurley et al. 2002 a další; *Fyzika a metody*) |
| zákryty a světelné křivky | geometrie zákrytů, zákryt disků s okrajovým ztemněním, členy elipsoidální proměnnosti, beamingu a odrazu, Rocheovy a rotující povrchy na sítích, prachové disky, doba letu světla, pulzace | tabulky okrajového ztemnění (Claret & Bloemen 2011) a koeficienty TLUSTY (Reeve & Howarth 2016); křivky propustnosti pásem; PHOEBE 2 (Prša et al. 2016) jako volitelný výpočetní modul pro vnitřní dvojhvězdy |
| okamžiky zákrytů a O−C | konjunkce nalezené v integraci, doba letu světla, O−C, fenomenologická minima (Mikulášek 2015) | – |
| radiální rychlosti | rychlosti v okamžiku vyzáření, střed světla s Rossiterovým–McLaughlinovým efektem | – |
| řetězec proložení | modely kroků 1 a 2, škálování chyb a jitter, vyřazování odlehlých bodů, apriorní rozdělení, faktory červeného šumu, ukotvení, předpovědi s nejistotami | nejmenší čtverce ze SciPy; emcee (Foreman-Mackey et al. 2013) pro MCMC; astropy pro převody času |
| ansámbly | losování, paralelní běhy, navazující řádky, váhy, efektivní počet členů, Wilsonovy intervaly, převažování | pandas a pyarrow jen pro zápis Parquet (volitelné) |
| strojové učení | příznaky, trénovací množiny, vyhodnocení, kalibrační grafy | sítě scikit-learn; `sbi` a PyTorch (volitelné) |
| stabilita | kritéria Mardlinga–Aarsetha a Vynatheyi et al. (2022) pro každou vnořenou trojici, N-částicový test duchem, mapy chaosu, rozptyl, teorie úniku | MEGNO z REBOUND; neuronové klasifikátory publikované v práci Vynatheya et al. (2023) (volitelné, stahují se zvlášť) |
| konfigurace | schéma, výchozí hodnoty, kontroly, chybová hlášení | pydantic (validační rámec), PyYAML |
| výstupy | struktura HDF5, produkty, 3D prohlížeč, zprávy HTML | h5py, matplotlib, plotly |
| webová aplikace, zpracovatel úloh | server, stránky výsledků, nápověda, fronta úloh, pozastavení a zastavení | standardní knihovna Pythonu; plotly.js v prohlížeči |

**Tabulka 6.** Co SWINGERsim implementuje sám (prostřední sloupec) a na které externí knihovny
a publikovaná data se spoléhá (pravý sloupec), úkol po úkolu. Pomlčka znamená, že kromě NumPy a SciPy
se neuplatňuje žádná externí složka.

## Zásady návrhu {#design}

- **Především dopředný model.** Jádro pouze předpovídá. Proložení, váhy ansámblu a strojové učení
  jsou vrstvy, které jádro volají a nikdy ho nemění.
- **Dva backendy, jeden výsledek.** N-částicový backend je referenční; sekulární backend je rychlá
  aproximace; oba vytvářejí stejnou strukturu, takže jeden může kontrolovat druhý.
- **Libovolná hierarchie jako strom.** Nic v řetězci nepředpokládá čtyři hvězdy: předvolby jsou
  stromy, diagnostiky procházejí strom, stabilita se testuje na každé vnořené trojici.
- **Ověřená fyzika.** Každá fyzikální složka je srovnána s publikovaným výsledkem (*Testování,
  část I*).
- **Reprodukovatelnost.** Každý soubor běhu ukládá svou úplnou konfiguraci a verze softwaru, který
  ho spočítal; ansámbl je určen svým semínkem; při navázání ansámblu se kontroluje otisk konfigurace,
  takže se nikdy nesmíchají řádky různých nastavení.
- **Bezpečné vzdálené výpočty.** Zpracovatel úloh spouští jen příkazy SWINGERsimu s argumenty, které
  jejich parser povoluje, a jen s cestami uvnitř složky SWINGERsimu; soubor úlohy nemůže spouštět
  příkazy shellu.
