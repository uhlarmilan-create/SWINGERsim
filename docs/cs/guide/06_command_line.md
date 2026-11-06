# Příkazový řádek {#command-line}

Vše, co SWINGERsim počítá, se spouští jedním příkazem, `swingersim`, za nímž následuje podpříkaz
a jeho argumenty. `swingersim --help` vypíše podpříkazy a `swingersim <sub-command> --help` volby
jednoho z nich. Webová aplikace i zpracovatel úloh sestavují své formuláře a příkazové řádky úloh ze
stejných definic, takže volba popsaná zde je k dispozici všude.

## Společné zvyklosti {#cli-conventions}

- Prvním argumentem většiny podpříkazů je konfigurační soubor (YAML). Relativní cesty uvnitř
  konfigurace (tabulky minim, světelné křivky, výstupní složky proložení) jsou relativní vůči
  konfiguračnímu souboru.
- `-o` / `--output` určuje výstupní soubor nebo složku; chybějící složky se vytvoří.
- Každý podpříkaz, který vytváří výsledky, zapíše také **zprávu HTML**: jeden samostatný soubor
  s příkazovým řádkem, konfigurací, vypsaným souhrnem, interaktivními grafy a obrázky, který se
  otevře bez připojení v libovolném prohlížeči. Zapisuje se vedle výstupu; `--html PATH` volí jiný
  název a `--no-html` ji vypne. `gui`, `worker`, `view`, `pulsation`, `applegate`, `mesa-grid`
  a `mesa-reduce` zprávu nezapisují.
- `-j` / `--processes` nastavuje počet paralelních procesů u příkazů, které mohou využít více jader;
  `0` znamená všechna jádra (jádra úlohy, spouští-li ji zpracovatel úloh).

| Podpříkaz | Účel | Hlavní vstup | Hlavní výstup |
| --- | --- | --- | --- |
| `run` | jedna simulace | konfigurace | souhrn, HDF5, CSV, PNG, 3D prohlížeč, zpráva |
| `view` | 3D prohlížeč uložených běhů | soubory nebo složky HDF5 | prohlížeč HTML |
| `stability` | kritéria stability a N-částicový test duchem | konfigurace | zpráva |
| `timings` | okamžiky zákrytů z dynamiky a jejich O−C | konfigurace | CSV, PNG |
| `minima` | minima změřená na syntetické světelné křivce | konfigurace | CSV, PNG |
| `lightcurve` | syntetická světelná křivka | konfigurace | CSV, PNG |
| `rv` | křivky radiálních rychlostí | konfigurace | CSV, PNG |
| `apsidal` | apsidální stáčení rozložené na jednotlivé příspěvky | konfigurace | tabulka |
| `pulsation` | periody pulzujících hvězd | konfigurace | vypsaná tabulka |
| `fit` | proložení naměřených minim, radiálních rychlostí, světelných křivek, sklonů | konfigurace s oddílem `fit` | parametry, rezidua, proložená konfigurace, předpovědi |
| `ensemble` | ansámbl Monte Carlo | konfigurace s oddílem `ensemble` | složka s řádky a souhrnem |
| `ensemble-report` | statistika dokončeného ansámblu, vážená měřeními | složka ansámblu | zpráva, obrázek |
| `emulator-train`, `emulator-report` | neuronové emulátory ansámblu | složka ansámblu | síť, souhrn aposteriorního rozdělení |
| `sbi` | simulační inference na ansámblu | složka ansámblu | vzorky aposteriorního rozdělení, korelační graf |
| `chaosmap` | mapa MEGNO v jednom nebo dvou parametrech | konfigurace | mapa, tabulka |
| `scatter` | rozptylové experimenty několika těles | konfigurace rozptylu | tabulka výsledků |
| `escape-theory` | statistický rozpad soustavy tří těles | čísla | rozdělení výsledků |
| `applegate` | může magnetická aktivita vysvětlit vlnu v O−C? | čísla | poměr energií |
| `validate` | srovnávací testy proti publikovaným výsledkům | identifikátory testů | verdikty, zpráva |
| `mesa-grid`, `mesa-reduce` | vlastní mřížky modelů hvězd programem MESA | soubor mřížky | běhy MESA, zredukovaná mřížka |
| `worker` | spouští soubory úloh ze složky | složka úloh | výsledky úloh |
| `gui` | webová aplikace | složka SWINGERsimu | – |

**Tabulka 17.** Všechny podpříkazy. Sloupce: název, k čemu slouží, co čte a co zapisuje. Oddíly níže
uvádějí volby každého z nich a příklad.

## Simulace {#cli-sims}

### `swingersim run` {#cli-run}

Integruje jednu konfiguraci po dobu `run.t_max_yr`, vypíše souhrn a s volbou `-o` zapíše soubor
běhu a jeho standardní produkty (kapitola [Výstupy](#out-run)).

| Volba | Význam |
| --- | --- |
| `config` | konfigurační soubor |
| `-o`, `--output` | výstupní soubor HDF5; produkty se zapíší vedle něj |
| `--backend nbody` nebo `secular` | nahradí `run.backend` z konfigurace |
| `--t-max-yr T` | nahradí `run.t_max_yr`; záporná hodnota integruje do minulosti |
| `--no-viewer` | nezapisovat 3D prohlížeč `<output>_viewer.html` |
| `--no-plots` | nezapisovat obrázky PNG `<output>_*.png` |

**Tabulka 18.** Volby příkazu `swingersim run`.

```bash
swingersim run configs/examples/example_2p2.yaml -o output/example.h5
swingersim run configs/examples/example_2p2.yaml --backend secular --t-max-yr -1000 -o output/past.h5
```

### `swingersim view` {#cli-view}

Zapíše interaktivní 3D prohlížeč uložených běhů: jednoho běhu, několika, nebo všech běhů ve
složkách. `run -o` ho zapisuje už sám; `view` ho sestaví znovu, např. s více snímky animace.

| Volba | Význam |
| --- | --- |
| `h5 …` | soubory HDF5 z příkazu `run`, nebo složky, v nichž se hledají |
| `-o`, `--output` | soubor HTML jednoho běhu (výchozí `<run>_viewer.html`) |
| `--frames N` | počet snímků animace (výchozí 150) |
| `--force` | u složek: znovu sestavit každý prohlížeč (výchozí: jen chybějící a starší než jejich běh) |
| `--cdn` | načíst plotly.js z internetu: menší soubor, který potřebuje připojení |

**Tabulka 19.** Volby příkazu `swingersim view`.

```bash
swingersim view output/example.h5 --frames 300
swingersim view output                # every run under output/ that has no up-to-date viewer
```

### `swingersim stability` {#cli-stability}

Vyhodnotí stabilitu každé vnořené trojice hierarchie kritérii Mardlinga & Aarsetha (2001)
a Vynatheyi et al. (2022), pro každou dílčí čtyřhvězdu neuronovými klasifikátory Vynatheyi et al.
(2023), jsou-li nainstalovány, a provede N-částicový test duchem podle Vynatheyi et al. (2023)
(*Fyzika a metody*, kapitola Stabilita a osudy).

| Volba | Význam |
| --- | --- |
| `--orbits N` | počet vnějších oběhů N-částicového testu (výchozí 100) |
| `--no-nbody` | vynechat N-částicový test |

**Tabulka 20.** Volby příkazu `swingersim stability`.

```bash
swingersim stability configs/examples/example_2p2_zlk.yaml
```

## Pozorovatelné veličiny {#cli-observables}

### `swingersim timings` {#cli-timings}

Najde každý zákryt vnitřních dvojhvězd v N-částicové integraci (konjunkci obou hvězd páru, s dobou
letu světla z těžiště páru) a zapíše okamžiky a jejich O−C vůči lineární efemeridě proložené
primárními zákryty.

| Volba | Význam |
| --- | --- |
| `--years Y` | konec řady [yr po $t = 0$] (povinná) |
| `--from-years Y0` | začátek řady [yr, ≤ 0]; záporné hodnoty přidají zákryty před $t = 0$, integrované do minulosti |
| `--bjd-t0 BJD` | BJD modelového času $t = 0$: přidá sloupec `bjd_obs` |
| `-o`, `--output` | soubor CSV (povinná) |
| `--binaries …` | názvy drah (výchozí: všechny dvojhvězdy ze dvou hvězd) |
| `--samples N` | počet vzorků na periodu dvojhvězdy (výchozí 40) |
| `--plot` | PNG s diagramem O−C |

**Tabulka 21.** Volby příkazu `swingersim timings`.

```bash
swingersim timings configs/examples/example_timings_czev4315.yaml --years 5 \
  -o output/timings_czev4315.csv --plot output/timings_czev4315.png
swingersim timings configs/squadra_test3/v907_sco_published.yaml --from-years -10 --years 10 \
  --bjd-t0 2459360.0 -o output/v907_timings.csv
```

### `swingersim lightcurve` {#cli-lightcurve}

Spočítá tok všech hvězd ve fotometrickém pásmu z N-částicových poloh: okrajové ztemnění, vzájemné
zákryty každé dvojice hvězd, elipsoidální proměnnost, dopplerovský beaming, odraz, dobu letu světla,
pulzace, rotující a Rocheovy hvězdy a prachové disky podle konfigurace (*Fyzika a metody*, kapitola
Pozorovatelné veličiny).

| Volba | Význam |
| --- | --- |
| `--days D` | délka řady [d] (povinná) |
| `--start-day S` | první okamžik [d po $t = 0$] (výchozí 0) |
| `--cadence-min C` | vzorkování [min] (výchozí 2) |
| `--band B` | fotometrické pásmo: `U`, `B`, `V`, `R`, `I`, `TESS`, `Gaia_G`, `Gaia_BP`, `Gaia_RP` (výchozí `TESS`) |
| `--exposure-min E`, `--supersample N` | expoziční doba [min] a počet dílčích expozic na bod, pro rozmazání zákrytů |
| `--noise-ppm N`, `--seed S` | bílý šum [ppm] a jeho semínko |
| `--ld FILE` | soubor YAML: hvězda → čtyři koeficienty okrajového ztemnění, místo tabulek |
| `--no-proximity` | bez elipsoidální proměnnosti, beamingu a odrazu |
| `--no-light-time` | zanedbat dobu letu světla |
| `--engine swingersim` nebo `phoebe` | světlo vnitřních dvojhvězd z vlastního modelu SWINGERsimu (výchozí) nebo z PHOEBE 2 |
| `--phoebe-atm`, `--phoebe-distortion`, `--phoebe-ntriangles`, `--no-irradiation` | nastavení PHOEBE: atmosféra (`auto`, `ck2004`, `blackbody`), tvar (`roche`, `sphere`), velikost sítě (výchozí 1500), bez odrazu |
| `-o`, `--output` | soubor CSV (povinná) |
| `--plot` | PNG se světelnou křivkou |

**Tabulka 22.** Volby příkazu `swingersim lightcurve`.

```bash
swingersim lightcurve configs/squadra_test3/v907_sco_published.yaml --days 30 --cadence-min 10 \
  -o output/v907_lc.csv --plot output/v907_lc.png
```

### `swingersim minima` {#cli-minima}

Spočítá syntetickou světelnou křivku a každé minimum na ní změří fenomenologickým profilem zákrytu
podle Mikuláška (2015), s chybou, tak jak by pozorovatel měřil skutečná data. Přijímá volby světelné
křivky z tabulky 22 (`--band`, `--cadence-min`, `--noise-ppm`, `--seed`, `--no-proximity`,
`--no-light-time`, volby PHOEBE) a navíc:

| Volba | Význam |
| --- | --- |
| `--years Y` | délka řady [yr] (povinná) |
| `-o`, `--output` | soubor CSV s minimy (povinná) |
| `--binaries …`, `--samples N` | názvy drah; počet vzorků na periodu dvojhvězdy |
| `--window W` | šířka prokládaného okna v jednotkách trvání zákrytu (výchozí 2,5) |
| `--season-days S`, `--seasons FILE` | proložit jeden tvar zákrytu pro každou sezónu této délky a zapsat tvary do CSV |
| `--blended` | měřit na celkovém světle (výchozí: každou dvojhvězdu na jejím vlastním světle plus třetím světle) |
| `--bjd-t0 BJD` | BJD_TDB času $t = 0$ pro časový sloupec (výchozí `fit.stage2.bjd_t0`, jinak J2000) |
| `--plot` | PNG s O−C minim a konjunkcí |

**Tabulka 23.** Volby příkazu `swingersim minima` nad rámec voleb příkazu `lightcurve`.

```bash
swingersim minima configs/benchmarks/qzcar_minima.yaml --years 2 --binaries A \
  -o output/qzcar_minima.csv --plot output/qzcar_minima.png
```

### `swingersim rv` {#cli-rv}

Křivky radiálních rychlostí hvězd z N-částicového modelu, v okamžicích vyzáření, s naměřenými
rychlostmi pro srovnání.

| Volba | Význam |
| --- | --- |
| `--start`, `--end` | první a poslední okamžik [BJD_TDB] (povinné) |
| `--bjd-t0` | BJD_TDB času $t = 0$ (výchozí `fit.stage2.bjd_t0`) |
| `--step` | krok křivek [d] (výchozí 0,01) |
| `--stars …` | hvězdy (výchozí všechny) |
| `--gamma` | systémová rychlost [km/s] (výchozí `fit.rv.gamma`) |
| `--measured FILE` | CSV s naměřenými rychlostmi (výchozí `fit.rv.file`) |
| `-o`, `--output` | výstupní složka (výchozí `output/rv`) |

**Tabulka 24.** Volby příkazu `swingersim rv`.

```bash
swingersim rv configs/squadra_test3/v907_sco_fit_lc_v2.yaml --start 2459700 --end 2459710 -o output/v907_rv
```

### `swingersim apsidal` {#cli-apsidal}

Stáčení apsid každé dvojhvězdy ze dvou hvězd rozložené na jednotlivé příspěvky: slapovou a rotační
deformaci obou hvězd podle Sterna (1939) pro synchronní, pseudosynchronní a nastavenou rotaci
a relativistické stáčení; volitelně celkové stáčení změřené v N-částicovém běhu, které přidává
poruchy od ostatních hvězd; a $k_2$, které by bylo potřeba pro cílovou periodu apsidálního stáčení.

| Volba | Význam |
| --- | --- |
| `--target A=194 B=98` | cílové periody apsidálního stáčení [yr] pro jednotlivé páry |
| `--nbody-years Y` | spustit také N-částicový model na Y let a změřit celkové stáčení |
| `-o`, `--output` | soubor CSV s tabulkou |

**Tabulka 25.** Volby příkazu `swingersim apsidal`.

```bash
swingersim apsidal configs/squadra_test3/v907_sco_published.yaml -o output/v907_apsidal.csv
```

Tabulka uvádí pro každý pár periodu, excentricitu, obě hodnoty $k_2$, slapové stáčení, rotační
stáčení pro tři předpoklady o rotaci, relativistické stáčení, celkové stáčení v deg/yr a periodu
apsidálního stáčení $U$ v letech.

### `swingersim pulsation` {#cli-pulsation}

Vypíše pro každou hvězdu s `stars.<m>.pulsation` její typ a mód, periodu hlavního módu, $dP/dt$
z vyvíjejícího se modelu hvězdy, její poloměr a $T_\mathrm{eff}$ a poznámky k platnosti.

```bash
swingersim pulsation configs/examples/example_pulsating_triple.yaml
# m1: delta_scuti (F)  P = 0.0711228 d  dP/dt = 5.35e-06 s/yr  R = 2.03 Rsun  Teff = 7931 K
# m3: lpv (F)  P = 334.992 d  dP/dt = 135 s/yr  R = 348 Rsun  Teff = 2867 K
```

## Proložení naměřených dat {#cli-fit}

### `swingersim fit` {#cli-fit-cmd}

Proloží model měřeními uvedenými v oddílu `fit`: krok 1 je analytický model O−C každého zákrytového
páru, krok 2 společný N-částicový model celé soustavy proložený všemi minimy najednou, volitelně
s radiálními rychlostmi, světelnými křivkami a sklony (postup v kapitole
[Pracovní postupy](#wf-fit), klíče v kapitole [Konfigurační soubor](#cfg-fit)).

| Volba | Význam |
| --- | --- |
| `config` | konfigurace s oddílem `fit` |
| `-o`, `--output` | výstupní složka (výchozí `fit.output.dir`) |
| `--stage 1` nebo `2` | 1: analytické O−C (výchozí bez `fit.stage2`); 2: společný N-částicový model (výchozí s ním) |
| `--processes N` | krok 2: jádra pro výpočet derivací (0 = všechna; výchozí `fit.stage2.processes`) |
| `--sampler none` nebo `emcee` | nahradí vzorkovač spouštěného kroku |
| `--anchor` | jen ukotvit periodu a fázi každého zákrytového páru k minimům a zapsat `anchored.yaml` |
| `--bjd-t0 BJD` | s `--anchor`: BJD_TDB času $t = 0$ (výchozí `fit.stage2.bjd_t0`) |
| `--no-plots` | nekreslit obrázky |

**Tabulka 26.** Volby příkazu `swingersim fit`. `--sampler`, `--processes` a `--stage` se na
`--anchor` nevztahují; `--processes` platí jen pro krok 2.

```bash
swingersim fit configs/examples/example_fit.yaml                       # stage 1, seconds
swingersim fit configs/examples/example_fit.yaml --sampler emcee        # stage 1 with MCMC
swingersim fit configs/squadra_test/tic120362137_fit.yaml --processes 0 # stage 2
swingersim fit configs/squadra_test/bu_cmi_fit.yaml --anchor            # anchoring only
```

### `swingersim applegate` {#cli-applegate}

Může cyklus magnetické aktivity hvězdy vysvětlit periodickou změnu O−C její dvojhvězdy (Applegate
1992)? Příkaz udá největší poloviční amplitudu O−C, kterou hvězda dokáže během periody modulace
pohánět v modelu tenké slupky podle Tiana et al. (2009), a pro pozorovanou amplitudu poměr potřebné
a dostupné energie (pod 1: možné).

| Volba | Význam |
| --- | --- |
| `--mass`, `--radius`, `--teff` | aktivní hvězda [M☉, R☉, K] |
| `--a`, `--pbin` | vzdálenost složek [R☉] a perioda [d] dvojhvězdy |
| `--pmod` | perioda modulace [yr] |
| `--oc` | pozorovaná poloviční amplituda O−C [s] (volitelná) |

**Tabulka 27.** Volby příkazu `swingersim applegate`; všechny kromě `--oc` jsou povinné.

```bash
swingersim applegate --mass 0.93 --radius 0.96 --teff 5040 --a 3.3 --pbin 0.522 --pmod 30.5 --oc 137.2
# largest O - C semi-amplitude the star can power over 30.5 yr: 1145 s (0.0132 d; ...)
# observed 137.2 s: dP/P = 1.79e-06, dE/E_sec = 0.0144 -> possible
```

Příkladem je V471 Tau s hodnotami podle Völschowa et al. (2016).

## Ansámbly a inference {#cli-ensembles}

### `swingersim ensemble` {#cli-ensemble}

Vylosuje členy popsané v oddílu `ensemble`, každého paralelně spustí dopředu a/nebo dozadu a zapíše
jeden souhrnný řádek na člena a směr (postup v kapitole [Pracovní postupy](#wf-ensemble), výstup
v kapitole [Výstupy](#out-ensemble)).

| Volba | Význam |
| --- | --- |
| `-o`, `--output` | výstupní složka (povinná) |
| `-j`, `--processes` | paralelní procesy (0 = všechna jádra) |
| `--n N` | nahradí `ensemble.n` |
| `--backend` | nahradí `run.backend` |
| `--resume` | navázat na přerušený ansámbl (dokončení členové se zachovají) |
| `--retry-errors` | s `--resume`: znovu spustit členy, které skončily chybou |

**Tabulka 28.** Volby příkazu `swingersim ensemble`.

```bash
swingersim ensemble configs/examples/example_v994_ensemble.yaml -o output/v994 --n 40 -j 0
swingersim ensemble configs/examples/example_v994_ensemble.yaml -o output/v994 --n 400 -j 0 --resume
```

### `swingersim ensemble-report` {#cli-ensemble-report}

Statistika dokončeného ansámblu: podíly osudů, rozdělení největších excentricit, vzájemných sklonů,
zákrytů, rychlostí změny sklonů a verdiktů stability, zvlášť dopředu a dozadu, s 68% Wilsonovými
intervaly; jsou-li zadána pozorování, je každá statistika vážená a vylosované parametry se shrnou
před podmíněním a po něm.

| Volba | Význam |
| --- | --- |
| `directory` | složka ansámblu |
| `--observed FILE` | soubor YAML: souhrnný sloupec → `{value, sigma}` nebo meze; nahradí `ensemble.observed` |
| `--prior` | ignorovat pozorování (nepodmíněná zpráva) |
| `--thresholds …` | uvést podíl členů s $1-e_\mathrm{max}$ pod těmito hodnotami (výchozí 0,1, 0,01, 0,001) |
| `--plot FILE` | přehledový obrázek |
| `--x`, `--y` | vylosovaný parametr a souhrnný sloupec posledního panelu obrázku |

**Tabulka 29.** Volby příkazu `swingersim ensemble-report`.

```bash
swingersim ensemble-report output/v994 --plot output/v994/overview.png --x orbits.out.i --y B_didt0_deg_yr
```

### `swingersim emulator-train` a `emulator-report` {#cli-emulator}

`emulator-train` natrénuje neuronový emulátor (hluboký soubor sítí – vícevrstvých perceptronů ze
scikit-learn), který předpovídá souhrnné sloupce ansámblu z jeho vylosovaných parametrů, nebo
klasifikátor sloupce typu ano/ne, a vypíše jeho skóre na odložených členech. `emulator-report`
losuje vstupy emulátoru z rozdělení v `ensemble.vary`, váží je podle `ensemble.observed` a vypíše
aposteriorní rozdělení každého vstupu (medián, 68% a 95% interval) – s miliony losování za
zanedbatelnou cenu.

| Volba | Význam |
| --- | --- |
| `ensemble`, `-o` | složka ansámblu a soubor modelu (`.pkl`) |
| `--inputs …` | vstupní sloupce (vylosované parametry) |
| `--targets …`, `--classes …` | sloupce pro regresi; binární sloupce nebo `column==value` (např. `fate==intact`) |
| `--log-inputs`, `--log-targets` | použít $\log_{10}$ |
| `--angle-inputs` | úhly ve stupních, použité jako cos a sin |
| `--scale-targets target=input^power` | prokládat cíl vydělený známou mocninnou závislostí |
| `--hidden`, `--models`, `--activation`, `--alpha`, `--learning-rate`, `--batch-size`, `--max-iter`, `--test-fraction`, `--seed` | nastavení sítě (výchozí: čtyři vrstvy po 50 jednotkách, 5 sítí, ReLU, L2 $10^{-4}$, rychlost učení $10^{-3}$, 2000 iterací, 20 % odloženo) |
| `--plot` | PNG s diagnostikou na trénovacích a odložených členech |
| `emulator-report model config`, `--n`, `--seed` | model, konfigurace ansámblu, počet apriorních losování (výchozí $10^6$) |

**Tabulka 30.** Volby příkazů emulátoru.

```bash
swingersim emulator-train output/v994 -o output/v994/emulator.pkl \
  --inputs orbits.out.i orbits.out.Omega orbits.B.Omega \
  --angle-inputs orbits.out.i orbits.out.Omega orbits.B.Omega \
  --targets B_didt0_deg_yr --plot output/v994/emulator.png
swingersim emulator-report output/v994/emulator.pkl my_v994_with_observed.yaml
```

Vstupy musejí být sloupce, které ansámbl losoval (`ensemble.vary`); `emulator-report` potřebuje
konfiguraci, jejíž oddíl `ensemble` má také položky `observed`. Emulátor je jen tak dobrý jako jeho
trénovací množina: při natrénování na 40 členech prvního příkladu `swingersim ensemble` výše je
skóre tohoto příkladu na odložených členech $R^2 = -0{,}25$, tj. horší než předpovídat průměr;
potřebuje tisíce členů.

### `swingersim sbi` {#cli-sbi}

Simulační inference: neuronový odhad aposteriorního rozdělení (balík `sbi`) se natrénuje na
dokončeném ansámblu a vrátí aposteriorní rozdělení zvolených parametrů pro naměřené hodnoty;
natrénovaný model se znovu použije pro nová měření se stejnými chybami.

| Volba | Význam |
| --- | --- |
| `directory` | složka ansámblu (trénovací množina) |
| `--observed FILE` | YAML: sloupec → `{value, sigma}` a/nebo `{min, max}` |
| `--params …`, `--log …` | parametry k odhadu; ty, které se odhadují v $\log_{10}$ |
| `--model FILE`, `--retrain` | soubor natrénovaného modelu (znovu použit, pokud odpovídá); trénovat znovu |
| `--epochs`, `-n` | nejvyšší počet epoch trénování (2000); počet vzorků aposteriorního rozdělení (20 000) |
| `--save-samples`, `--plot` | CSV se vzorky; korelační graf PNG |
| `--coverage`, `--diagnostics` | kalibrační test na 20 % odložených členů; PNG se ztrátovou funkcí, pokrytím a kalibrací |

**Tabulka 31.** Volby příkazu `swingersim sbi` (potřebuje skupinu `sbi`, tabulka 8).

```bash
cat > output/v994/observed.yaml <<'EOF'
B_didt0_deg_yr: {value: 0.05, sigma: 0.02}
EOF
swingersim sbi output/v994 --observed output/v994/observed.yaml \
  --params orbits.out.i orbits.out.Omega --model output/v994/sbi_model.pt \
  --plot output/v994/corner.png --coverage
```

Pozorovaná hodnota v tomto příkladu je jen ilustrativní; aposteriorní rozdělení je jen tak dobré jako
počet členů, na nichž byla síť natrénována (*Fyzika a metody*, kapitola Ansámbly a inference).

## Stabilita, chaos a rozptyl {#cli-dynamics}

### `swingersim chaosmap` {#cli-chaosmap}

Spustí jednu N-částicovou integraci v každém bodě mřížky jednoho nebo dvou konfiguračních parametrů
a spočítá ukazatel chaosu MEGNO (Mean Exponential Growth factor of Nearby Orbits; 2 pro pravidelný
pohyb, větší pro chaos).

| Volba | Význam |
| --- | --- |
| `-o`, `--output` | výstupní složka (povinná) |
| `--x PATH`, `--x-range MIN MAX`, `--nx N`, `--log-x` | první parametr (např. `orbits.out.P`), jeho rozsah, počet bodů (30), logaritmické rozestupy |
| `--y`, `--y-range`, `--ny`, `--log-y` | volitelný druhý parametr |
| `--orbits N` | délka v periodách nejširší dráhy (výchozí 1000) |
| `--t-max T` | délka v letech (nahradí `--orbits`) |
| `--max-wall S` | sekundy na bod mřížky (výchozí 300) |
| `-j`, `--resume` | paralelní procesy; navázat na přerušenou mapu |

**Tabulka 32.** Volby příkazu `swingersim chaosmap`.

```bash
swingersim chaosmap configs/validation/quarles18_circumbinary.yaml -o output/chaosmap_cb \
  --x orbits.out.a --x-range 1.8 4.0 --nx 23 --y orbits.out.e --y-range 0 0.5 --ny 11 --t-max 3000 -j 0
```

### `swingersim scatter` {#cli-scatter}

Rozptylové experimenty několika těles: dvojhvězda nebo osamocená hvězda se setká s jinou dvojhvězdou
nebo osamocenou hvězdou; setkání se integrují, dokud není výsledek jasný (průlet, výměna, ionizace,
vázaná trojhvězda nebo čtyřhvězda, srážka). Konfigurace má vlastní oddíl `scatter`
([Konfigurační soubor](#cfg-scatter)).

| Volba | Význam |
| --- | --- |
| `-o`, `--output` | výstupní složka (povinná) |
| `-j`, `--processes`, `--n`, `--resume` | paralelní procesy, počet setkání (nahradí `scatter.n`), navázání |

**Tabulka 33.** Volby příkazu `swingersim scatter`.

```bash
swingersim scatter configs/examples/example_scatter_2p2.yaml -o output/scatter_2p2 --n 200 -j 0
```

### `swingersim escape-theory` {#cli-escape}

Statistická teorie rozpadu nehierarchické soustavy tří těles (Stone & Leigh 2019): která hvězda
unikne a jakou excentricitu bude mít zbylá dvojhvězda, pro zadané hmotnosti, dvojhvězdu, rychlost
v nekonečnu a srážkový parametr (nebo moment hybnosti); volitelně srovnaná se složkou `scatter`
s běhy dvojhvězda–osamocená hvězda.

```bash
swingersim escape-theory --masses 1 1 1 --a 1 --v-inf 1 --b 2
```

Volby: `--masses m1 m2 m3` (hvězda 1 dvojhvězdy, hvězda 2 dvojhvězdy, osamocená hvězda, M☉), `--a`
[au], `--e`, `--v-inf` [km/s], `--b` [au] nebo `--L0-ratio`, `--alpha`, `--criterion AE|SE`,
`--compare FOLDER`.

## Ověřování a vlastní modely hvězd {#cli-other}

### `swingersim validate` {#cli-validate}

Reprodukuje publikované výsledky a pro každý srovnávací test vypíše PASS, EXPLAINED nebo FAIL;
`--all` zahrne i pomalé testy (každý několik minut); identifikátory testů vyberou jen některé. Zpráva
se zapíše do `output/validation/`. Srovnávací testy jsou předmětem dokumentu *Testování, část I*.

```bash
swingersim validate            # the fast benchmarks
swingersim validate --all      # all of them
```

### `swingersim mesa-grid` a `mesa-reduce` {#cli-mesa}

Spočítají vlastní mřížku modelů hvězd programem MESA ze souboru mřížky (`configs/mesa/*.yaml`)
a zredukují její profily na veličiny, které SWINGERsim potřebuje ([Instalace](#install-mesa)).
Volby: `--root` (složka SWINGERsimu), `-j` (paralelní běhy MESA), `--test` (jen jedna vývojová dráha
hvězdy 1 M☉).

## Spouštění a řízení {#cli-control}

### `swingersim worker` {#cli-worker}

Spouští soubory úloh vložené do složky ([Instalace](#install-worker),
[Pracovní postupy](#wf-worker)).

| Volba | Význam |
| --- | --- |
| `jobdir` | složka úloh (výchozí `<root>/jobs`) |
| `--root` | složka SWINGERsimu (výchozí `.`) |
| `--poll S` | sekundy mezi kontrolami (výchozí 10) |
| `--once` | spustit čekající úlohy a skončit |
| `--cores N` | počet jader, která mohou úlohy využívat společně (výchozí všechna) |

**Tabulka 34.** Volby příkazu `swingersim worker`.

### `swingersim gui` {#cli-gui}

Spustí webovou aplikaci; její volby jsou v tabulce 11.
