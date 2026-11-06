# Pracovní postupy {#workflows}

Tato kapitola popisuje typické úlohy od začátku do konce. Každý postup začíná tabulkou vstupů, kroků
a výstupů a končí krátkým řešeným příkladem. Příkazy jsou ty z kapitoly
[Příkazový řádek](#command-line); každý krok lze provést i ve webové aplikaci (kapitola
[Webová aplikace](#gui)).

## Dopředný model publikované soustavy {#wf-forward}

**Otázka.** Článek uvádí periody, zákryty, hmotnosti a případně vnější dráhu vícenásobné soustavy.
Co tato soustava dělá a co pozorovatelé uvidí příště?

| Vstupy | Kroky | Výstupy |
| --- | --- | --- |
| publikované hmotnosti, poloměry, teploty; elementy každé dráhy; epocha, k níž se vztahují; naměřené apsidální stáčení nebo změny sklonu pro srovnání | 1. napsat konfiguraci; 2. zkontrolovat stabilitu; 3. spustit N-částicově a pro dlouhé doby sekulárně; 4. spočítat okamžiky zákrytů a světelnou křivku; 5. porovnat s měřeními | soubory běhů, okna zákrytů, diagram O−C, světelná křivka, verdikt stability; rozdíly oproti měřením |

**Tabulka 35.** Postup dopředného modelu. Levý sloupec uvádí, co potřebujete, prostřední kroky,
pravý to, co dostanete.

1. **Shromážděte parametry** i s jejich zdroji: vnitřní periody, excentricity, sklony, poloměry
   (z řešení světelných křivek), hmotnosti (spektroskopické, jsou-li k dispozici), vnější dráhu
   (z LTTE, radiálních rychlostí nebo astrometrie), periody apsidálního stáčení a naměřené rychlosti změn.
   Poznamenejte si epochu elementů: stane se časem $t = 0$.
2. **Napište konfiguraci.** Vyjděte z příkladu v `configs/` (katalog publikovaných soustav je
   v `configs/catalog/`). Každou hodnotu, která není naměřená, označte komentářem (`# ASSUMED`,
   `# ESTIMATED`). Zadejte poloměry (`stars.<m>.R`), jsou-li známy; jinak se použijí poloměry
   hlavní posloupnosti nulového stáří, které jsou pro vyvinuté hvězdy příliš malé. Neznámé uzly:
   ponechte `Omega = 0` pro zákrytový pár a vzájemný sklon nastavte pomocí `i` a `Omega` ostatních
   drah. U těsných dvojhvězd zapněte `physics.gr`, `tides` a `rotation`; určují apsidální stáčení
   a mohou potlačit cykly ZLK.
3. **Zkontrolujte stabilitu**: `swingersim stability config.yaml`. Nestabilní nominální soustava
   znamená, že některá hodnota je chybná nebo soustava není taková, jak se předpokládalo.
4. **Spusťte** N-částicový backend na několik period precese a prohlédněte si stránku výsledku nebo
   prohlížeč. Vývoj za miliony let viz [níže](#wf-evolution).
5. **Co by viděli pozorovatelé?** `swingersim timings` dá O−C každé dvojhvězdy s dobou letu světla
   a dynamickými efekty; `--bjd-t0` přiřadí každému minimu jeho BJD a záporné `--from-years` přidá
   minima před epochou elementů. `swingersim lightcurve` dá světelnou křivku ve fotometrickém pásmu.
   Nesouhlasí-li naměřená minima s modelem, proložte je ([další oddíl](#wf-fit)).

**Řešený příklad.** V907 Sco je zákrytová dvojhvězda s periodou 3,78 d a třetí hvězdou na dráze
s periodou 142 d; její zákryty během 20. století ustaly a znovu začaly. Konfigurace
`configs/squadra_test3/v907_sco_published.yaml` přebírá řešení Zascheho et al. (2023):

```bash
swingersim stability configs/squadra_test3/v907_sco_published.yaml
swingersim run configs/squadra_test3/v907_sco_published.yaml -o output/v907/run.h5
swingersim timings configs/squadra_test3/v907_sco_published.yaml --from-years -100 --years 100 \
  -o output/v907/timings.csv
```

Kontrola stability dává poměry 1,67 (Mardling–Aarseth) a 2,2 (Vynatheya et al. 2022) a stabilní
test duchem na 100 vnějších oběhů. Běh na 200 let (27 s) ukazuje vzájemný sklon 25,6–26,3°, periodu
precese uzlů 63 let a zákryty vnitřního páru během 50,5 % času, v oknech dlouhých asi 30 let
(obrázky 9–11); okamžiky zákrytů ukazují tatáž okna (obrázek 13). Stránky výsledků těchto tří úloh
jsou ty, které ukazuje kapitola [Webová aplikace](#gui).

## Prokládání naměřených dat {#wf-fit}

**Otázka.** Které parametry soustavy reprodukují naměřené okamžiky minim všech zákrytových párů –
a jsou-li k dispozici, také radiální rychlosti, světelné křivky a sklony naměřené v různých
epochách?

| Vstupy | Kroky | Výstupy |
| --- | --- | --- |
| konfigurace s počátečními hodnotami; tabulka minim; volitelně radiální rychlosti, segmenty světelných křivek, sklony; BJD času $t = 0$; apriorní rozdělení z publikovaných řešení | 1. připravit tabulky; 2. krok 1 (analytické O−C každého páru); 3. ukotvení; 4. krok 2 (společný N-částicový model) se škálováním chyb, vyřazením odlehlých bodů, apriorními rozděleními, červeným šumem; 5. zkontrolovat rezidua a varování; 6. MCMC pro spolehlivé chyby; 7. předpovědi a proložená konfigurace | proložené parametry s chybami a korelacemi, rezidua každého měření, předpovězená minima s nejistotou, `stage2_fitted.yaml` |

**Tabulka 36.** Postup prokládání („řetězec proložení“). Kroky jsou popsány níže; každá volba je
klíčem oddílu `fit` ([Konfigurační soubor](#cfg-fit)).

### Data {#wf-fit-data}

**Okamžiky minim** jsou tabulka CSV s jedním řádkem na minimum a sloupci `pair` (název zákrytové
dráhy), `type` (`p` primární, `s` sekundární), `time`, `sigma` (obojí ve dnech), `source` a `method`.
Přijímají se také názvy sloupců souborů naměřených minim z webové aplikace (`binary`, `kind`,
`bjd_tdb`, `sigma_d`, `ref`). Časy jsou v BJD_TDB; `fit.time_system: HJD_UTC` nebo `JD_UTC` je
převede (s pomocí `fit.coordinates`). *Primární* minimum je zákryt první hvězdy páru (`m1` dráhy
`inner`) bez ohledu na jeho hloubku; je-li primární hvězdou pozorovatele ta druhá, prohoďte typy nebo
hvězdy.

**Radiální rychlosti** (`fit.rv`): sloupce `star`, `time` (BJD_TDB), `rv` a `sigma` (km/s)
a `source`. Každý zdroj může mít vlastní nulový bod (`fit.rv.offsets`, volný jako
`rv.offset.<source>`); systémová rychlost je `rv.gamma`. U polodotykových a kontaktních dvojhvězd
použije `centre: light` místo rychlosti těžiště tokem váženou rychlost viditelného povrchu, včetně
Rossiterova–McLaughlinova efektu během zákrytu.

**Světelné křivky** (`fit.lc`): sloupce `time`, `flux`, `sigma` a volitelně `segment`. Každý segment
dostane vlastní měřítko, nulový bod a polynomický trend, které se při každém vyhodnocení řeší
lineárně; s `scale: false` plynou hloubky z hvězd tak, jak jsou zadány, takže světelná křivka omezuje
sklon a poloměry.

**Sklony v epochách** (`fit.inclination`): sloupce `mid_hjd` (nebo `time`), `i_deg`, `sigma_deg`
a volitelně `limit` a `source`, např. sklony z řešení světelných křivek archivních dat. Porovnává se
jen $\lvert\cos i\rvert$, protože zákryty neprozradí smysl oběhu. Sloupec `limit` činí řádek
jednostranným: `upper`, `lower`, `no_eclipse` (model pak musí být pod svou mezí zákrytů) nebo
`eclipse` (musí být nad ní).

### Krok 1: analytické O−C každého páru {#wf-fit-stage1}

Krok 1 popisuje O−C každého páru analytickými členy: efektem světelného času vnější dráhy (Irwin
1952), apsidálním stáčením excentrického páru a kvadratickým členem (`fit.stage1.terms` nebo pro
jednotlivé páry `fit.stage1.pairs`). S `shared_outer: true` sdílejí všechny páry soustavy 2+2 jednu
vnější dráhu (společnou periodu, epochu a excentricitu, opačná pericentra). Krok je rychlý (sekundy)
a dává počáteční hodnoty a první náhled: vnější periodu, amplitudu a hmotnostní funkci, periodu
apsidálního stáčení. `sampler: emcee` přidá běh MCMC (Foreman-Mackey et al. 2013).

### Ukotvení {#wf-fit-anchor}

Krok 2 konverguje jen z počátku, který minima správně očísluje: chyba vnitřní periody nastřádaná za
dobu pokrytou daty musí zůstat hluboko pod polovinou periody. Publikované elementy
z fotodynamického řešení se v N-částicovém modelu s mírně odlišnou fyzikou často rozcházejí s daty.
*Ukotvení* (anchoring) upraví jen periodu $P$ a střední anomálii $M_0$ každého zákrytového páru tak,
aby modelová minima padla na naměřená; `swingersim fit --anchor` zapíše výsledek jako
`anchored.yaml` a `fit.stage2.anchor: true` to před proložením provede automaticky.

### Krok 2: společný N-částicový model {#wf-fit-stage2}

Krok 2 počítá okamžiky všech minim všech párů z jedné N-částicové integrace celé soustavy – se
stejnou fyzikou jako `run` a `timings` – dopředu a dozadu od $t = 0$ (BJD `fit.stage2.bjd_t0`)
a prokládá parametry uvedené v `fit.stage2.free` váženou metodou nejmenších čtverců. Efekt světelného
času, apsidální stáčení a dynamické poruchy jsou v modelu současně; žádný člen není třeba volit. Volnými
parametry jsou elementy drah (`orbits.<o>.P`, `M0`, `e`, `omega`, `Omega`, `i`, nebo `ecosw`
a `esinw` u téměř kruhových drah, a `Pdot`), hmotnosti (`masses.<m>`), vlastnosti hvězd
(`stars.<m>.R`, `Teff`, `apsidal_k2`, `rg2`) a nulové body radiálních rychlostí (`rv.gamma`,
`rv.offset.<source>`). Derivace se počítají konečnými diferencemi, dvě N-částicové integrace na
volný parametr a iteraci, paralelně na `processes` jádrech. Metodu popisuje dokument *Fyzika
a metody*, kapitola Prokládání pozorovaných dat.

Co minima omezují: periody, fáze, excentricity a argumenty pericentra všech drah a prostřednictvím
amplitudy efektu světelného času a dynamických poruch také vzájemný poměr hmotností složek. Sklony
zákrytových párů a absolutní hmotnosti vyžadují světelné křivky a radiální rychlosti; ponechte je
pevné na hodnotách takového řešení, prokládejte je spolu s těmito daty, nebo jim dejte apriorní
rozdělení. Bez astrometrie je měřitelný jen rozdíl uzlů dvou drah: uzel jedné dráhy ponechte pevný.

**Chyby a jitter.** Minima z mnoha zdrojů mají realistické chyby jen zřídka.
`error_scale: per_source` dá každému zdroji (sloupec `source`) vlastní korekci, takže jeho
redukované $\chi^2$ je 1, nikdy pod uvedenými chybami; `global` používá jeden faktor pro všechny.
Korekce je faktorem, kterým se násobí chyby (`error_model: scale`), nebo minimální chybou přičtenou
kvadraticky (`error_model: jitter`), vhodnou pro přesná měření okamžiků s optimistickými chybami.
Zdroje s méně než `min_source_minima` minimy sdílejí jednu korekci. Proložení se pak opakuje ze
svého řešení, dokud se korekce neustálí, nejvýše však tolikrát, kolik udává `refits`.

**Odlehlé body.** `clip_sigma: 5` po každém proložení odstraní minima, která se odchylují o více než
pět robustních směrodatných odchylek svého zdroje (nejvýše podíl `max_clip_fraction` z nich),
a proloží znovu; v tabulce reziduí zůstávají se stavem `outlier`. Případně robustní ztrátová funkce
(`loss: soft_l1`, `huber` nebo `cauchy`, od `f_scale` σ výše) sníží jejich váhu, aniž by je
odstranila.

**Apriorní rozdělení.** `priors: {orbits.out.e: [0.47, 0.05]}` přidá k $\chi^2$ pro volný parametr
gaussovský člen $((x-\mu)/\sigma)^2$, např. z publikovaného řešení světelné křivky nebo radiálních
rychlostí. Apriorní rozdělení brání slabě omezeným parametrům (vnější dráze viděné jen přes jeden
zákrytový pár, dráze efektu světelného času pokryté jen jednou) v toulání k řešením, která prokládají
minima stejně dobře, ale odporují jiným datům. Zpráva uvádí odchylku každého apriorního rozdělení
v jednotkách σ; několik σ znamená skutečné napětí mezi daty a apriorním rozdělením.

**Červený šum.** Rezidua světelné křivky bývají korelovaná (přístrojové trendy, skvrny, oscilace),
takže formální chyby jsou příliš malé. `red_noise: true` u položky světelné křivky vynásobí chyby
každého segmentu faktorem β metody časového průměrování (Pont et al. 2006; Winn et al. 2008): o kolik
víc se rozptylují rezidua průměrovaná přes 20–120 min, než by se rozptyloval bílý šum. Uplatňuje se
mezi jednotlivými proloženími v řadě, takže vyžaduje `error_scale` nebo `clip_sigma`.

**MCMC.** Chyby z kovarianční matice předpokládají přímé údolí $\chi^2$. Uvádí-li zpráva silně
korelované parametry ($\lvert r\rvert > 0{,}95$), spustí `sampler: emcee` afinně invariantní MCMC
(Goodman & Weare 2010; Foreman-Mackey et al. 2013) kolem řešení s `walkers` × `mcmc_steps`
vyhodnoceními modelu; první polovina se zahodí. S `mcmc_tau_target: 50` se řetězec prodlužuje, dokud
není dlouhý alespoň 50 integrovaných autokorelačních časů, jak doporučují autoři emcee, nebo dokud
nedosáhne `mcmc_max_steps`. Stránka proložení ukazuje poměr počtu kroků k autokorelačnímu času.

**Předpovědi.** `fit.output.predict: {years: 5}` zapíše příští minima každého páru s jejich
nejistotou: sady parametrů se losují ze vzorků MCMC nebo z kovarianční matice a každá se
zintegruje, takže chyby zahrnují korelace parametrů.

**Kontrola výsledku.** Čtěte v tomto pořadí: varování; $\chi^2$ od počátečního k proloženému modelu
a redukované $\chi^2$; rms každého druhu dat; střední reziduum primárních a sekundárních minim
(rozdíl několika σ znamená chybné $e\cos\omega$); soustavné vzory v reziduích (parametr, který by měl
být volný, je pevný, nebo proložení uvízlo v lokálním minimu); parametry na svých mezích; odchylky od
apriorních rozdělení.

**Řešený příklad.** Proložení V907 Sco v dokumentu *Testování, část III* kombinuje 58 minim z TESS,
172 radiálních rychlostí ze dvou zdrojů, 2562 bodů dvou segmentů světelné křivky z TESS a 14 sklonů
naměřených mezi lety 1900 a 2022, s 16 volnými parametry a 3 apriorními rozděleními:

```bash
swingersim fit configs/squadra_test3/v907_sco_fit_lc_v2.yaml
```

Jeho konfigurace ukotví vnitřní pár (`anchor: true`), opraví chyby jednotlivých zdrojů jitterem,
odstraní odlehlé body za 5σ a na světelné křivky uplatní faktory červeného šumu. Proložení sníží
$\chi^2$ z 11 920 na 1052 za 4384 vyhodnocení modelu; rms minim je 2,6 min, radiálních rychlostí
3,85 km/s a světelných křivek 0,666 ppt (obrázky 19–23).

## Ansámbly: co nevíme {#wf-ensemble}

**Otázka.** Některé parametry soustavy neznáme (neviděný průvodce, orientace vnější dráhy, smysl
oběhu). Které hodnoty jsou možné, co dovolují měření a co z toho plyne pro budoucnost a minulost?

| Vstupy | Kroky | Výstupy |
| --- | --- | --- |
| nominální konfigurace; rozdělení neznámých parametrů (`ensemble.vary`); naměřené hodnoty nebo meze (`ensemble.observed`); požadavky na stabilitu | 1. napsat oddíl `ensemble`; 2. spustit několik set členů lokálně, víc na zpracovateli úloh; 3. podmínit měřeními; 4. přečíst $n_\mathrm{eff}$, podmíněná rozdělení a výsledky; 5. otevřít nejlepší členy jako konfigurace | `summary.csv` (jeden řádek na člena a směr), zpráva, vážená rozdělení a podíly, nejlépe odpovídající členy |

**Tabulka 37.** Postup ansámblu.

![**Obrázek 25.** Myšlenka ansámblu, schematicky. Vlevo: neznámý parametr (zde vnější perioda ve
dnech) se mnohokrát vylosuje ze svého rozdělení (histogram losování). Uprostřed: každý člen se
zintegruje; křivky ukazují excentricitu vnitřní dvojhvězdy v závislosti na čase v letech pro několik
členů. Vpravo: rozdělení odpovědi přes všechny členy (zde největší dosažená excentricita). Vážení
měřeními popsané níže promění levý histogram v aposteriorní rozdělení a pravý v podmíněnou
předpověď.](../../figures/concepts/ensemble_idea.png)

1. **Co nevíme** (obrázek 25). Losujte jen to, co je neznámé; naměřené hodnoty ponechte nebo je
   losujte v rámci jejich chyb (`normal`). Typické volby: neznámá orientace je `isotropic`
   (rovnoměrně v $\cos i$); neznámý uzel nebo fáze je `angle` (rovnoměrně 0–360°); parametr měřítka,
   jako vnější perioda nebo hmotnost průvodce, je `loguniform`; zákrytový pár s neznámým smyslem oběhu
   dostane `mirror: true` ($i$ nebo $180° - i$). `require: [stable_ma01]` znovu vylosuje členy, které
   kritérium stability označí za nestabilní. Každou volbu si zapište: výsledek na ní závisí.
2. **Spusťte** dopředu i dozadu (`directions: [forward, backward]`): rychlosti změn „nyní“ se
   nejlépe měří z obou stran $t = 0$ a minulost je polovinou odpovědi.
3. **Váhy podle pozorování.** Každá položka `ensemble.observed` je souhrnný sloupec ansámblu buď
   s měřením `{value, sigma}`, s mezemi `{min, max}`, nebo s požadovanou textovou hodnotou
   `{equals: intact}`. Člen dostane váhu $w = \exp(-\chi^2/2)$ s
   $\chi^2 = \sum ((x - \mathrm{value})/\sigma)^2$ a $w = 0$, porušuje-li mez; vážené členy
   aproximují aposteriorní rozdělení losovaných parametrů (vážení podle důležitosti s apriorním
   rozdělením jako návrhovým). Současné rychlosti změn jsou ve výchozím stavu *centrální* – průměr
   lineárního proložení dopředu a dozadu –, tak jak je měří pozorovatel; předponou `forward:`,
   `backward:` nebo `central:` před názvem sloupce lze volit.
4. **Kolik zbývá.** Efektivní počet členů $n_\mathrm{eff} = (\sum w)^2/\sum w^2$ (Kish 1965) říká,
   kolika stejně váženým členům odpovídá vážený ansámbl. Pod zhruba 20 jsou podmíněná rozdělení
   šumem: spusťte více členů, nebo zužte losované rozsahy na oblast, kterou měření připouštějí.
5. **Evidence.** Střední váha členů, $Z = \langle w \rangle$, je pravděpodobnost měření za hypotézy
   ansámblu (evidence). Poměr $Z$ dvou ansámblů srovnává dvě hypotézy (např. 2+2 proti 2+(2+1)) –
   ale jen tehdy, jsou-li společné parametry losovány ze stejných apriorních rozdělení; jinak poměr
   měří apriorní rozdělení. `ensemble-report` vypisuje $Z$ s monte-carlovskou chybou
   $\mathrm{std}(w)/\sqrt{n}$ vedle $n_\mathrm{eff}$; když vyhovují méně než dva členy, upozorní,
   že se má uvádět horní mez $Z < 1/n$. Totéž několika řádky Pythonu:

   ```python
   import numpy as np
   from swingersim.config import load_config
   from swingersim.ensemble import load_rows, observation_weights

   cfg = load_config("configs/squadra_test3/qx_cas_ensemble.yaml")
   rows = load_rows("output/squadra_test3/qx_cas/ensemble")
   w = np.array(list(observation_weights(rows, cfg.ensemble.observed).values()))
   print("evidence Z =", w.mean(), " n_eff =", w.sum() ** 2 / (w ** 2).sum())
   ```

6. **Využijte členy.** Na stránce výsledku lze členy s největší vahou otevřít jako konfigurace
   a spustit se světelnou křivkou nebo okamžiky zákrytů. Pro rychlé aposteriorní rozdělení nových
   měření natrénujte na ansámblu emulátor nebo neuronový odhad aposteriorního rozdělení
   (`emulator-train`, `sbi`).

**Řešený příklad.** Zákrytová dvojhvězda QX Cas (6,0 d) mění svůj sklon o −0,23 ± 0,08 °/yr;
příčinou je neviděné třetí těleso. `configs/squadra_test3/qx_cas_ensemble.yaml` losuje jeho hmotnost
(0,3–6 M☉, log-rovnoměrně), periodu (30–11 000 d, log-rovnoměrně), excentricitu (0–0,6) a orientaci
jeho dráhy, spustí 2000 členů na 30 let dopředu i dozadu a váží je podle současného sklonu
(82 ± 1,5°), jeho centrální rychlosti změny a přežití soustavy v obou směrech:

```bash
swingersim ensemble configs/squadra_test3/qx_cas_ensemble.yaml -o output/qx_cas_ensemble -j 0
swingersim ensemble-report output/qx_cas_ensemble
```

Podmíněný ansámbl má $n_\mathrm{eff} = 186{,}7$; zúží vnější periodu z losovaných 118–4429 d
(16–84 %) na 193–1028 d, hmotnost a úhly ponechá téměř neomezené a předpovídá, že u 19,5 % soustav,
které nyní zakrývají, zákryty do 30 let skončí (obrázky 15–18).

## Vývoj za miliony let {#wf-evolution}

**Otázka.** Jak se hvězdy a dráhy mění během života soustavy: kdy hvězda vyplní svůj Rocheův lalok
a co následuje?

| Vstupy | Kroky | Výstupy |
| --- | --- | --- |
| hmotnosti, dráhy, stáří a metalicita soustavy; fyzika hvězd a jejich interakce | 1. zvolit model hvězd a stáří; 2. zkontrolovat hvězdy krátkým během; 3. zapnout vývoj, slapy, vítr, přenos hmoty; 4. spustit se sekulárním backendem dopředu nebo dozadu; 5. přečíst události a porovnat světelné křivky v různém stáří | běhy s vyvíjejícími se hvězdami (`_stars.png`), události historie přenosu hmoty, splynutí, soustava v pozdějším stáří |

**Tabulka 38.** Postup vývoje.

1. **Model hvězd.** `stellar_model: {kind: mist, age_myr: …, feh: …}` použije vývojové dráhy MIST
   od fáze před hlavní posloupností do konce asymptotické větve obrů (AGB); `kind: sse` rychlejší
   vzorce SSE až do zapálení helia. Všechny hvězdy jsou stejně staré; `stars.<m>.age_offset_myr`
   dá hvězdě zdánlivé stáří (obnažený dárce vypadá starší, omlazený příjemce mladší). Krátký běh
   (malé `t_max_yr`) vypíše v souhrnu $T_\mathrm{eff}$, poloměr, fázi a stáří každé hvězdy.
2. **Fyzika.** `physics.evolution: true` nechá hvězdy během běhu sledovat jejich vývojové dráhy;
   `tides: true` s `tidal_lag: hurley2002` bere slapovou disipaci z obálky každé hvězdy;
   `mass_loss` přidá hvězdný vítr; `mass_transfer` sleduje přetečení Rocheova laloku (stabilní
   přenos, kontaktní dvojhvězdy, společné obálky, přenos z třetí hvězdy); `mergers` pokračuje v běhu
   po splynutí dvou hvězd; `magnetic_braking` zpomaluje rotaci hvězd s konvektivními obálkami.
   Modely popisuje dokument *Fyzika a metody*, kapitoly Modely a vývoj hvězd, Slapy a rotace
   a Přenos hmoty, společná obálka a splynutí.
3. **Spusťte** s `run.backend: secular` (N-částicový backend nedokáže sledovat miliony let těsné
   dvojhvězdy); záporné `t_max_yr` integruje zpět k zrodu hvězd. Bez `mass_transfer` se běh zastaví,
   jakmile hvězda vyplní svůj Rocheův lalok; zastaví se také na konci vývojových drah a souhrn řekne
   proč (`EVENT`).
4. **Podívejte se na soustavu v pozdějším stáří**: konfigurace s pozdějším stářím a s drahami
   v tomto čase dá její světelnou křivku a okamžiky zákrytů.

Mějte na paměti soustavnou nejistotu vývojových drah hvězd; závisí-li závěr na stáří, losujte
`stellar_model.age_myr` v ansámblu. Řešený příklad je v kapitole [Řešené příklady](#ex-evolution).

## Stabilita, mapy chaosu a rozptyl {#wf-stability}

**Otázka.** Je soustava stabilní? Kde v prostoru jejích neznámých parametrů je stabilní nebo
chaotická? Mohla vzniknout zachycením?

| Vstupy | Kroky | Výstupy |
| --- | --- | --- |
| konfigurace; jeden nebo dva parametry a jejich rozsahy; pro rozptyl obě srážející se soustavy | 1. `stability` (sekundy); 2. dlouhý `run`; 3. `chaosmap` přes nejisté parametry; 4. `scatter` a `escape-theory` pro vznik | poměry a verdikty stability, mapy MEGNO s mezí stability, podíly výsledků a účinné průřezy |

**Tabulka 39.** Postup posouzení stability.

1. **Kritéria.** `swingersim stability` použije kritérium Mardlinga & Aarsetha (2001) a jeho
   zpřesnění podle Vynatheyi et al. (2022) na každou vnořenou trojici (dráhu a těleso, které ji
   ruší); poměr nad 1 znamená stabilitu. N-částicový test duchem integruje soustavu po řadu vnějších
   oběhů a porovnává ji s nepatrně posunutou kopií (Vynatheya et al. 2023). Kritéria jsou blízko 1
   jen přibližná; rozhodne dlouhý běh.
2. **Mapy chaosu.** `swingersim chaosmap` zintegruje mřížku konfigurací a spočítá MEGNO (2 pro
   pravidelný pohyb, roste pro chaos); úniky a srážky jsou označeny. Mapa ukazuje, kde soustava leží
   vůči nestabilním oblastem a rezonancím.
3. **Rozptyl.** `swingersim scatter` posílá dvojhvězdu proti jiné dvojhvězdě nebo osamocené hvězdě
   s řadou srážkových parametrů a zaznamenává výsledky; `escape-theory` dává pro srovnání
   statistickou předpověď podle Stonea & Leigha (2019).

**Řešený příklad.** Mapa chaosu cirkumbinární dráhy je v kapitole [Řešené příklady](#ex-chaos),
rozptylový experiment dvojhvězda–dvojhvězda v téže kapitole ([Řešené příklady](#ex-scatter)).

## Dlouhé výpočty na zpracovateli úloh {#wf-worker}

Ansámbly s $10^4$ členy, N-částicové ansámbly, proložení v kroku 2 s MCMC a mřížky MESA trvají hodiny
až dny. Běží jako **úlohy** na zpracovateli úloh ([Instalace](#install-worker)). Úloha je soubor YAML
v `jobs/inbox/`:

```yaml
command: ensemble                      # any sub-command except gui and worker
config: configs/squadra_test3/qx_cas_ensemble.yaml
output: output/qx_cas_ensemble
args: {n: 2000, processes: 0, resume: true}
note: QX Cas, unseen third body, 2000 members forward and backward
```

| Klíč | Význam |
| --- | --- |
| `command` | podpříkaz |
| `config` (a další poziční argumenty podle svého názvu: `ensemble`, `directory`, `h5`, `ids`, `model`) | poziční argumenty příkazu |
| `output` | výstupní argument (`-o`, nebo složka u `ensemble-report` a `view`) |
| `args` | volby podle svého dlouhého názvu, s `-` zapsaným jako `_` (`retry_errors: true`, `t_max_yr: 3000`) |
| `cores` | počet jader, která úloha obsadí (číslo nebo `all`); výchozí podle `processes` |
| `note` | volný text, zobrazený ve webové aplikaci |

**Tabulka 40.** Klíče souboru úlohy. Každá cesta musí ležet uvnitř složky SWINGERsimu; zpracovatel
úloh spouští jen příkazy SWINGERsimu s argumenty, které jejich parser povoluje.

Pravidla, na kterých v praxi záleží:

- Jednojádrové úlohy (`run`, `lightcurve`, `timings`, `minima`, `fit` v kroku 1) běží vedle sebe až
  do počtu jader; úloha, která potřebuje všechna jádra, běží sama. Každá úloha dostane počet jader,
  pro který byla zarezervována (proměnná prostředí `SWINGERSIM_CORES`), takže `processes: 0` znamená
  jádra úlohy, nikoli celý počítač.
- U ansámblů, rozptylu a map chaosu používejte `resume: true`: přerušená úloha pokračuje tam, kde se
  zastavila. Úloha přerušená zastavením zpracovatele nebo restartem počítače se automaticky vrátí do
  fronty. `retry_errors: true` znovu spustí jen členy, které skončily chybou.
- Soubor s jinou příponou (např. `.yaml.hold`) zůstane ve frontě, aniž by se spustil.
- Klesne-li volná paměť pod 8 % celkové (nejméně 1,5 GB), zpracovatel zastaví úlohu, která používá
  nejvíc paměti, a označí ji jako neúspěšnou s vysvětlením; její dosavadní výsledky zůstanou.
- Po každé aktualizaci SWINGERsimu zpracovatele restartujte.

## Reprodukovatelnost {#wf-reproducibility}

Každý výsledek lze dohledat až k jeho vstupu a znovu spočítat:

- soubor HDF5 běhu ukládá úplnou konfiguraci (se všemi doplněnými výchozími hodnotami), verze
  SWINGERsimu, Pythonu, NumPy, SciPy, h5py, REBOUND, REBOUNDx a emcee a čas svého vzniku; tabulky CSV
  a tabulky proložení nesou tytéž údaje o původu ve svých hlavičkových řádcích;
- složka ansámblu ukládá konfiguraci a její otisk; `--resume` odmítne smíchat řádky změněné
  konfigurace a vylosované členy určuje `ensemble.seed`;
- proložení zapíše proloženou konfiguraci s její epochou (`fit_epoch_bjd_tdb`), historii optimalizace
  a s MCMC i vzorky (`seed` určuje walkery);
- zpracovatel úloh uchovává každý soubor úlohy, jeho záznam a stav v `jobs/done/` nebo
  `jobs/failed/` a webová aplikace uchovává zmrazenou kopii každé odeslané konfigurace
  v `jobs/submitted/`.

Chcete-li výsledek sdílet, sdílejte konfiguraci, datové tabulky, na které odkazuje, a verzi
SWINGERsimu; se semínky v konfiguraci ho může kdokoli spočítat znovu. *Methods paragraph* na stránce
výsledku napíše popis výpočtu do článku.
