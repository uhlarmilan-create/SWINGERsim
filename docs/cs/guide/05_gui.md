# Webová aplikace {#gui}

`swingersim gui` spustí na vašem počítači malý webový server a otevře v prohlížeči stránku, z níž
lze ovládat vše, co SWINGERsim umí: psát konfigurace, volit, co se má počítat, spouštět a řídit
úlohy a číst výsledky na stránkách připravených pro každý druh výsledku. Stránka je v angličtině.
Sama nic nepočítá: zapisuje soubory úloh, které provádí zpracovatel úloh (kapitola
[Instalace](#install-worker)), a čte soubory, které úlohy zapíší. Výsledek spočtený z příkazového
řádku proto vypadá stejně jako výsledek spuštěný ze stránky a prohlížeč lze během výpočtu úloh
zavřít.

Snímky obrazovky v této kapitole zachycují skutečnou relaci. Soustavou je V907 Sco, zákrytová
dvojhvězda s periodou 3,78 d a třetí hvězdou na dráze s periodou 142 d, jejíž zákryty přicházejí
a odcházejí (Zasche et al. 2023); její běhy, světelná křivka, okamžiky zákrytů a proložení jsou
tytéž jako v dokumentu *Testování, část III*. Ansámbl je ansámbl QX Cas z téhož dokumentu.

## Spuštění aplikace {#gui-start}

```bash
swingersim gui                                  # http://127.0.0.1:8765, opens the browser
swingersim gui --port 8800 --no-browser         # another port, no browser
swingersim gui --root ~/swingersim              # another SWINGERsim folder
swingersim gui --host 0.0.0.0                   # reachable from other computers
```

| Volba | Výchozí hodnota | Význam |
| --- | --- | --- |
| `--root` | aktuální složka, nebo nejbližší nadřazená složka, která obsahuje `configs/` | složka SWINGERsimu, jejíž `configs/`, `jobs/` a `output/` se používají |
| `--host` | `127.0.0.1` | adresa, na které server naslouchá; `0.0.0.0` ho zpřístupní ze sítě |
| `--port` | `8765` | port serveru |
| `--no-browser` | vypnuto | neotevírat prohlížeč |

**Tabulka 11.** Volby příkazu `swingersim gui`. S `--host 0.0.0.0` může kdokoli v síti číst
výstupy, psát konfigurace a spouštět úlohy; používejte to jen v důvěryhodné síti. Server poskytuje
jen soubory pod `--root`.

## Rámec každé stránky {#gui-frame}

Lišta nahoře nese název programu a čtyři karty – **Configuration** (konfigurace), **Jobs** (úlohy),
**Outputs** (výstupy) a **Tools** (nástroje) –, které představují čtyři stránky aplikace. Práce
postupuje jedním směrem: na stránce *Configuration* popíšete soustavu a zaškrtnete, co se má
počítat, *Jobs* provádí výpočty, *Outputs* ukazuje, co daly; *Tools* obsahuje kalkulačky, které
nepotřebují konfiguraci. Vpravo na liště otevírá **?** panel nápovědy k aktuální stránce a **◐**
přepíná mezi světlým a tmavým barevným schématem. Pod lištou začíná každá stránka modrým rámečkem,
který dvěma větami říká, k čemu stránka slouží; **×** ho zavře. Každý oddíl, pole, výpočet, graf
a varování má vlastní **?**, které otevře panel nápovědy s vysvětlením (obrázek 5).

## Configuration {#gui-config}

Stránka Configuration (konfigurace) slouží k psaní a úpravě konfiguračních souborů a k určení,
co se s nimi má počítat (obrázek 4).

![**Obrázek 4.** Stránka Configuration s otevřenou konfigurací V907 Sco. Vlevo: seznam konceptů
(zde prázdný). Nahoře: panel nástrojů – *New…*, *Open…*, volič topologie (`custom`), *Add stars and
orbits* a štítek stavu (*valid*, zelený). Pod ním lišta názvu: koncept se uloží jako
`jobs/drafts/v907_sco_published.yaml`; *Save as …*, *Save as new* a *Delete*; řádek pod ní uvádí,
odkud konfigurace pochází. Uprostřed: formulář, zde s oddíly *hierarchy* (sbalený), *masses* (tři
hvězdy) a *orbits* (dráha `inner` s poli periody, velké poloosy, excentricity a sklonu; každé pole má
**?**, jednotku a pod sebou výchozí hodnotu, **×** ho odstraní). Vpravo: text YAML téže konfigurace
(*preview*) a horní část panelu *Quick look*.](../../figures/guide/gui_config.png)

**Drafts** (koncepty). Seznam vlevo ukazuje konfigurace v `jobs/drafts/`. Kliknutím se koncept
otevře. Konfigurace odjinud (příklad v `configs/`, proložená konfigurace, člen ansámblu) se stane
konceptem, jakmile se uloží.

**Panel nástrojů.** *New…* (nová) založí konfiguraci, prázdnou nebo z jedné z pěti šablon (zákrytová
dvojhvězda se třetím tělesem, čtyřhvězda 2+2, čtyřhvězda 3+1, ansámbl pro neviděného průvodce
a proložení okamžiků minim). *Open…* (otevřít) načte libovolnou konfiguraci pod `configs/`, koncept
nebo odeslanou úlohu. Volič **topology** (topologie) nastavuje hierarchii (předvolby z obrázku 1
nebo `custom`); *Add stars and orbits* (přidat hvězdy a dráhy) doplní hmotnosti a dráhy, které
topologie potřebuje. Štítek vpravo ukazuje stav konfigurace: *valid* (platná), *not checked*
(nezkontrolovaná) nebo počet chyb.

**Lišta názvu.** Název konceptu (`jobs/drafts/<name>.yaml`). *Save* (uložit) zapíše koncept pod
tímto názvem (dokud se název liší od otevřeného souboru, tlačítko zní *Save as …*), *Save as new*
(uložit jako nový) se zeptá na nový název a nikdy nepřepíše existující koncept, *Delete* (smazat)
přesune koncept do `jobs/drafts/.deleted/`. Tečka vedle názvu označuje neuložené změny; úlohy vždy
používají uloženou verzi.

**Formulář.** Formulář se generuje ze schématu konfigurace, takže vždy nabízí každý klíč s jeho
jednotkou, výchozí hodnotou a jednořádkovým popisem; je to totéž schéma, jaké používá příkazový
řádek (kapitola [Konfigurační soubor](#configuration)). Nepoužité oddíly ukazují *not set*
(nenastaveno) a uvádějí, které výpočty je používají (např. *stars – used by Light curve, Eclipse
times, Tides*). Oddíl se otevře kliknutím; *new key* (nový klíč) a *Add* (přidat) přidají položku do
slovníku (hvězdu do `masses`, dráhu do `orbits`); **×** pole odstraní. Každá změna se okamžitě ověří
stejnými kontrolami jako na příkazovém řádku a chyby se zobrazí u pole a vypíší pod textem YAML.

**YAML.** Tatáž konfigurace jako text. *preview* (náhled) ji zobrazí; *edit as text* (upravit jako
text) ji změní na editor a *Apply text to the form* (převzít text do formuláře) načte upravený text
zpět do formuláře. Komentáře v textu formulář neuchovává; komentované soubory uchovávejte pod
`configs/`.

![**Obrázek 5.** Panel nápovědy otevřený tlačítkem **?** oddílu *orbits*. Zakrývá pravou třetinu
stránky a má části *What it is* (co to je), *When to use it* (kdy to použít), *What it affects* (co
to ovlivňuje), *Typical values* (typické hodnoty) a *Example* (příklad); *Insert the example into
the configuration* zkopíruje příklad do formuláře. Formulář vlevo zůstává použitelný; **×** panel
zavře.](../../figures/guide/gui_help.png)

Panel nápovědy (obrázek 5) je vestavěná příručka. U konfiguračního klíče vysvětluje, co klíč
znamená, kdy je potřeba, co ovlivňuje, typické hodnoty, příklad a případně literaturu; u výpočtu
říká, co se počítá, co konfigurace potřebuje, jak dlouho to trvá a co ukazuje stránka výsledků;
u grafu vysvětluje, jak ho číst; u varování říká, co znamená a co dělat.

![**Obrázek 6.** Pravý sloupec stránky Configuration pro V907 Sco, rozdělený na dvě části umístěné
vedle sebe. Levá část, shora: panel *Quick look* s hierarchií `((m1 + m2) + m3)`, poměry stability
jediné vnořené trojice (MA01 1,67, V22 2,20, „stable with a margin“, stabilní s rezervou)
a geometrií zákrytů páru `inner` (sklon 76,6° proti mezi zákrytů 72,8°: „eclipses at t = 0“,
zákryty v t = 0); poté *What to compute* se šesti otázkami *What do I want to know?* – vybrána je
otázka *What would a photometer see?* a její odpověď je v modrém rámečku – a položky *Run*
(nezaškrtnutá) a *Light curve* (zaškrtnutá, s volbami pásmo, začátek, délka, kadence a šum). Pravá
část: *Eclipse times (O − C)*, zaškrtnutá, přepnutá na *as an observer* se svými volbami (rozsah,
pásmo, kadence, šum; *shapes per season* a *also a PNG figure* jako doplňkové soubory); *Stability
check* a *Chaos map* (nezaškrtnuté); *Pulsation periods*, *Ensemble*, *Fit the measurements*
a *Scattering experiment* zašedlé s uvedeným důvodem (konfigurace nemá pulzující hvězdu ani oddíl
`ensemble`, `fit` či `scatter`); dole *Add to the inbox* a *Add and run* (nedostupné, dokud není
koncept uložen) a *From the command line*.](../../figures/guide/gui_config_compute.png)

**Quick look** (rychlý náhled; obrázek 6 vlevo nahoře) se počítá z konfigurace už během psaní,
dříve než se cokoli spustí: hierarchie v závorkovém zápisu, poměry stability každé vnořené trojice
podle Mardlinga & Aarsetha (2001, MA01) a Vynatheyi et al. (2022, V22) – nad 1 stabilní, od 1,5
zeleně – a pro každý pár dvou hvězd jeho sklon ve srovnání s mezí zákrytů
$\lvert\cos i\rvert < (R_1+R_2)/a$ a zaplnění Rocheových laloků jeho hvězd. Varuje včas, např. že
pozorovaná zákrytová dvojhvězda má při daném sklonu příliš malé poloměry na to, aby zakrývala.

**What to compute** (co počítat) uvádí všechny výpočty, které konfigurace umožňuje. Každá
zaškrtnutá položka se stane jednou úlohou uloženého konceptu s názvem `<number>_<draft>_<item>`; její
výsledky se objeví na stránce Outputs pod názvem konfigurace. Položky, které konfigurace neumožňuje,
jsou zašedlé a říkají, co chybí. Otázky pod *What do I want to know?* (co chci vědět) zaškrtnou
položky, které na ně odpovídají, a jednou větou řeknou jak (tabulka 12).

| Otázka | Zaškrtne položky | Poznámka |
| --- | --- | --- |
| Is it stable, and for how long? (Je stabilní, a jak dlouho?) | Stability check, Run | kontrola dá během minuty verdikt tří kritérií; dlouhý běh ukáže, co se skutečně stane |
| When will it eclipse? (Kdy bude zakrývat?) | Run, Eclipse times (exact) | běh ukáže geometrii zákrytů v čase, okamžiky zákrytů každý zákryt a jeho O−C |
| What would a photometer see? (Co by viděl fotometr?) | Light curve, Eclipse times (as an observer) | světelnou křivku a minima na ní změřená jako u skutečných dat |
| I have measurements – what system made them? (Mám měření – jaká soustava je vytvořila?) | Fit the measurements | otevře k vyplnění oddíl `fit` |
| Something unseen perturbs it – what can it be? (Něco neviděného ji ruší – co to může být?) | Ensemble | otevře k vyplnění oddíl `ensemble` |
| How did it look in the past? (Jak vypadala v minulosti?) | Run | se záporným `run.t_max_yr` |

**Tabulka 12.** Otázky panelu *What do I want to know?*. Každý řádek obsahuje otázku (vlevo),
výpočty, které vybere (uprostřed), a nápovědu, která se s ní zobrazí (vpravo).

| Položka | Příkaz | Volby v panelu | Potřebuje |
| --- | --- | --- | --- |
| Run | `run` | backend (prázdné: podle konfigurace) | hmotnosti a dráhy |
| Light curve | `lightcurve` | pásmo, začátek [d po $t = 0$], délka [d] (výchozí 27), kadence [min] (výchozí 2), šum [ppm]; obrázek PNG | poloměry a teploty (z `stars` nebo z modelu hvězd) |
| Eclipse times (O − C), exact | `timings` | rozsah [yr] (výchozí 10); obrázek PNG | poloměry pro mez zákrytů |
| Eclipse times (O − C), as an observer | `minima` | rozsah [yr] (výchozí 3), pásmo, kadence, šum; tvary po sezónách, obrázek PNG | jako světelná křivka |
| Stability check | `stability` | počet vnějších oběhů N-částicového testu (výchozí 100) | hmotnosti a dráhy |
| Chaos map | `chaosmap` | parametr x, rozsah x, počet bodů x (30); volitelně parametr y, rozsah, počet bodů | hmotnosti a dráhy |
| Pulsation periods | `pulsation` | – | hvězdu s `pulsation` |
| Ensemble | `ensemble` | počet členů (prázdné: `ensemble.n`) | oddíl `ensemble` |
| Fit the measurements | `fit` | krok (prázdné: 2, pokud existuje `fit.stage2`) | oddíl `fit` |
| Scattering experiment | `scatter` | – | konfiguraci rozptylu |

**Tabulka 13.** Položky panelu *What to compute*. Sloupce: popisek v panelu (Light curve – světelná
křivka, Eclipse times – okamžiky zákrytů, Stability check – kontrola stability, Chaos map – mapa
chaosu, Pulsation periods – periody pulzací, Fit the measurements – proložení měření, Scattering
experiment – rozptylový experiment), příkaz, který úloha spouští (kapitola
[Příkazový řádek](#command-line)), volby zobrazené po zaškrtnutí položky a co musí konfigurace
obsahovat.

Oba druhy okamžiků zákrytů odpovídají na různé otázky. *Exact* (přesné) bere okamžiky konjunkce
z integrace, včetně doby letu světla: skutečné O−C modelu, bez světelné křivky a bez šumu, za
sekundy výpočtu. *As an observer* (jako pozorovatel) spočítá syntetickou světelnou křivku a každé
minimum na ní změří fenomenologickým modelem Mikuláška (2015), s chybou: O−C stejného druhu jako ze
skutečných minim, za minuty. Liší se o to, co posouvá změřené minimum od konjunkce (nesouměrné
zákryty excentrických drah, světlo ostatních hvězd, pulzace, šum a vzorkování).

*Add to the inbox* (přidat do fronty) zapíše úlohy jako čekající (`.yaml.hold`), *Add and run*
(přidat a spustit) je zařadí do fronty a spustí zpracovatele úloh. *From the command line* (z
příkazového řádku) ukáže příkazy, jimž zaškrtnuté položky odpovídají, pro použití v terminálu nebo
ve skriptu.

## Jobs {#gui-jobs}

Stránka Jobs (úlohy) je ovládacím panelem zpracovatele úloh a fronty úloh (obrázek 7).

![**Obrázek 7.** Stránka Jobs během běhu úlohy. Nahoře: lišta zpracovatele – *worker running*
(zpracovatel běží, zeleně), počítač (`vm`), jeho 2 procesory, čas poslední kontroly, *cores in use:
1 of 2* – a *Stop worker*. Vlevo: *Drafts* s konceptem `v907_sco` (zaškrtávací políčko; *run* říká,
která položka by se spustila), *To the inbox* a *Run*. Vpravo, shora: *Inbox* se dvěma čekajícími
úlohami (sloupce úloha, příkaz, jádra, konfigurace, výstup, stav) a tlačítka *Select one-core jobs*,
*Run selected*, *Hold*, *Remove* a ↻; *Running* s úlohou `010_v907_sco_long_run` – N-částicový běh
na 3000 let – s ukazatelem průběhu (*t = 784.5 of 3000 yr · elapsed 1 min 1 s · about 2 min 52 s
left*), spouštěným příkazem a tlačítky *Pause*, *Stop* a *Intermediate results*; *Finished* s kartami
*Done (7)*, *Failed (0)*, *Cancelled (0)*, *Removed (0)* a jedním řádkem na každou dokončenou úlohu
s jejím příkazem, časem dokončení, dobou trvání a tlačítky *HTML report*, *Results* a *Run
again*.](../../figures/guide/gui_jobs.png)

**Lišta zpracovatele.** Zda zpracovatel úloh běží, na kterém počítači, s kolika jádry, kdy se
naposledy ohlásil a kolik jader používají běžící úlohy. *Start worker* (spustit zpracovatele) ho
spustí na pozadí, takže běží dál i po zavření prohlížeče nebo webového serveru; *Stop worker*
(zastavit zpracovatele) ho zastaví a běžící úloha se vrátí do fronty (ansámbl s `resume` pokračuje
tam, kde se zastavil). Zpracovatel spuštěný se starší verzí SWINGERsimu je označen: po aktualizaci ho
restartujte.

**Drafts** (koncepty). Zaškrtněte koncepty a stiskněte *To the inbox* (do fronty; čekající úlohy)
nebo *Run* (spustit; zařazené do fronty, zpracovatel spuštěn). Je to totéž jako *Add to the inbox*
a *Add and run* na stránce Configuration.

**Inbox** (fronta). Čekající (`.yaml.hold`) a zařazené (`.yaml`) úlohy. Sloupec *cores* (jádra)
říká, zda úloha potřebuje jedno jádro, nebo všechna. *Run selected* (spustit vybrané) spustí jen
zaškrtnuté úlohy: ostatní zařazené úlohy se pozdrží a stránka uvede které. *Select one-core jobs*
(vybrat jednojádrové úlohy) zaškrtne všechny úlohy, které potřebují jedno jádro, aby běžely vedle
sebe. *Hold* (pozdržet) pozdrží zaškrtnuté úlohy, *Remove* (odebrat) je přesune do `jobs/cancel/`,
↻ obnoví seznam.

**Running** (běžící). Jeden blok na každou běžící úlohu: její příkaz a jádra, ukazatel průběhu
s odhadem zbývajícího času (doba pozastavení se nezapočítává), příkazový řádek a tři tlačítka.
*Pause* (pozastavit) zastaví procesy úlohy (SIGSTOP) a *Continue* (pokračovat) je nechá běžet dál
přesně tam, kde byly. *Stop* (zastavit) úlohu ukončí, smaže její dosavadní výsledky (jen uvnitř
`output/`) a přesune úlohu do `jobs/cancel/`. *Intermediate results* (průběžné výsledky) otevře to,
co úloha dosud zapsala: dlouhý N-částicový běh zapisuje svůj stav každých 20 s, ansámbl zapisuje
jeden řádek za každého dokončeného člena.

**Finished** (dokončené). Hotové, neúspěšné, zrušené a odebrané úlohy na čtyřech kartách. *HTML
report* otevře zprávu úlohy, *Results* (výsledky) její stránku výsledků na stránce Outputs, *Run
again* (spustit znovu) vrátí kopii úlohy do fronty. Neúspěšná úloha ukazuje konec svého záznamu
(logu) a u běžných chyb (chybějící soubor, neplatná konfigurace, nedostatek paměti, mez chyby
energie) také to, co chyba znamená a co dělat.

Úlohy odeslané ze stránky, které mohou využít více jader (`ensemble`, `scatter`, `chaosmap`,
`fit` v kroku 2, příkazy MESA), dostanou `processes: 0` (všechna jádra úlohy), není-li nastavena
jiná hodnota. Jednotlivý běh, světelná křivka nebo sada okamžiků zákrytů je jeden výpočet na jednom
jádře.

## Outputs {#gui-outputs}

Každá dokončená, neúspěšná nebo běžící úloha je **výsledkem**. Stránka Outputs (výstupy) je uvádí
jako karty a každý otevře na stránce připravené pro jeho druh (obrázek 8).

![**Obrázek 8.** Seznam výsledků. Vlevo: *Configurations* – výsledky seskupené podle konfigurace,
z níž byly spočteny, s jejich počtem (*all configurations 7*, `example_fit` 1, `qx_cas_ensemble` 1,
`v907_sco_fit_lc_v2` 1, `v907_sco_published` 4), a *everything under output/*. Nahoře: nadpis
*Results*, vyhledávací pole a ↻. Pod ním filtrační štítky *All*, *Runs*, *Ensembles*, *Light curves*,
*Timings & minima*, *Fits*, *Benchmarks* a *Other*. Každá karta uvádí svůj druh (barevný štítek),
číslo úlohy a datum, svůj stav (*done*), konfiguraci, jeden řádek o tom, co se počítalo, malý graf
(O−C zákrytů, sklon vnitřní dvojhvězdy během běhu, světelnou křivku, rezidua proložení naměřené
minus modelové, podmíněné rozdělení parametru ansámblu) a dva až čtyři údaje přečtené z výstupu,
např. *Eclipsing (fraction of run) inner 50.5 %*. *Open* otevře stránku
výsledku.](../../figures/guide/gui_outputs.png)

Údaje na kartách se čtou z výstupních souborů; uvádějí, co výstup obsahuje (rozsahy, podíly,
události, poměry stability, efektivní počty členů), a soustavu nikdy neinterpretují. Barvy označují
kvalitu čísla jednotně na všech stránkách (tabulka 14).

| Veličina | zelená | jantarová | červená |
| --- | --- | --- | --- |
| největší relativní chyba energie $\lvert\Delta E/E\rvert$ běhu (u disipativních běhů se neoznačuje) | $< 10^{-6}$ | $< 10^{-3}$ | $\geq 10^{-3}$ |
| poměr stability (Mardling–Aarseth, Vynatheya et al. 2022) | $\geq 1{,}5$ | 1–1,5 | $< 1$ |
| efektivní počet členů $n_\mathrm{eff}$ ansámblu | $\geq 100$ | 20–100 | $< 20$ |
| redukované $\chi^2$, délka MCMC, apriorní rozdělení proložení | jak uvádí stránka proložení | | |

**Tabulka 14.** Barevný kód stránek výsledků. Každý řádek je veličina, sloupce jsou rozsahy
zobrazené zeleně (v pořádku), jantarově (prohlédněte si to) a červeně (bez kontroly číslu
nevěřte). Rozsahy jsou vypsány pod každým rámečkem *Data quality*.

Každá stránka výsledku má stejný rámec: záhlaví s názvem konfigurace, druhem výsledku, jeho stavem,
řádkem o tom, co se počítalo, a poznámkou úlohy; hlavní údaje v rámečcích; rámeček *Data quality*
(kvalita dat) vpravo; modrý rámeček *What am I looking at?* (na co se dívám); části připravené pro
tento druh výsledku; *Other results of this configuration* (další výsledky této konfigurace);
*Compute next* (spočítat dále); *For a paper* (do článku); a *Advanced* (pokročilé). *Compute next*
nabízí, co obvykle následuje (světelnou křivku nebo O−C běhu, další členy ansámblu, proloženou
konfiguraci proložení), a otevře stránku Configuration se zaškrtnutou položkou; zašedlé akce říkají,
proč nejsou k dispozici. *For a paper* vysvětluje, že každý obrázek lze uložit jako vektorovou
grafiku (ikona fotoaparátu nad ním) a jeho data jako CSV (ikona diskety), a *Methods paragraph*
napíše odstavec, který popisuje, jak byl výsledek spočten, i s potřebnými citacemi. *Advanced*
obsahuje soubory úlohy, její záznam, její konfiguraci a soubor úlohy, volné grafy libovolných
sloupců, všechny časové řady běhu a 3D prohlížeč.

### Stránka výsledku běhu {#gui-run}

![**Obrázek 9.** Horní část stránky výsledku běhu (V907 Sco, N-částicový běh na 200 let
s relativitou, slapy a rotací). Záhlaví: `v907_sco_published`, *run*, *finished*, „3 stars · orbits
out, inner · nbody · GR, tides, rotation · 0 → 200 yr · 4001 outputs · 39 s“. Rámečky údajů:
*Stopping event* none (ukončující událost: žádná), *Hierarchy changes* 0 (změny hierarchie),
*Eclipsing (fraction of run)* inner 50,5 % (zakrývá, podíl běhu), *Mutual inclination out–inner*
25,6–26,3° (vzájemný sklon), *Lowest stability ratio* 1,67 (MA01) (nejnižší poměr stability),
*Max |ΔE/E|* $1{,}3\times10^{-11}$. *Data quality*: zachování energie $1{,}3\times10^{-11}$, moment
hybnosti $4{,}1\times10^{-7}$, vývoj hvězd ne, backend nbody. *The model system* (modelová soustava):
souhrnné rámečky (hierarchie, dynamická stabilita, zákryty dráhy inner v $t = 0$); *Hierarchy* –
strom s drahami jako obdélníky a hvězdami jako kroužky obarvenými podle teploty; *On the sky* (na
obloze) – dráhy obou složek každé dráhy promítnuté na oblohu (x, y v au, uzel $\Omega$ měřený od x k y; modře: out,
oranžově: inner); *From above* (shora) – totéž při pohledu od pólu roviny oblohy, pozorovatel nahoře;
*Pair inner on the sky, stars to scale* (pár inner na obloze, hvězdy v měřítku) – obě hvězdy
zákrytového páru nakreslené se svými poloměry podél zdánlivé dráhy, se sklonem, mezí zákrytů, $a$
a srážkovým parametrem. Dole: tabulky drah ($P$, $e$, $i$) a vnořené trojice
($a_\mathrm{in}/a_\mathrm{out}$, vzájemný sklon, poměry MA01 a V22).](../../figures/guide/gui_run_top.png)

Horní část stránky běhu (obrázek 9) odpovídá na otázku „co je to za soustavu a podařil se výpočet?“.
Rámečky údajů říkají, zda běh skončil předčasně, zda se změnila hierarchie, které dvojhvězdy zakrývají
a po jakou část běhu, rozsah vzájemných sklonů, nejnižší poměr stability a chybu energie. *The model
system* kreslí konfiguraci v $t = 0$ čtyřmi způsoby. V pohledu na oblohu pár zakrývá, když se jeho
dvě hvězdy v konjunkci na obloze překrývají; ve zvětšeném panelu jsou hvězdy nakresleny v měřítku,
takže lze geometrii zákrytu přečíst přímo (zde je sklon 76,6° těsně nad mezí 72,85° a srážkový
parametr je 0,785 z $R_1+R_2$: tečné, částečné zákryty).

![**Obrázek 10.** Prostřední část stránky běhu. *System view* (pohled na soustavu; tmavý panel):
hvězdy v čase nastaveném posuvníkem, při pohledu ze Země (*as seen from Earth*) nebo kolmo na vnější
dráhu; vpravo detail páru `inner` s hvězdami v měřítku a popiskem *out of eclipse* (mimo zákryt)
nebo *eclipse now* (právě zákryt); *Play*, časový posuvník a rychlost přehrávání (*inner-orbit
speed* nebo *outer-orbit speed*). *Hierarchy*: nastavená hierarchie a zda se během běhu změnila.
*Stars* (hvězdy): hmotnost, poloměr, $T_\mathrm{eff}$ a největší zaplnění Rocheova laloku každé
hvězdy na začátku a na konci. *Interactive viewer* (interaktivní prohlížeč): otevře nebo uloží 3D
prohlížeč běhu jako jeden soubor HTML. *Eclipses over the run* (zákryty během běhu): sklon
dvojhvězdy `inner` k obloze (modrá křivka, stupně) v závislosti na čase (roky), její mez zákrytů
(čárkovaná čára, 72,85°) a časová okna, v nichž může zakrývat (modré pruhy nahoře); zde 50,5 % běhu
ve čtyřech oknech.](../../figures/guide/gui_run_system.png)

*System view* (obrázek 10 nahoře) běh animuje: polohy pocházejí z uložených oskulačních elementů,
mezi výstupy s keplerovským pohybem. *Eclipses over the run* je ústředním pozorovacím výsledkem:
sklon vnitřní dvojhvězdy osciluje mezi asi 53° a 89° s periodou precese uzlů (v tomto běhu změřeno
63 let) a dvojhvězda zakrývá jen tehdy, je-li sklon nad její mezí zákrytů
$\arccos[(R_1+R_2)/a]$. Pruhy nahoře jsou tato okna; mezera je doba, kdy se dráha odvrátila a zákryty
chybějí.

![**Obrázek 11.** Spodní část stránky běhu. *Orbits over the run* (dráhy během běhu): karty
*Eccentricity*, *Inclination to the sky*, *Mutual inclinations*, *Semi-major axis* a *Energy error*
(excentricita, sklon k obloze, vzájemné sklony, velká poloosa, chyba energie); jeden malý panel na
každou dráhu, každý s vlastním svislým rozsahem (min–max vypsané vpravo nahoře); zde excentricita
dráhy `inner` (0,0041–0,0131) a dráhy `out` (0,4663–0,4712) v závislosti na čase v letech.
*Secular dynamics: theory vs. the run* (sekulární dynamika: teorie a běh): pro každou vnitřní dráhu
vzájemný sklon v $t = 0$, rozsah excentricity podle kvadrupólové teorie a podle běhu, odhadnutá
a změřená perioda precese, změřená perioda apsidálního stáčení a perioda očekávaná pro samotnou
dvojhvězdu a poměry MA01 a V22. *Other results of this configuration* (okamžiky zákrytů, kontrola
stability, světelná křivka) a *Compute next* (Light curve, Eclipse times přesné nebo jako
pozorovatel, Stability check; *Ensemble around it* a *Pulsation periods* zašedlé s uvedeným
důvodem; *Change and run again*).](../../figures/guide/gui_run_orbits.png)

*Orbits over the run* (obrázek 11) ukazuje oskulační elementy vzorkované ve výstupech; rychlé
oscilace vnější excentricity (oranžově) jsou odezvou na vnitřní dvojhvězdu při každém vnějším oběhu,
nikoli sekulární změnou. Tabulka *Secular dynamics* srovnává běh s kvadrupólovou teorií středovanou
přes oběh: shoda říká, že se soustava chová tak, jak předpovídá analytická teorie; velké rozdíly
ukazují na oktupólové, rezonanční nebo nehierarchické efekty (*Fyzika a metody*, kapitola
Dynamika).

### Světelné křivky {#gui-lightcurve}

![**Obrázek 12.** Stránka světelné křivky (V907 Sco, pásmo TESS, 30 d s kadencí 10 min). *Light
curve*: relativní tok (1 = všechny hvězdy mimo zákryt) v závislosti na čase ve dnech po $t = 0$,
s tlačítky *show:* pro celkové světlo, dvojhvězdu *inner (m1 + m2)* a každou hvězdu a s volbou
*relative flux* (relativní tok) nebo *magnitude* (magnituda); lišta pod grafem vybírá časový
rozsah. *Folded on the period of a binary* (složeno s periodou dvojhvězdy): světlo dráhy `inner`
složené s její periodou (3,77487 d, perioda z konfigurace; fáze 0 = její nejhlubší bod), s volbou *light of this binary*
(světlo této dvojhvězdy) nebo *total light* (celkové světlo). *Light of each star and binary*
(světlo každé hvězdy a dvojhvězdy): podíl na světle mimo zákryt a největší ztráta vlastního světla
(inner 98,3 % a 6,18 %; m1 46,5 % a 12,1 %; m2 52,1 % a 10,0 %; m3 1,69 % a 0,0413 %). *Compare with
a measured light curve* (srovnání s naměřenou světelnou křivkou): soubor CSV (čas v BJD nebo ve
dnech, tok nebo magnituda) vykreslený přes model, s BJD modelového času
$t = 0$.](../../figures/guide/gui_lightcurve.png)

Stránka světelné křivky (obrázek 12) rozkládá světlo na jednotlivé části, takže každý rys lze
přiřadit hvězdě, která ho způsobuje. Složená křivka ukazuje primární a sekundární zákryty dvojhvězdy
s její vlastní periodou; rozptyl kolem středního tvaru pochází od ostatních hvězd a ze změn dráhy
během okna. Naměřená data se vykreslí přes model po zadání BJD modelového času $t = 0$, protože model
nemá kalendář: jeho čas začíná v $t = 0$, v epoše elementů.

### Okamžiky zákrytů a minima {#gui-timings}

![**Obrázek 13.** Stránka okamžiků zákrytů (přesné okamžiky z dynamiky, V907 Sco, 100 let před
a po $t = 0$). *O − C of every eclipsing binary* (O − C každé zákrytové dvojhvězdy): O−C v minutách
vůči lineární efemeridě proložené primárními zákryty, v závislosti na čase v letech od $t = 0$; plné
body: primární zákryty (zákryty m1), prázdné kroužky: sekundární zákryty, posunuté o svůj střední
posun (−8,44 min, uvedený v legendě). Zákryty existují jen ve čtyřech oknech; v každém okně O−C
klesá zhruba o 200 min. *Per binary* (po dvojhvězdách): počty primárních a sekundárních zákrytů
(8953, 8924), perioda efemeridy (3,776279318 d), rozsah O−C (−104,2 až 103,8 min), jeho polovina
a srážkový parametr. *Compare with measured minima* (srovnání s naměřenými minimy): soubor CSV
naměřených minim a BJD modelového času $t = 0$. *Other results* a *Compute
next*.](../../figures/guide/gui_timings.png)

Obrázek 13 ukazuje, jak se dynamika projevuje v diagramu O−C. Lineární efemerida proložená všemi
zákryty nedokáže sledovat pomalý drift periody zákrytů uvnitř každého okna, který způsobuje precese
dráhy vůči zornému paprsku; uvnitř okna probíhají primární a sekundární zákryty téměř souběžně, takže
apsidální stáčení téměř kruhové dvojhvězdy je malé. U dvojhvězdy s dráhou projevující se efektem
světelného času ukazuje O−C naopak vlnu s vnější periodou. Stránka minim změřených *as an observer*
má stejné uspořádání a přidává chybu každého minima, tvar zákrytů (hloubku, šířku a špičatost
Mikuláškova profilu, po sezónách) a posun mezi změřeným minimem a geometrickou konjunkcí.

### Kontroly stability a další výsledky {#gui-other}

![**Obrázek 14.** Stránka kontroly stability, příklad stránek příkazů bez vlastního uspořádání.
Rámečky údajů: *inner in out: stability ratio* MA01 1,67 · V22 2,2 („with m3, i_mut 26.2° · above 1:
stable“) a *N-body ghost test* (N-částicový test duchem) stable („100 outer orbits, largest
divergence $9.60\times10^{-4}$“). *Report* (zpráva): text, který úloha vypsala – vnořená trojice
s $a_\mathrm{in}/a_\mathrm{out}$, excentricitami, vzájemným sklonem, vnějším poměrem hmotností
a oběma poměry stability – a výsledek N-částicového testu. *Other results of this configuration*,
*Compute next* (*Open the configuration*, *Tools*), *For a paper* a
*Advanced*.](../../figures/guide/gui_stability.png)

Příkazy bez vlastní stránky (kontroly stability, mapy chaosu, rozptylové experimenty, periody
pulzací, emulátory) ukazují text, který vypsaly, a údaje v něm rozpoznané (obrázek 14), se svými
soubory pod *Advanced*.

### Stránka výsledku ansámblu {#gui-ensemble}

![**Obrázek 15.** Horní část stránky ansámblu (QX Cas: 2000 členů, dopředu i dozadu, každý
N-částicově na 30 let; losována hmotnost a dráha třetí hvězdy). Záhlaví: 3 hvězdy, backend, délka
běhu, „7 varied parameters · 4 observations · required: stable_ma01“. Rámečky údajů: *Members*
(členy) 2 000, *Effective members (n_eff)* (efektivní počet členů) 186,7 („conditioned on 4
observations“, podmíněno 4 pozorováními), *Most common fate* (nejčastější osud) intact 100 %,
*Stable (Mardling–Aarseth)* 100 %, *Eclipsing (time average)* inner 94 %. *Data quality*:
$n_\mathrm{eff}$ s pozorováními (zeleně), backend, členy s chybou
(0).](../../figures/guide/gui_ensemble_top.png)

Stránka ansámblu (obrázky 15–18) odpovídá na otázku „které soustavy jsou možné a které připouštějí
měření?“. Modelová soustava nahoře (není zobrazena) je nominální konfigurace.

![**Obrázek 16.** Podmínění a parametry ansámblu. *Observations used for conditioning* (pozorování
použitá k podmínění): jeden řádek na každou položku `ensemble.observed` se zaškrtávacím políčkem
*use* (použít), veličinou (`inner_i0`, `inner_didt0_deg_yr`, `fate` dozadu, `fate` dopředu), směrem
(dopředu, centrálně, dozadu) a podmínkou (82 ± 1,5, −0,23 ± 0,08, = intact); efektivní počet členů se
zaškrtnutými pozorováními (186,7) a *Re-weight* (převážit). *Parameters: drawn vs. conditioned*
(parametry: vylosované a podmíněné): jeden histogram na každý losovaný parametr a pro počáteční
vzájemný sklon, šedě = vylosované (apriorní rozdělení), modře = vážené pozorováními; nad každým medián
a rozsah 16–84 % před → po; panely jsou seřazeny podle toho, jak moc pozorování rozsah zúžila (první
je vnější perioda, její rozsah ×0,46; logaritmická osa u log-rovnoměrného
losování).](../../figures/guide/gui_ensemble_params.png)

Tabulka pozorování (obrázek 16) je interaktivní: zrušte zaškrtnutí měření a stiskněte *Re-weight*,
abyste viděli, čím přispívá; nic se znovu neintegruje. Váha člena je $w = \exp(-\chi^2/2)$ pro
gaussovská měření (1 pro přesnou shodu, 0,61 při 1σ, 0,14 při 2σ) a 0, když člen porušuje mez;
efektivní počet členů $n_\mathrm{eff} = (\sum w)^2/\sum w^2$ (Kish 1965) říká, kolik členů nese
váhu. Pod zhruba 20 jsou podmíněné histogramy jen šumem a jsou skryty. Kde je modrý histogram užší
než šedý, měření parametr omezují: zde dávají přednost vnějším periodám 190–1030 d před losovanými
120–4400 d a úhly nechávají téměř neomezené.

![**Obrázek 17.** Výsledky ansámblu vážené pozorováními. Karty *Forward (future)* (dopředu,
budoucnost) a *Backward (past)* (dozadu, minulost). *Fates* (osudy): intact 100 % (99,5–100 %, 68%
Wilsonův interval). *Stability*: podíly stabilních podle Mardlinga–Aarsetha a Vynatheyi et al.
(2022). *Orbits*: pro každou dráhu největší excentricita (medián a 16–84 %), podíl s $1-e < 0{,}1$,
současná rychlost změny sklonu $di/dt$ v $t = 0$, časově středovaný podíl zákrytů a podíl členů,
jejichž zákryty během běhu skončí (19,5 %, medián 23 let). *Mutual inclinations*: v $t = 0$ (84,7°,
44–135°) a podíl členů, které překročí 90°. Histogramy (šedě vylosované, modře vážené): největší $e$
dráhy inner a dráhy out, podíl zákrytů dráhy inner, nejnižší poměr stability
MA01.](../../figures/guide/gui_ensemble_outcomes.png)

Výsledky (obrázek 17) jsou odpovědí ansámblu: podíly s jejich 68% Wilsonovými intervaly (Wilson
1927), které se chovají dobře i pro podíly blízké 0 nebo 1, a mediány s rozsahy 16–84 %, vše vážené.
Karta *Backward* dává totéž pro běhy do minulosti.

![**Obrázek 18.** *Members with the largest weight* (členy s největší vahou). Šedý rámeček tabulku
vysvětluje: každý řádek je jedna vylosovaná soustava; *weight* (váha) říká, jak dobře reprodukuje
zaškrtnutá pozorování; *share* (podíl) je její část celkové váhy; sloupce stability jsou poměry MA01
a V22 (zeleně ≥ 1,5, jantarově 1–1,5); „59 members have a weight above 0.5: they reproduce the
observations equally well“ (59 členů má váhu nad 0,5: reprodukují pozorování stejně dobře).
Sloupce: číslo člena, váha, $\chi^2$, podíl, vylosované hodnoty ($m_3$, vnější $P$, $e$, $i$,
$\Omega$, $\omega$, $M_0$), osud dopředu/dozadu, oba poměry stability a *open as
configuration*.](../../figures/guide/gui_ensemble_members.png)

*Open as configuration* (otevřít jako konfiguraci; obrázek 18) promění člena v obyčejný koncept s jeho
vylosovanými hodnotami, připravený ke spuštění, se světelnou křivkou nebo diagramem O−C. Odpovědí
ansámblu je rozsah parametrů, nikoli jeden člen: členy s podobnými vahami jsou stejně dobrými
příklady.

### Stránka výsledku proložení {#gui-fit}

![**Obrázek 19.** Horní část stránky proložení v kroku 2 (V907 Sco: 58 minim, 172 radiálních
rychlostí, 2562 bodů světelné křivky, 14 sklonů, 3 apriorní rozdělení, 16 volných parametrů, 4384
vyhodnocení modelu). Rámečky údajů: *χ² start → fitted* (χ² počáteční → proložené) 11 920 → 1052
(„reduced 0.3765 – below 0.5: the errors used are larger than the scatter, the parameter errors are
conservative“: redukované 0,3765 – pod 0,5: použité chyby jsou větší než rozptyl, chyby parametrů
jsou konzervativní), *priors > 2σ off* (apriorní rozdělení odchýlená o více než 2σ) 1
(`stars.m2.R` −4,3σ), *inner: rms of the minima* 2,6 min (na začátku 3,5 min), *radial velocities:
rms* 3,85 km/s, *light curves: rms* 0,666 ppt (2 segmenty), *inclinations: χ²* 222 (14 epoch). *Data
quality*: redukované χ² a dvě varování (faktory chyb se po 4 proloženích stále měnily; silně
korelovaná dvojice parametrů). *The model system*: jako u běhu, s počáteční konfigurací čárkovaně
a proloženou plně.](../../figures/guide/gui_fit_top.png)

Stránka proložení (obrázky 19–23) ukazuje model proložený měřeními a vše, co je potřeba k jeho
posouzení. Nejprve čtěte varování: říkají, co vyžaduje pozornost.

![**Obrázek 20.** Varování, parametry a O−C proložení. *Check these first* (nejprve zkontrolujte):
varování proložení s jejich **?**. *Parameters* (parametry): jeden řádek na každý volný parametr s jeho
počáteční hodnotou, proloženou hodnotou, lineární chybou z kovarianční matice, jednotkou, apriorním
rozdělením (střední hodnota ± σ) a odchylkou proložené hodnoty od apriorního rozdělení a s posunem
(proložená − počáteční)/chyba; odchylky od apriorního rozdělení nad 2σ jsou červeně. *Eclipse timing
(O − C) with the model* (okamžiky zákrytů, O − C, s modelem): naměřená minima (body: primární plné,
sekundární prázdné, s chybovými úsečkami) a každý zákryt proloženého modelu (plně: primární,
čárkovaně: sekundární), jako O−C v minutách vůči jedné lineární efemeridě za roky
2021–2025.](../../figures/guide/gui_fit_params.png)

V panelu O−C na obrázku 20 apsidální stáčení od sebe vzdaluje křivky primárních a sekundárních minim,
zatímco efekt světelného času a poruchy od třetí hvězdy posouvají obě společně; zde leží sekundární
křivka asi 15–25 min pod primární a obě nesou 142denní vlnu vnější dráhy.

![**Obrázek 21.** Rezidua minim a radiálních rychlostí. *Times of minima* (okamžiky minim): naměřené
minus proložený model v minutách v závislosti na čase, primární plně, sekundární prázdně, křížky pro
minima vynechaná jako odlehlá; tlačítka přepínají mezi proloženým a počátečním modelem; chybové
úsečky jsou chyby, které proložení použilo (s minimální chybou zdroje). Tabulky: rms podle zdrojů
(sektory TESS) a střední reziduum primárních a sekundárních minim s jejich rozdílem (0,435 min).
*Radial velocities* (radiální rychlosti): naměřené body a model ve stejných časech pro m1 a m2,
v závislosti na fázi dvojhvězdy (karty *phase of its own binary*, *against time*, *phase of inner*,
*phase of out*); modře CHIRON 2022, oranžově Lacy (1999); šedé značky jsou model, který zahrnuje pohyb
kolem třetí hvězdy, a proto není jednou sinusoidou; dole rezidua v km/s.](../../figures/guide/gui_fit_minima_rv.png)

Rozdíl mezi středními rezidui primárních a sekundárních minim (obrázek 21) je citlivou kontrolou:
rozdíl několika σ znamená, že model klade jeden druh zákrytů do špatné fáze, obvykle kvůli
chybnému $e\cos\omega$ páru.

![**Obrázek 22.** Světelné křivky a sklony proložení. *Light curves* (světelné křivky): jedna karta
na každý segment (TESS S39, S66), *fitted* (proložený) nebo *+ start model* (i počáteční model);
naměřený tok (body) a kalibrovaný model (čára; měřítko, nulový bod a trend řešené pro každý segment)
v oknech kolem zákrytů, v závislosti na BJD − 2459361; dole rezidua v promilích (ppt); tabulka uvádí
počet bodů a rms každého segmentu. *Inclinations against time* (sklony v čase): sklon páru
v N-částicovém modelu (modře: proložený, oranžově: počáteční) a naměřené epochy (černé body
s chybovými úsečkami) od roku 1900 do 2022; tabulka uvádí každou epochu s naměřenou hodnotou, její
chybou, počáteční a proloženou modelovou hodnotou, mezí zákrytů (pro jednostranná data),
normovaným reziduem a zdrojem (DASCH, APPLAUSE, TESS, FRAM).](../../figures/guide/gui_fit_lc_inc.png)

Panel sklonů (obrázek 22) ukazuje vlastnost, kterou dokáže proložit jen dynamický model: zákryty
V907 Sco během 20. století zmizely a znovu se objevily, protože se mění sklon její dráhy,
a N-částicový model celé soustavy naměřené epochy sleduje.

![**Obrázek 23.** Konec stránky proložení. *Predicted minima* (předpovězená minima; nahoře,
oříznuto): příští minima v BJD_TDB a jako kalendářní datum, s nejistotou σ v minutách, dobou trvání
zákrytu a srážkovým parametrem. *Other fits of this system* (další proložení této soustavy): každá
složka proložení vedle této s redukovaným χ², počtem volných parametrů, rms minim a použitými daty;
*Compare the parameters with* (porovnat parametry s) přidá do tabulky parametrů druhé proložení.
*Report*: textová zpráva proložení. *Compute next*: *Open the fitted configuration* (otevřít
proloženou konfiguraci; běh, světelná křivka, okamžiky zákrytů řešení), *Change the fit and run
again* (změnit proložení a spustit znovu), *HTML report* (zde zašedlé) a *Could a residual wave be
magnetic activity?* (může být zbytková vlna magnetickou aktivitou?; Applegateův test). *For a paper*
a *Advanced*.](../../figures/guide/gui_fit_bottom.png)

Předpovědi (obrázek 23) se počítají z mnoha sad parametrů vylosovaných ze vzorků MCMC nebo
z kovarianční matice, s jednou N-částicovou integrací pro každou, takže jejich chyby zahrnují
korelace parametrů. *Other fits of this system* umožňuje na první pohled porovnat varianty (se
světelnými křivkami a bez nich, s apriorními rozděleními nebo radiálními rychlostmi a bez nich).

## Tools {#gui-tools}

![**Obrázek 24.** Stránka Tools (nástroje) s vybranou kalkulačkou Applegateova mechanismu. Vlevo:
nástroje seskupené jako *Calculators* (kalkulačky: Applegateův mechanismus, teorie rozpadu soustavy
tří těles), *Checks* (kontroly: ověření proti publikovaným výsledkům), *Results* (výsledky:
interaktivní prohlížeče běhů) a *Stellar models* (modely hvězd: vlastní modely programem MESA); u
každého příkaz, který spouští. Vpravo: nápověda k nástroji (co dělá, kdy ho použít, co potřebuje, co
dostanete, literatura), jeho formulář (hmotnost, poloměr, $T_\mathrm{eff}$, $a$, $P_\mathrm{bin}$,
$P_\mathrm{mod}$ aktivní hvězdy a dvojhvězdy; pozorovaná poloviční amplituda O−C), název úlohy
a tlačítka *Run now*, *Add to the inbox*, *Add and run*. Dole: *Where the former commands are* (kde
jsou dřívější příkazy) – pro každý příkaz stránka a položka, kde se
nachází.](../../figures/guide/gui_tools.png)

Stránka Tools (obrázek 24) obsahuje to, co nepatří ani ke konfiguraci, ani k výsledku. Každý nástroj
má formulář vygenerovaný z voleb svého příkazu (kapitola [Příkazový řádek](#command-line)). *Run
now* (spustit hned) spustí příkaz okamžitě a ukáže, co vypíše (nejvýše 5 minut); delší výpočty jdou
do fronty. Tabulka dole uvádí, kde se který příkaz v aplikaci nachází (tabulka 15).

| Příkaz | Kde v aplikaci |
| --- | --- |
| `run`, `lightcurve`, `timings`, `minima`, `stability`, `chaosmap`, `ensemble`, `pulsation`, `scatter`, `fit` | Configuration › What to compute |
| `view` | Outputs › běh › Interactive viewer; všechny běhy složky: Tools › Interactive viewers |
| `ensemble-report` | Outputs › ansámbl (textová zpráva pod Advanced) |
| `emulator-train`, `sbi` | Outputs › ansámbl › Compute next |
| `emulator-report` | Outputs › emulátor › Compute next |
| `validate`, `applegate`, `escape-theory`, `mesa-grid`, `mesa-reduce` | Tools |
| `rv`, `apsidal`, `worker`, `gui` | jen příkazový řádek (zpracovatel úloh se spouští na stránce Jobs) |

**Tabulka 15.** Kde se ve webové aplikaci nachází každý příkaz SWINGERsimu. Vlevo: podpříkaz;
vpravo: stránka a místo na ní.

## Složky používané aplikací {#gui-folders}

| Složka | Obsah | Zapisuje |
| --- | --- | --- |
| `configs/` | konfigurace dodávané se SWINGERsimem a vaše vlastní komentované soubory | vy |
| `jobs/drafts/` | konfigurace napsané v aplikaci, členy ansámblů otevřené jako konfigurace | Configuration, Outputs |
| `jobs/submitted/` | zmrazené kopie konceptů v podobě, v jaké byly odeslány | Configuration, Jobs |
| `jobs/inbox/`, `running/`, `done/`, `failed/`, `cancel/` | fronta úloh (tabulka 9) | aplikace a zpracovatel úloh |
| `output/<job>/` | výsledky úlohy odeslané z aplikace | příkaz úlohy |
| `output/.gui_cache/` | uložené karty, stránky a prohlížeče; znovu se sestaví, když se změní úloha nebo její výstup | aplikace |

**Tabulka 16.** Složky, které webová aplikace čte a zapisuje, všechny uvnitř složky SWINGERsimu
zadané volbou `--root`.
