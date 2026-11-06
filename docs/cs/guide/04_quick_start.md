# Rychlý start {#quick-start}

Tato kapitola vás za několik minut dovede od nainstalovaného SWINGERsimu k prvním výsledkům, jednou
přes webovou aplikaci a jednou z příkazového řádku. Obě cesty počítají totéž: ukázkovou čtyřhvězdu
2+2 `configs/examples/example_2p2.yaml` za 200 let.

## Ukázková soustava {#qs-system}

```yaml
topology: "2+2"
masses: {m1: 1.2, m2: 1.0, m3: 0.9, m4: 0.8}        # solar masses
orbits:
  A:   {P: 6.74, e: 0.05, i: 85.0, Omega: 0.0,  omega: 90.0,  M0: 0.0}
  B:   {P: 2.68, e: 0.02, i: 60.0, Omega: 40.0, omega: 10.0,  M0: 120.0}
  out: {P: 500,  e: 0.3,  i: 70.0, Omega: 10.0, omega: 120.0, M0: 45.0}
physics: {gr: true}                                  # general relativity on
run: {t_max_yr: 200, n_output: 2000}                 # 200 years, 2000 snapshots
```

Čtyři hvězdy podobné Slunci tvoří dvě dvojhvězdy, A (perioda ve dnech: 6,74) a B (2,68), které kolem
sebe obíhají po vnější dráze `out` (500 d, excentricita 0,3). Každá dráha má periodu `P` ve dnech,
excentricitu `e` a tři úhly ve stupních: sklon `i` k rovině oblohy ($i = 90°$ znamená pohled z boku),
délku výstupného uzlu `Omega` a argument periastra `omega`; `M0` je střední anomálie na počátku, tj.
kde na dráze se hvězdy nacházejí v čase $t = 0$. Hvězdy nemají zadané poloměry, a proto dostanou
poloměry odpovídající hlavní posloupnosti nulového stáří.

## Cesta přes webovou aplikaci {#qs-gui}

```bash
swingersim gui
```

Příkaz spustí malý webový server a otevře v prohlížeči adresu `http://127.0.0.1:8765`.

1. Na stránce **Configuration** (konfigurace) vyberte v seznamu *Open…* (otevřít) soubor
   `configs/examples/example_2p2.yaml`. Formulář i text YAML ukazují soustavu; štítek vpravo nahoře
   hlásí *valid* (platná).
2. Napište název konceptu (např. `my_2p2`) a stiskněte **Save** (uložit). Konfigurace je nyní
   konceptem v `jobs/drafts/`.
3. V panelu **What to compute** (co počítat) zaškrtněte **Run** (běh) (a chcete-li, také
   **Eclipse times (O − C)**, okamžiky zákrytů) a stiskněte **Add and run** (přidat a spustit). Úlohy
   se zařadí do fronty a na pozadí se spustí zpracovatel úloh.
4. Stránka **Jobs** (úlohy) ukazuje běžící úlohu s ukazatelem průběhu. Asi po minutě se přesune do
   stavu *Finished* (dokončeno).
5. Na stránce **Outputs** (výstupy) otevřete kartu běhu. Stránka výsledků ukazuje modelovou soustavu,
   kdy která dvojhvězda zakrývá, dráhové elementy v čase a poměry stability.

Stránky jsou popsány panel po panelu v kapitole [Webová aplikace](#gui); obrázky 4, 7, 8 a 9 v ní
ukazují, co uvidíte.

## Cesta přes příkazový řádek {#qs-cli}

```bash
swingersim run configs/examples/example_2p2.yaml -o output/example.h5
```

Zhruba po minutě a půl program vypíše souhrn (běh této konfigurace v návodu trval 80 s na jednom
jádře). Nejdůležitější řádky jsou:

| Řádek souhrnu | Co vám říká |
| --- | --- |
| `out a = 1.94..1.9418 AU e = 0.2995..0.3007` | vnější dráha si zachovává velikost (asi dvojnásobek vzdálenosti Země–Slunce) i tvar |
| `i_mut A-B 25.59..46.67 deg` | úhel mezi rovinami A a B se kývá mezi 26° a 47°: dvojhvězdy se vůči sobě naklánějí |
| max \|dE/E\| = 3.23e-13 | energie se zachovává na 13 platných číslic: integrace je přesná |
| `hierarchy: configured [[1,2],[3,4]]`, `no change detected` | žádná hvězda nebyla vymrštěna ani nezměnila partnera |
| `precession period : 283 yr estimated` | rovina dvojhvězdy A se jednou otočí asi za 283 let |
| `stability ratio Mardling-Aarseth = 3.67 (stable)` | dráhy jsou od sebe dost daleko pro dlouhý život (nad 1 znamená stabilní) |
| `eclipses (observer along z): A 2 %, B 43 %` | B zakrývá během 43 % běhu, A jen krátce na začátku |
| `m1 R = 1.14 Rsun, …, Teff = 6144 K, main sequence` | hvězdy: poloměr, teplota, vývojová fáze |

**Tabulka 10.** Klíčové řádky souhrnu ukázkového běhu. Levý sloupec cituje řádek tak, jak je
vypsán, pravý sloupec říká, co znamená. Každý řádek souhrnu je vysvětlen v kapitole
[Výstupy](#out-summary).

Běh také zapsal vedle souboru HDF5 své standardní produkty:

```bash
ls output/example*
# example.h5  example.html  example_summary.txt  example_timeseries.csv
# example_orbits.png  example_inclinations.png  example_conservation.png  example_viewer.html
```

Otevřete ve webovém prohlížeči `output/example.html` (zprávu) a `output/example_viewer.html`
(animovaný 3D pohled na běh). Oba jsou samostatné soubory, které fungují bez SWINGERsimu a bez
připojení k internetu.

## Změňte jednu věc {#qs-change}

Typický způsob práce je změnit jeden parametr nebo zapnout či vypnout jeden efekt a porovnat
výsledky. Připravená varianta naklání dvojhvězdu B o 65° vůči vnější dráze:

```bash
swingersim run configs/examples/example_2p2_zlk.yaml -o output/zlk.h5
```

Excentricita B se nyní každých 14 let kývá mezi 0,01 a 0,81: to jsou cykly ZLK (von Zeipel–Lidov–
Kozai), výměna excentricity a sklonu dráhy rušené vzdáleným skloněným tělesem
([Řešené příklady](#ex-zlk)).

## Okamžiky zákrytů téže soustavy {#qs-timings}

```bash
swingersim timings configs/examples/example_2p2.yaml --years 20 -o output/example_timings.csv \
  --plot output/example_oc.png
```

Příkaz najde každý zákryt A a B v příštích 20 letech, včetně doby letu světla přes vnější dráhu,
a zapíše jejich O−C vůči lineární efemeridě (tabulku a obrázek). Pokračujte kapitolou
[Pracovní postupy](#workflows) pro úplné úlohy, nebo kapitolou [Řešené příklady](#examples) pro další
soustavy.
