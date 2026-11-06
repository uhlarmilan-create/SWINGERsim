# Řešení problémů a časté otázky {#troubleshooting}

Nejčastější problémy, seřazené podle toho, kde se objevují. Každá tabulka uvádí příznak, jeho
obvyklou příčinu a nápravu. Na obecné otázky odpovídá konec kapitoly.

## Konfigurace je odmítnuta {#ts-config}

SWINGERsim před spuštěním zkontroluje celou konfiguraci a zastaví se se zprávou, která jmenuje
příslušný klíč. Ve webové aplikaci se stejná zpráva zobrazí u pole.

| Zpráva (zkrácená) | Příčina | Náprava |
| --- | --- | --- |
| `give exactly one of P (days) or a (AU)` | dráha má zadané obě veličiny, nebo žádnou | ponechte jednu |
| `hierarchy: unknown child 'x'` | název v `hierarchy` není ani hvězda (`m1`, …), ani dráha | zkontrolujte pravopis; hvězdy se číslují od 1 bez mezer |
| `hierarchy: needs exactly one outermost orbit` | dvě dráhy nejsou spojeny | přidejte dráhu, která je spojuje |
| `topology … needs orbits [...]`, `masses must be exactly [...]` | dráhy nebo hmotnosti neodpovídají hierarchii | ve webové aplikaci použijte *Add stars and orbits*, nebo porovnejte s tabulkou 46 |
| `stellar_model.kind = mist needs 'age_myr'` | vyvíjející se modely potřebují stáří | zadejte `stellar_model.age_myr` |
| `physics.rotation needs physics.tides` | rotace a slapy sdílejí jeden model deformace hvězdy | přidejte `tides: true` |
| `physics.evolution: stars [...] have a fixed radius R or Teff` | při vývoji plynou poloměry a teploty z hmotnosti a stáří | odstraňte `R` a `Teff`, nebo vývoj vypněte |
| `physics.mass_transfer needs physics.evolution: true and stellar_model.kind sse or mist` | historie přenosu hmoty pocházejí z vývojových drah | zapněte vývoj se `sse` nebo `mist` |
| `physics.gas needs run.backend: secular` | plynná fáze je implementována v sekulárním backendu | nastavte `run.backend: secular` |
| `run.megno assumes constant masses` | MEGNO se ztrátou nebo přenosem hmoty | jedno z nich vypněte |
| `give value and sigma, or min and/or max (or equals)` | položka `observed` je neúplná | pozorování je hodnota s chybou, mez, nebo textová hodnota |
| `fit.stage2.priors.<name>: not a free parameter` | apriorní rozdělení pevného parametru nemá žádný účinek | uvolněte parametr, nebo apriorní rozdělení odstraňte |
| `fit.lc red_noise: true … needs fit.stage2.error_scale … or clip_sigma` | faktory červeného šumu se uplatňují mezi proloženími v řadě | přidejte `error_scale` nebo `clip_sigma` |
| `Extra inputs are not permitted` | překlep v klíči; u klíčů záleží na velikosti písmen | porovnejte s kapitolou [Konfigurační soubor](#configuration) |

**Tabulka 76.** Běžné zprávy při odmítnutí konfigurace.

## Běh je pomalý nebo skončí předčasně {#ts-run}

| Příznak | Příčina | Náprava |
| --- | --- | --- |
| N-částicový běh trvá hodiny | každý oběh nejkratší dráhy potřebuje desítky kroků | zůstává-li hierarchie neporušená, použijte sekulární backend; zkraťte `t_max_yr`; pro sekulární backend nainstalujte Numba (`fast`) |
| brzy `EVENT: collision` | dvě hvězdy se dotknou v periastru | u vysokých excentricit očekávané; zkontrolujte poloměry a zda slapy a relativita potlačují cykly ZLK |
| `EVENT: roche_overflow` | hvězda vyplní svůj Rocheův lalok | chcete-li ho sledovat, zapněte `physics.mass_transfer`, nebo berte čas události jako odpověď |
| `EVENT: stellar_evolution_end` | hvězda dosáhla konce své vývojové dráhy | běh nemůže pokračovat za hranice modelů |
| `EVENT: wall_time` | dosaženo `run.max_wall_s` | zvyšte hodnotu, nebo použijte sekulární backend |
| `EVENT: energy_error` | $\lvert\Delta E/E\rvert$ překročilo `run.max_energy_error`, obvykle při velmi těsném přiblížení | hledejte téměř srážky; ověřte nástroji pro stabilitu |
| sekulární backend odmítne běh se ztrátou hmoty | index ztráty hmoty Ψ překračuje `max_psi`: adiabatická odezva středovaná přes oběh neplatí | použijte N-částicový backend |

**Tabulka 77.** Problémy během běhu.

## Výsledky vypadají chybně {#ts-results}

| Příznak | Příčina | Náprava |
| --- | --- | --- |
| žádné zákryty, ačkoli je dvojhvězda pozorována jako zákrytová | výchozí poloměry hlavní posloupnosti nulového stáří jsou pro daný sklon příliš malé | zadejte naměřené `stars.<m>.R` nebo vyvinutý model hvězd; Quick look ukazuje mez zákrytů |
| apsidální stáčení je příliš rychlé | polytropické $k_2$ modelu `zams` je pro hvězdy 2–3 M☉ příliš velké | použijte `mist` nebo `sse`, nebo zadejte naměřené `apsidal_k2` |
| cykly ZLK dosahují $e \approx 1$ a hvězdy se srazí | hmotné body: žádná slapová ani relativistická precese | zapněte `gr`, `tides` a `rotation` |
| oskulační elementy kmitají | jsou to okamžité elipsy, rušené při každém oběhu | normální; sekulární backend uvádí středované hodnoty |
| rychlost změny se liší od publikované | jiné časové okno nebo definice | nastavte `run.rate_window_yr` zhruba na polovinu pozorovacího intervalu |
| primární a sekundární zákryty jsou oproti článku prohozené | primární zákryt je ve SWINGERsimu zákryt první hvězdy páru | prohoďte hvězdy nebo typy |
| O−C modelu má vůči naměřeným minimům špatnou fázi | model nemá kalendář | zadejte BJD času $t = 0$ (`--bjd-t0`) |
| pokaždé jiný výsledek ansámblu | změnilo se semínko | zafixujte `ensemble.seed` |

**Tabulka 78.** Výsledky, které vypadají chybně.

## Ansámbly {#ts-ensembles}

| Příznak | Příčina | Náprava |
| --- | --- | --- |
| přerušený ansámbl začíná od nuly | chybí `--resume` | spusťte týž příkaz s `--resume` |
| členy skončily stavem `error` | konfigurace, která pro některá losování selhává | přečtěte zprávu v `rows.jsonl`; po opravě použijte `--resume --retry-errors` |
| podmíněný výsledek stojí na několika málo členech | přesná pozorování, široké apriorní rozdělení | spusťte více členů, zužte apriorní rozdělení, nebo použijte emulátor či SBI |
| podíly se širokými intervaly | příliš málo členů | více členů zužuje 68% interval jako $1/\sqrt{n}$ |
| `require` odmítá většinu losování | apriorní rozdělení dává převážně nestabilní soustavy | zužte losované rozsahy |
| `observed columns not in the summary` | překlep v názvu sloupce, nebo sloupec, který běh nevytváří | porovnejte s tabulkou 71 |

**Tabulka 79.** Problémy s ansámbly.

## Proložení naměřených dat {#ts-fits}

| Příznak | Příčina | Náprava |
| --- | --- | --- |
| krok 2 přeskočí na špatný cyklus: rezidua celých period, $\chi^2$ nad počáteční hodnotou | počátek je od řešení příliš daleko, takže se minima očíslují jinak | nejprve ukotvěte (`anchor: true` nebo `--anchor`); v prvním proložení uvolněte jen `P` a `M0` a zbytek ve druhém, spuštěném z jeho `stage2_fitted.yaml`; zužte `bounds` |
| varování „strongly correlated parameters“ (silně korelované parametry) | data je neoddělí (např. $M_0$, $\omega$ a $P$ vnější dráhy pokryté méně než dvakrát) | velké korelované chyby jsou správné; ověřte je pomocí `sampler: emcee`, nebo jeden parametr zafixujte |
| krok 2 je pomalý | každé vyhodnocení integruje celý rozsah dat; každá iterace potřebuje na každý volný parametr dvě další | zvyšte `processes`; vynechte nepřesná minima (`max_sigma_d`); prokládejte méně parametrů |
| několik minim má rezidua mnoha σ | chybný typ nebo pár v tabulce, nebo odlehlé body | zkontrolujte tyto řádky v `stage2_residuals.csv`; `clip_sigma: 5`, robustní `loss` nebo `exclude` |
| redukované $\chi^2$ v řádu desítek až tisíců | uvedené chyby jsou hluboko pod skutečným rozptylem | `error_scale: per_source` s `error_model: jitter` |
| varování „the errors / outliers still changed after N fits“ (chyby či odlehlé body se po N proloženích stále měnily) | korekce se v rámci `refits` neustálily | zvyšte `refits`, nebo další proložení spusťte z `stage2_fitted.yaml` |
| rezidua primárních a sekundárních minim se liší o několik σ | chybný vektor excentricity páru (lokální minimum) | začněte znovu s $e\cos\omega$ opačného znaménka, nebo uvolněte `ecosw` a `esinw` |
| parametr na své mezi | řešení pravděpodobně není skutečné minimum, jeho chyba nemá smysl | rozšiřte `bounds` nebo přidejte apriorní rozdělení; ověřte, zda ho data omezují |
| proložená vnější dráha nesouhlasí s publikovanou | minima ji omezují jen slabě | přidejte apriorní rozdělení z publikovaného řešení; zpráva uvádí odchylku každého z nich v σ |
| rychlost apsidálního stáčení z kroku 1 je daleko od očekávané hodnoty, vysoké redukované $\chi^2$ | minima vzorkovaná v pravidelných intervalech cyklů připouštějí aliasy rychlosti | přidejte minima v jiných cyklech, nebo rychlost začněte či omezte blízko očekávané hodnoty ([Řešené příklady](#ex-fit1)) |
| `sampler: emcee needs the package emcee` | není nainstalována skupina `fit` | `pip install -e ".[fit]"` |

**Tabulka 80.** Problémy s proložením.

## Zpracovatel úloh a webová aplikace {#ts-worker}

| Příznak | Příčina | Náprava |
| --- | --- | --- |
| úlohy zůstávají ve frontě | zpracovatel úloh neběží, nebo jsou úlohy pozdrženy (`.yaml.hold`) | spusťte zpracovatele (stránka Jobs nebo `swingersim worker`); zařaďte úlohy do fronty |
| úloha běží se starým kódem | zpracovatel byl spuštěn před aktualizací | zastavte ho (stránka Jobs nebo prázdný soubor `jobs/stop`) a spusťte znovu |
| úloha selhala se zprávou o nedostatku paměti | zpracovatel zastavil úlohu, která používala nejvíc paměti | méně paralelních `processes`, méně výstupů, kratší běh |
| soubor ve frontě je ignorován | jeho název nekončí na `.yaml` | přejmenujte ho |
| grafy webové aplikace se nezobrazují | plotly není nainstalováno a není připojení k internetu | `pip install -e ".[viz]"` |
| stránka není dostupná z jiného počítače | server naslouchá na 127.0.0.1 | `--host 0.0.0.0` (jen v důvěryhodné síti) |

**Tabulka 81.** Problémy se zpracovatelem úloh a webovou aplikací.

## Časté otázky {#faq}

**Kolik hvězd může soustava mít?** Libovolná hierarchie vnořených dvojic; otestovány jsou soustavy
až do deseti hvězd. Výpočetní cena roste s počtem krátkých vnitřních drah, ne s počtem hvězd.

**Umí SWINGERsim proložit moji světelnou křivku?** V kroku 2 příkazu `swingersim fit` prokládá
segmenty světelných křivek spolu s okamžiky minim, s kulovými, rotujícími nebo Rocheovými hvězdami
a s kalibrací pro každý segment; jde o fotodynamické proložení celé soustavy, nikoli o náhradu
podrobného modelu povrchu jedné dvojhvězdy (skvrny, podrobné atmosféry).

**Lze počítat zpět v čase?** Ano, se záporným `run.t_max_yr` nebo s
`ensemble.directions: [backward]`. Disipativní procesy (slapy, vítr, přenos hmoty) nejsou vratné
v čase a chaotické soustavy mají zpětně smysl jen jako ansámbly (*Fyzika a metody*, kapitola
Omezení).

**Který backend použít?** N-částicový pro vše, co měří pozorovatel, a pro doby do tisíců let;
sekulární pro miliony let a velké ansámbly ([Dva backendy](#backends)).

**Který model hvězd použít?** `zams` pro rychlé dynamické studie nevyvinutých hvězd; `mist` pro
skutečné soustavy a pro vše, co stárne; `sse` pro rychlý vývoj s přenosem hmoty až do zapálení
helia. Naměřené poloměry a teploty v `stars` mají vždy přednost.

**Proč se liší okamžiky zákrytů „exact“ a „as an observer“?** První jsou konjunkce z dynamiky,
druhé minima změřená na syntetické světelné křivce. Nesouměrné zákryty, světlo ostatních hvězd,
pulzace, šum a vzorkování posouvají změřené minimum od konjunkce. Naměřená minima porovnávejte
s druhými.

**Jak přesný je SWINGERsim?** Každá fyzikální složka je srovnána s publikovanými výsledky
v dokumentu *Testování, část I*; aplikace na skutečné soustavy jsou v částech II–IV.

**Jak citovat výsledky?** Citujte články použitých metod (uvádí je *Methods paragraph* na stránce
výsledku) a uveďte konfiguraci a verzi SWINGERsimu, kterou zaznamenává každý výstupní soubor.
