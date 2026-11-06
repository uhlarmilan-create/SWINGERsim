# Výstupy {#outputs}

Tato kapitola popisuje každý soubor, který SWINGERsim zapisuje, význam jeho sloupců, jak hodnoty
interpretovat a jak soubory číst v Pythonu. Všechny tabulky jsou soubory CSV, jejichž první řádky
mohou začínat znakem `#`: zaznamenávají verze softwaru, jednotky, časový systém a u některých
i konfiguraci; první řádek, který nezačíná `#`, je hlavička.

## Souhrn běhu {#out-summary}

`swingersim run` vypíše souhrn, který se také uloží jako `<name>_summary.txt`. Příklad 3 z kapitoly
[Řešené příklady](#ex-tides) (soustava 2+2 se slapy a rotací) dává ve zkrácené podobě:

```text
integrated to t = 300 yr (1500 outputs, osculating elements)
  out    a = 2.442..2.4431 AU  e = 0.3998..0.4002
  A      a = 0.041764..0.041767 AU  e = 0.0494..0.0502
  i_mut A-B           8.44..  25.08 deg
  max |dE/E| = 2.22e-11
  max |dL/L| = 1.83e-07
hierarchy: configured [[1,2],[3,4]]
  no change detected
orbit pairs:
  A in out: i_mut0 = 11.6 deg, g1_0 = 98 deg
    quadrupole theory : e = 0.048..0.050, i_mut = 11.6..11.6 deg, e-cycle 323 yr (outside ZLK window)
    simulation        : e = 0.049..0.050, i_mut = 11.6..13.6 deg
    precession period : 1.29e+03 yr estimated, not completed in the run, closest periastron = 2.78 x (R1+R2)
    stability ratio Mardling-Aarseth = 7.97 (stable), Vynatheya+22 = 10.17 (stable), eps_oct = 0.00054, eps_SA = 0.0017
    GR apsidal period ~ 2.31e+03 yr
    tidal+rotational apsidal period 64 yr, ZLK time 949 yr -> faster than ZLK, ZLK suppressed
binaries A and B: P_A/P_B = 0.6207 (nearest 3:5, +3.45 %), precession ratio A/B = 1.97 estimated
  mutual inclination of A and B: 8.44..25.08 deg, analytic limit 104.21 deg (Vokrouhlicky 2016, circular A and B)
  eclipses (observer along z): A 100 %, B 91 %, both at once 91 % of the run
secular coupling beta (Hamers & Lai 2017): A 1.45, B 0.365
stars:
  m1  R = 1.7 Rsun, k2 = 0.0144, rg2 = 0.0754, lag = 0 s, P_rot = 1.8 -> 1.8 d, obliquity 0.0 -> 0.1 deg, Roche lobe filled up to 0.51
  binary A (m1+m2): apsidal period 62.29 yr from the binary itself, tides+rotation 64.01 yr, GR 2311 yr; measured 59.43 yr (with the other stars)
  wall time 152.9 s
```

| Řádek | Obsah |
| --- | --- |
| `integrated to` | dosažený čas, počet výstupů a druh elementů (oskulační pro N-částicový backend, střední pro sekulární); kratší než `t_max_yr`, pokud běh zastavila událost |
| `EVENT` | jen tehdy, když se běh zastavil předčasně: `collision` (s hvězdami), `roche_overflow`, `hierarchy_change`, `energy_error`, `wall_time`, `e_reached`, `stellar_evolution_end`, `stellar_birth`, `gw_unresolved` nebo `mass_transfer: <kind>` (viz níže) |
| `MERGER`, `MASS TRANSFER STARTED BY THE DYNAMICS` | s `physics.mergers` nebo `mass_transfer.dynamic_onset`: co se stalo a jak běh pokračoval |
| jeden řádek na dráhu | minimum..maximum velké poloosy [au] a excentricity |
| `i_mut X-Y` | minimum..maximum vzájemného sklonu dvou drah [deg] |
| max \|dE/E\|, max \|dL/L\| | největší relativní chyba energie a celkového momentu hybnosti; míra numerické přesnosti (s disipací se energie mění fyzikálně a poznámka to uvede) |
| `hierarchy` | nastavená hierarchie a každá potvrzená změna: čas, druh, stará → nová struktura |
| `X in Y` | každá vnitřní dráha ve své nadřazené dráze: počáteční vzájemný sklon a argument periastra měřený od uzlové přímky |
| `quadrupole theory` | předpověď kvadrupólových rovnic středovaných přes oběh: rozsahy $e$ a $i_\mathrm{mut}$ a perioda cyklu excentricity; zda dráha leží uvnitř okna ZLK (39,2–140,8°) |
| `simulation` | tytéž rozsahy naměřené v běhu (maximum $e$ zahrnuje vrcholy mezi výstupy) a počet dokončených cyklů ZLK |
| `precession period` | odhadnutá a naměřená perioda precese uzlů; nejbližší periastrum v jednotkách součtu poloměrů hvězd |
| `stability ratio`, `eps_oct`, `eps_SA` | poměry stability podle Mardlinga & Aarsetha (2001) a Vynatheyi et al. (2022) pro tuto vnořenou trojici (nad 1 stabilní); síla oktupólu; parametr jednoduchého středování – nad zhruba 0,01 ztrácejí středované odhady přesnost (Luo et al. 2016) a poznámka to uvede |
| `GR apsidal period`, `tidal+rotational apsidal period` | relativistická a (s `physics.tides`) slapová a rotační perioda apsidálního stáčení ve srovnání s časovou škálou ZLK; poznámka uvede, kdy jsou cykly ZLK potlačeny |
| `binaries A and B` | jádra 2+2: poměr period a nejbližší soudělnost malých celých čísel, poměr period precese; rozsah vzájemného sklonu A a B a jeho analytická mez (Vokrouhlický 2016) |
| `eclipses` | podíl běhu, během něhož každá dvojhvězda, a obě současně, zakrývá pro pozorovatele podél $z$ |
| `secular coupling beta` | parametr $\beta$ podle Hamerse & Laie (2017) pro každou vnitřní dráhu soustavy 2+2 |
| `stars` | poloměr, $k_2$, $I/(MR^2)$, zpoždění, rotační perioda a sklon osy rotace na začátku a na konci, největší zaplnění Rocheova laloku; u vyvíjejících se hvězd stáří, poloměr, $T_\mathrm{eff}$ a fáze na začátku a na konci |
| `binary X` | každá dvojhvězda ze dvou hvězd: perioda apsidálního stáčení očekávaná ze samotné dvojhvězdy a naměřená v běhu (která zahrnuje poruchy od ostatních hvězd) |
| `tidal time scales` | se slapovým zpožděním: časy zakulacení, rozpadu dráhy a synchronizace, pseudosynchronní rotační perioda |
| `MEGNO` | s `run.megno`: ukazatel chaosu a Ljapunovův čas |
| `wall time` | doba výpočtu |

**Tabulka 65.** Řádky souhrnu běhu. Vlevo: začátek řádku; vpravo: jeho obsah.

Druhy událostí `mass_transfer` jsou u dvojhvězdy `common_envelope`, `contact`,
`contact_merger_darwin`, `contact_merger_l2`, `common_envelope_merger`, `post_common_envelope`
a `donor_stripped`; u přenosu z třetí hvězdy `triple_common_envelope`, `inner_roche_overflow`,
`dynamical_instability` a `donor_stripped`.

## Soubory běhu {#out-run}

`swingersim run config.yaml -o out/sys.h5` zapíše vedle sebe:

| Soubor | Obsah |
| --- | --- |
| `sys.h5` | všechna data běhu (viz níže) |
| `sys.html` | zpráva HTML: příkazový řádek, konfigurace, souhrn, interaktivní grafy, obrázky |
| `sys_summary.txt` | souhrn a seznam souborů běhu |
| `sys_timeseries.csv` | jeden řádek na okamžik výstupu (nejvýše 20 000 řádků; u delších běhů každý k-tý výstup, v souboru HDF5 všechny) |
| `sys_viewer.html` | interaktivní 3D prohlížeč (vyžaduje plotly) |
| `sys_orbits.png` | $a/a_0$, $e$ a sklon k obloze každé dráhy v čase |
| `sys_inclinations.png` | vzájemné sklony v čase |
| `sys_conservation.png` | relativní chyby energie a momentu hybnosti |
| `sys_stars.png` | hmotnost, poloměr a $T_\mathrm{eff}$ hvězd v čase (jen u vyvíjejících se hvězd) |

**Tabulka 66.** Standardní produkty běhu. Každý běh spuštěný zpracovatelem úloh nebo webovou
aplikací zapisuje tytéž soubory; `--no-viewer`, `--no-plots` a `--no-html` vynechají prohlížeč,
obrázky a zprávu.

Sloupce souboru `sys_timeseries.csv` (názvy drah jako v konfiguraci; úhly ve stupních, na rozdíl od
radiánů v souboru HDF5):

| Sloupec | Obsah |
| --- | --- |
| `t_yr` | čas [yr] od $t = 0$ (záporný u běhů dozadu) |
| `a_<orbit>_au`, `e_<orbit>` | velká poloosa a excentricita každé dráhy |
| `i_<orbit>_deg`, `i_inv_<orbit>_deg` | sklon k obloze a k invariantní rovině (rovině kolmé k celkovému momentu hybnosti) |
| `Omega_<orbit>_deg`, `Omega_inv_<orbit>_deg`, `omega_<orbit>_deg` | uzel na obloze a na invariantní rovině, argument periastra |
| `imut_<o1>-<o2>_deg` | vzájemný sklon každé dvojice drah |
| `dE_over_E`, `dL_over_L` | relativní chyby energie a momentu hybnosti |
| `M_<star>_msun`, `R_<star>_rsun`, `Teff_<star>_k` | hmotnost, poloměr a $T_\mathrm{eff}$ (jen u vyvíjejících se hvězd) |

**Tabulka 67.** Sloupce souboru CSV s časovými řadami.

### Soubor HDF5 {#out-hdf5}

HDF5 je standardní binární formát pro numerická data, čitelný z Pythonu (`h5py`) a z nástrojů, jako
jsou HDFView, MATLAB nebo Julia.

| Cesta | Obsah | Jednotka |
| --- | --- | --- |
| `time_yr` | okamžiky výstupů | yr |
| `energy_error`, `angmom_error` | relativní chyby energie a celkového momentu hybnosti (dráhy + rotace) | – |
| `orbits/<orbit>/a`, `e` | velká poloosa, excentricita | au, – |
| `orbits/<orbit>/inc`, `Omega`, `omega`, `M` | sklon k obloze, uzel, argument periastra, střední anomálie | **rad** |
| `orbits/<orbit>/inc_inv`, `Omega_inv` | sklon a uzel vůči invariantní rovině | **rad** |
| `mutual_inclination/<o1>-<o2>` | vzájemný sklon | **rad** |
| `stars/<m>` (atributy) | použité hodnoty: `radius_rsun`, `apsidal_k2`, `rg2`, `lag_s`, `mass_msun`, `mass_msun_config`, zaplnění Rocheova laloku (`roche_fill_max`, jeho čas) a z modelu hvězdy `luminosity_lsun`, `teff_k`, `age_myr`, `phase`, hmotnost a poloměr obálky, doba obratu konvekce | R☉, –, –, s, M☉ … |
| `stars/<m>/P_rot_d`, `obliquity_deg` | rotační perioda a úhel mezi osou rotace a dráhou v každém výstupu (s `physics.tides`) | d, deg |
| `stars/<m>/radius_rsun`, `luminosity_lsun`, `teff_k`, `phase`, `age_myr`, `mass_msun` | vyvíjející se hvězda v každém výstupu (s `physics.evolution` / `mass_loss`); fáze 0 před hlavní posloupností, 1 hlavní posloupnost, 2 Hertzsprungova mezera, 3 větev obrů, 4 hoření helia v jádře, 5 raná AGB, 6 tepelně pulzující AGB | R☉, L☉, K, –, Myr, M☉ |
| `megno` | MEGNO v každém výstupu (s `run.megno`) | – |
| `extremes/<orbit>` (atributy) | `e_max`, `t_e_max_yr`, `q_min_au`, `t_q_min_yr`: extrémy nalezené mezi výstupy | –, yr, au |
| `hierarchy/detected` | rozpoznaná hierarchie v každém výstupu, např. `[[1,2],[3,4]]` | text |
| `hierarchy` (atributy) | `configured`, `final`, `transitions_json` (potvrzené změny: čas, druh, z, do), `detected_changes_json` | – |
| `diagnostics/<inner orbit>` (atributy) | diagnostiky každé dvojice drah (tabulka 69) | – |

**Tabulka 68.** Obsah souboru HDF5 běhu. Úhly jsou uloženy v radiánech.

Atributy souboru: `config_yaml` (úplná konfigurace se všemi výchozími hodnotami, aby bylo možné běh
zopakovat), `swingersim_version` a verze Pythonu, NumPy, SciPy, h5py, REBOUND, REBOUNDx a emcee,
`created_utc`, `wall_s`, `units_json` (zvyklosti: juliánské roky od $t = 0$, au, M☉, radiány,
soustava oblohy), `backend`, `event`, `event_time_yr` a hvězdy a stáří při události, `dissipative`,
`gas_phase`, `mass_loss` (true, když se energie fyzikálně nezachovává), `binaries_json` (pro každou
dvojhvězdu ze dvou hvězd: očekávané a naměřené periody apsidálního stáčení, slapové časové škály, podíl
zákrytů), `system_json` (poměry period a soudělnosti jader 2+2, analytická mez vzájemného sklonu,
podíly zákrytů, vazební parametry, indexy ztráty hmoty, historie přenosu hmoty a společné obálky
každé dráhy, počty setkání s míjejícími hvězdami) a `lyapunov_time_yr` (s `run.megno`).

**Zápis hierarchie.** Závorky seskupují vázané dvojice: `[[1,2],[3,4]]` je 2+2, `[[[1,2],3],4]` je
3+1; nevázané části jsou odděleny znakem `|`, např. `[[1,2],3] | 4` je trojhvězda s unikající
hvězdou. Druhy změn jsou `exchange` (výměna: stejný typ, jiní partneři), `reconfiguration`
(přestavba: typ se mění, např. 2+2 → 3+1), `ejection` (vymrštění: jedna hvězda se stane nevázanou),
`breakup` (rozpad: několik nevázaných částí) a `capture` (zachycení: nevázané části se stanou
vázanými).

| Pole | Význam |
| --- | --- |
| `i_mut0_deg`, `g1_0_deg` | počáteční vzájemný sklon vůči nadřazené dráze; počáteční argument periastra od uzlu na místní invariantní rovině |
| `in_zlk_window`, `e_max_quadrupole` | zda $i_\mathrm{mut,0}$ leží mezi 39,2° a 140,8°; odhad pro testovací částici $\sqrt{1-\frac53\cos^2 i_0}$ |
| `t_zlk_yr`, `t_prec_yr`, `t_gr_yr`, `gr_suppresses_zlk` | časová škála ZLK, perioda precese uzlů, relativistická perioda apsidálního stáčení a zda je relativita rychlejší než ZLK |
| `epsilon_oct`, `epsilon_sa` | síla oktupólu (Naoz 2016) a parametr jednoduchého středování (Luo et al. 2016) |
| `stability_ratio`, `stability_ratio_v22` | poměry podle Mardlinga–Aarsetha a Vynatheyi et al. (2022) |
| `e_min`, `e_max`, `i_mut_min_deg`, `i_mut_max_deg` | rozsahy v běhu |
| `e_max_true`, `t_e_max_yr`, `q_min_au`, `q_min_over_radii` | největší excentricita včetně vnitřních vzorků, její čas; nejmenší periastrum, také v jednotkách $R_1+R_2$ |
| `quad_e_min`, `quad_e_max`, `quad_i_min_deg`, `quad_i_max_deg`, `quad_period_yr` | předpověď středovaných kvadrupólových rovnic (Naoz et al. 2013) |
| `n_zlk_cycles`, `zlk_period_yr`, `prec_period_measured_yr` | dokončené cykly excentricity a jejich perioda; perioda precese naměřená vůči invariantní rovině |

**Tabulka 69.** Diagnostiky dvojice drah (atributy `diagnostics/<inner orbit>`).

## Ansámbly {#out-ensemble}

`swingersim ensemble config.yaml -o folder` zapíše do složky:

| Soubor | Obsah |
| --- | --- |
| `ensemble.yaml` | úplná konfigurace a její otisk (navázaný ansámbl musí používat tutéž) |
| `rows.jsonl` | jeden řádek na dokončeného člena a směr, zapisovaný průběžně, jak členy končí; `--resume` z něj pokračuje |
| `summary.csv` | tytéž řádky jako tabulka, seřazené podle člena a směru |
| `summary.parquet` | tatáž tabulka pro pandas (s nainstalovanou skupinou `ensemble`) |
| `runs/member_<n>_<direction>.h5` | úplný běh každého člena (jen se `save_runs: true`) |
| `report.html` | zpráva |

**Tabulka 70.** Soubory ansámblu.

| Sloupec | Význam |
| --- | --- |
| `member`, `direction` | číslo člena; `forward` nebo `backward` |
| `masses.m3`, `orbits.out.P`, … | vylosované hodnoty parametrů z `vary` |
| `fate` | `intact` (žádná událost, hierarchie beze změny), změna hierarchie (`hierarchy_change`, `exchange`, `ejection`, `breakup`, `reconfiguration`, `capture`), událost (`collision`, `merger`, `roche_overflow`, `energy_error`, `wall_time`, `e_reached`, `gw_unresolved`) nebo `error` |
| `event`, `event_time_yr`, `t_end_yr` | událost, její čas, konec běhu |
| `n_transitions`, `first_transition_yr`, `final_hierarchy` | změny hierarchie |
| `<orbit>_e0`, `_e_min`, `_e_max`, `_log1me_min` | počáteční, nejmenší, největší excentricita; $\log_{10}(1-e_\mathrm{max})$ |
| `<orbit>_t_e_max_yr`, `_q_min_au`, `_a_min`, `_a_max` | čas největší excentricity, nejbližší periastrum, rozsah $a$ |
| `<orbit>_t_1me_lt_1e-1_yr` … `_1e-4_yr` | první okamžik, kdy $1-e$ kleslo pod $10^{-1}$ … $10^{-4}$ (prázdné, pokud nikdy) |
| `<orbit>_U_measured_yr` | perioda apsidálního stáčení naměřená za běh (prázdná, pokud se periastrum posunulo příliš málo) |
| `imut_<o1>_<o2>_0`, `_min`, `_max` | počáteční, nejmenší, největší vzájemný sklon [deg] |
| `<binary>_eclipse_fraction`, `eclipse_fraction_both`, `eclipse_fraction_all` | podíl běhu, během něhož dvojhvězda (obě dvojhvězdy prvního jádra 2+2; všechny dvojhvězdy) zakrývá |
| `<binary>_i0`, `_didt0_deg_yr`, `_dtiltdt0_deg_yr` | současný sklon k obloze; rychlost jeho změny; rychlost změny $\lvert i - 90°\rvert$ (záporná: natáčení k pohledu z boku, hlubší zákryty); prokládá se jednostranně od $t = 0$ po dobu `run.rate_window_yr` |
| `<binary>_eclipsing0`, `_eclipse_change_yr` | 1, pokud dvojhvězda nyní zakrývá; první okamžik, kdy se stav zákrytů změní |
| `<star>_roche_fill_max`, `roche_first`, `roche_first_kind` | největší zaplnění Rocheova laloku; s `roche_overflow` hvězda, která svůj lalok vyplnila první, a zda je jejím průvodcem hvězda, nebo podsoustava |
| `mass_transfer_*` | skalární položky přenosu hmoty a společné obálky v běhu |
| `ma01_min`, `v22_min` | nejmenší poměr stability všech vnořených trojic v $t = 0$ |
| `mlp_stable`, `ghost_stable`, `ghost_delta_max` | s `mlp` / `ghost_test`: 1 stabilní, 0 nestabilní |
| `max_energy_error`, `wall_s` | kvalita a výpočetní cena běhu |

**Tabulka 71.** Sloupce souhrnu ansámblu. Tyto názvy se používají v `ensemble.observed`.

`swingersim ensemble-report` vypíše, zvlášť dopředu a dozadu, podíl každého osudu, rozdělení
největších excentricit a podíly s $1-e$ pod 0,1, 0,01 a 0,001, vzájemné sklony, současné rychlosti
změn a to, kdy zákryty končí nebo začaly, podíly zákrytů a verdikty stability; podíly jsou uvedeny
s 68% Wilsonovým intervalem (Wilson 1927). S pozorováními je každé číslo vážené, záhlaví uvádí
efektivní počet členů a vylosované parametry jsou uvedeny jako medián a rozsah 16–84 % apriorní →
podmíněný.

## Okamžiky zákrytů, minima, světelné křivky a radiální rychlosti {#out-observables}

**`swingersim timings`** zapíše jeden řádek na zákryt:

| Sloupec | Význam |
| --- | --- |
| `binary`, `kind` | název dráhy; `primary` (zakrývána je první hvězda dvojhvězdy) nebo `secondary` |
| `epoch` | číslo cyklu lineární efemeridy (u sekundárních zákrytů poločíselné) |
| `t_yr`, `t_obs_yr` | okamžik středu zákrytu [yr od $t = 0$], bez doby letu světla a s ní |
| `o_c_d` | pozorovaný minus vypočtený [d] vůči lineární efemeridě primárních zákrytů |
| `impact` | promítnutá vzdálenost ve středu zákrytu / $(R_1+R_2)$: 0 centrální, blízko 1 tečný |
| `T0_yr`, `P_d` | proložená lineární efemerida |
| `bjd_obs` | s `--bjd-t0`: BJD minima |

**Tabulka 72.** Sloupce souboru CSV s okamžiky zákrytů. Souhrn na obrazovce uvádí počet zákrytů,
periodu, rozsah O−C a intervaly bez zákrytů; s `physics.activity` uvádí také energetickou kontrolu
každé aktivní hvězdy.

**`swingersim minima`** zapíše jeden řádek na změřené minimum:

| Sloupec | Význam |
| --- | --- |
| `binary`, `kind`, `epoch` | jako u okamžiků zákrytů |
| `time`, `t_min_d`, `sigma_d` | okamžik minima (BJD_TDB a dny po $t = 0$) z profilu podle Mikuláška (2015) a jeho směrodatná chyba |
| `t_geom_d` | geometrická konjunkce s dobou letu světla [d] |
| `o_c_d`, `o_c_geom_d` | O−C minim a konjunkcí vůči jejich vlastním lineárním efemeridám [d] |
| `depth_mag`, `width_d`, `rms_mag` | hloubka a šířka tohoto zákrytu, rezidua jeho proložení |

**Tabulka 73.** Sloupce souboru CSV s minimy. `--seasons` zapíše tvar pro každou sezónu (`t_mid_d`,
`n`, `depth_mag`, `width_d`, `gamma`, `corr`, `rms_mag`).

**`swingersim lightcurve`** zapíše jeden řádek na okamžik: `t_day` (dny po $t = 0$, čas
pozorovatele), `flux` (celkový tok normovaný na součet hvězd mimo zákryt), `flux_m1`, … (příspěvek
každé hvězdy; u pulzujících hvězd včetně pulzace, u hvězd za diskem včetně jeho absorpce)
a `flux_disc_<name>` (světlo disků s emisí nebo rozptylem). S `--engine phoebe` jsou sloupci
`flux_A`, `flux_B`, … (každá vnitřní dvojhvězda), osamocené hvězdy a `flux_mutual` (světlo ztracené
při zákrytech mezi různými dvojhvězdami). Hlavičkové řádky uvádějí pásmo, časovou konvenci a BJD
času $t = 0$, je-li znám.

**`swingersim rv`** zapíše `rv_model.csv` (sloupce `bjd_tdb`, `rv_m1`, `rv_m2`, …, km/s) a obrázek
`rv_model.png` s naměřenými rychlostmi, jsou-li zadány.

**`swingersim apsidal`** zapíše pro každý pár periodu, excentricitu, $k_2$, poloměry a rotační
periody obou hvězd, slapové stáčení, rotační stáčení pro synchronní, pseudosynchronní a nastavenou
rotaci, relativistické stáčení, celkové stáčení (deg/yr) a periodu apsidálního stáčení $U$ (yr).

## Proložení {#out-fit}

`swingersim fit` zapisuje do `fit.output.dir`. Krok 1:

| Soubor | Obsah |
| --- | --- |
| `stage1_parameters.csv`, `.tex` | proložené parametry s chybami 1σ (`<pair>.T0`, `.P`, `.Q`, `.e`, `.w0`, `.dw` [rad/cyklus], `.asini` [au], `.P3`, `.T3`, `.e3`, `.w3`; `out.*` pro společnou vnější dráhu) a odvozené hodnoty: `.K_s` (poloviční amplituda efektu světelného času), `.f_m3` (hmotnostní funkce), `.P3_yr`, `.U_yr` (perioda apsidálního stáčení), `.Pa_d` (anomalistická perioda), `mass_ratio_A_over_B` |
| `stage1_oc.csv` | každé minimum: pár, typ, epocha, okamžik, sigma, zdroj, O−C vůči lineární efemeridě, model, reziduum [d] |
| `stage1_summary.json` | $\chi^2$, počet stupňů volnosti, parametry, odvozené hodnoty, škálování chyb, odlehlé body, efemeridy, verze softwaru |
| `stage1_predictions.csv` | s `output.predict`: budoucí minima a jejich rozptyl |
| `stage1_oc.png`, `stage1_corner.png` | O−C s modelem a rezidui a jednotlivými členy zvlášť; korelační graf vzorků MCMC |

**Tabulka 74.** Soubory proložení v kroku 1.

Krok 2:

| Soubor | Obsah |
| --- | --- |
| `stage2_params.csv` | každý volný parametr: `parameter`, `start`, `value`, `error` (1σ; s MCMC polovina rozsahu 16–84 %), `unit`, `prior_mean`, `prior_sigma`, `error_linear` (z kovarianční matice), `p16`, `p50`, `p84` (percentily MCMC) |
| `stage2_covariance.csv` | kovarianční matice lineárního modelu |
| `stage2_residuals.csv` | každé minimum: `pair`, `type`, `bjd_tdb`, `sigma_d` (použitá chyba), `model_start_bjd`, `model_bjd`, `o_minus_model_min`, `source`, `sigma_quoted_d` (chyba z tabulky), `status` (`fitted` nebo `outlier`) |
| `stage2_rv.csv` | každá rychlost: `star`, `bjd_tdb`, `rv`, `sigma`, `sigma_quoted`, `model_start`, `model`, `o_minus_model`, `source`, `status` |
| `stage2_lc.csv` | každý bod světelné křivky: `source`, `segment`, `bjd_tdb`, `flux`, `sigma`, `sigma_quoted`, `model_start`, `model`, `o_minus_model` |
| `stage2_inclination.csv` | každá epocha sklonu: `label`, `pair`, `time`, `i_obs`, `sigma`, `limit`, `i_start`, `i_model`, `i_eclipse_limit`, `z` (normované reziduum), `source` |
| `stage2_fitted.yaml` | konfigurace s proloženými hodnotami a `fit_epoch_bjd_tdb`, bez oddílu `fit`; její hlavička uvádí $\chi^2$, rms pro každý pár a epochu; běží tak, jak je |
| `stage2_report.txt` | $\chi^2$ od počátku k proložení, rms pro každý pár a druh dat, parametry s chybami, silně korelované dvojice ($\lvert r\rvert > 0{,}95$), faktory chyb nebo jittery jednotlivých zdrojů, apriorní rozdělení a odchylky od nich, proložení v řadě a odlehlé body; tentýž text se vypíše |
| `stage2_history.json` | $\chi^2$ a parametry při každém vyhodnocení modelu |
| `stage2_samples.csv` | s MCMC: vzorky po zahoření (burn-in) |
| `stage2_predictions.csv`, `.png` | předpovězená minima: `pair`, `type`, `bjd_tdb`, `sigma_min`, `bjd_16`, `bjd_84`, `duration_h`, `impact`, pozorovací okno a se souřadnicemi i časy v UTC |
| `stage2_oc.png`, `stage2_rv.png`, `stage2_lc.png`, `stage2_inclination.png` | obrázky O−C, rychlostí, světelných křivek a sklonů s počátečním a proloženým modelem a rezidui |
| `anchored.yaml` | ukotvená počáteční konfigurace (s `anchor` nebo `--anchor`) |
| `report.html` | zpráva |

**Tabulka 75.** Soubory proložení v kroku 2.

## Rozptylové experimenty a mapy chaosu {#out-scatter}

`swingersim scatter` zapíše `scatter.yaml` (konfiguraci a kritickou rychlost), `rows.jsonl` (jeden
řádek na setkání: srážkový parametr, výsledek, stav, konečná hierarchie, koncový čas, chyba energie,
vázané dráhy na konci, stabilita vázané podsoustavy, srážející se hvězdy, celková energie a moment
hybnosti, vylosované počáteční podmínky) a `summary.txt` (počet, podíl s 68% Wilsonovým intervalem
a účinný průřez každého výsledku).

`swingersim chaosmap` zapíše `chaosmap.json` (mřížku a konfiguraci), `grid.jsonl` (jeden řádek na
bod mřížky: indexy, hodnoty, `megno`, `lyapunov_yr`, `t_end_yr`, `status` – `ok`, `escape`,
`collision` nebo `wall_time`), `summary.txt` (počty pravidelných, chaotických a unikajících bodů
a mez stability podél x) a `map.png` (MEGNO na mřížce, oříznuté na 10; úniky a srážky vykresleny
jako 10). `lyapunov_yr` (a `lyapunov_time_yr` běhu) je Ljapunovův čas $1/\lambda$: funkce REBOUND
`sim.lyapunov()` vrací sklon středního MEGNO, $\langle Y\rangle \approx \lambda t/2$, takže SWINGERsim uvádí
$1/(2\,\mathrm{sklon})$. Počáteční tečný vektor MEGNO má pevné semínko, takže mapa vyjde při každém běhu stejně.

## Strojové učení {#out-ml}

`emulator-train` zapíše model (`.pkl`) a vedle něj `model.json` (sloupce, skóre, nastavení) a s
`--plot` i jeho diagnostiku; `emulator-report` vypíše aposteriorní rozdělení každého vstupu
s efektivní velikostí vzorku. `sbi` vypíše aposteriorní rozdělení (medián, 68% a 95% interval)
a volitelně zapíše síť (`--model`), vzorky (`--save-samples`), korelační graf a kalibrační
diagnostiku.

## Jak interpretovat hodnoty {#out-interpret}

- **Oskulační elementy** N-částicového backendu kmitají na časových škálách oběhu (typicky pod 1 %
  v $a$ a 0,01 v $e$); sledujte trendy a extrémy za mnoho vnějších oběhů. Sekulární backend dává
  střední elementy bez těchto zákmitů.
- **Velké poloosy** zůstávají v dynamice hmotných bodů v průměru stálé (se ztrátou hmoty zůstává
  stálé $aM$). Drift spolu s rostoucí chybou energie ukazuje na numerický problém; náhlý skok na
  těsné přiblížení.
- **Sklon k obloze** rozhoduje o tom, zda dvojhvězda zakrývá; **vzájemný sklon** rozhoduje
  o dynamice (precese, ZLK).
- **Ω a ω** mohou přeskočit přes 360°; před vykreslením je rozviňte.
- **NaN** znamená, že dráha v daném výstupu nebyla vázaná, např. po rozpadu.
- **Periody apsidálního stáčení** naměřené v běhu zahrnují poruchy od ostatních hvězd a mohou se silně
  lišit od hodnoty izolované dvojhvězdy.
- **Zaplnění Rocheova laloku** blízko 1 znamená, že hvězda v periastru téměř vyplňuje svůj lalok.
- **Chyba energie.** Bez relativity je normální $10^{-14}$–$10^{-12}$; s relativitou se kontroluje
  energie 1PN a typické je $10^{-13}$–$10^{-11}$; hodnoty blížící se `max_energy_error` ($10^{-5}$)
  znamenají, že výsledku nelze věřit. S disipací, plynem nebo ztrátou hmoty se energie mění
  fyzikálně.
- **Primární zákryt** je ve SWINGERsimu zákryt první hvězdy dvojhvězdy bez ohledu na hloubku;
  pozorovatel obvykle nazývá primárním hlubší zákryt.
- **Modelový čas** začíná v $t = 0$, v epoše elementů; naměřené BJD se porovnávají po zadání BJD
  času $t = 0$ (`--bjd-t0`, `fit.stage2.bjd_t0`).

## Čtení souborů v Pythonu {#out-python}

```python
import h5py
import numpy as np
import pandas as pd

with h5py.File("output/example.h5") as f:
    t = f["time_yr"][:]
    e_A = f["orbits/A/e"][:]
    i_AB = np.degrees(f["mutual_inclination/A-B"][:])      # radians in the file
    print(f.attrs["event"] or "no event")
    print(f.attrs["config_yaml"])

ts = pd.read_csv("output/example_timeseries.csv", comment="#")
ens = pd.read_csv("output/v994/summary.csv")
fwd = ens[ens.direction == "forward"]
print(fwd.fate.value_counts(normalize=True))
res = pd.read_csv("output/example_fit/stage1_oc.csv", comment="#")
```
