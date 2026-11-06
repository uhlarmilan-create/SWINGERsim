# Konfigurační soubor {#configuration}

Konfigurace je textový soubor YAML, který popisuje jednu soustavu a to, co se s ní má dělat. Tato
kapitola je jejím úplným přehledem: každý oddíl, každý klíč s jeho typem, výchozí hodnotou,
přípustnými hodnotami a významem, tak jak je program kontroluje. Klíč, který zde není uveden, je
odmítnut, stejně jako hodnota mimo svůj rozsah; chybové hlášení příslušný klíč jmenuje.

## Jednotky a zvyklosti {#cfg-units}

| Veličina | Jednotka v souboru |
| --- | --- |
| hmotnosti | M☉ |
| periody drah a hvězd | dny |
| velké poloosy, poloměry disků | au |
| poloměry hvězd | R☉ |
| úhly (sklon, uzel, argument pericentra, střední anomálie, sklon osy rotace) | stupně |
| časy běhů, časové škály | roky; stáří v Myr |
| okamžiky měření | BJD_TDB (pokud `fit.time_system` neurčuje jinak) |
| rychlosti | km/s |

**Tabulka 44.** Jednotky konfiguračního souboru. Vnitřně SWINGERsim pracuje v M☉, au a letech.

**Vztažná soustava.** Sklony a uzly se měří od roviny oblohy; osa $z$ míří podél zorného paprsku
směrem od pozorovatele, takže $i = 90°$ znamená pohled z boku (zákryty možné) a $i = 0°$ pohled
shora. Hodnoty $i$ nad 90° popisují tutéž dráhu obíhanou v opačném smyslu; na tom záleží u vzájemných
sklonů, ne u zákrytů. Každá dráha je relativní dráhou své druhé složky kolem první, přičemž obě
složky jsou zastoupeny svými těžišti (Jacobiho souřadnice). Zadané elementy jsou oskulační elementy
v čase $t = 0$, v epoše konfigurace.

**Názvy.** Hvězdy jsou `m1`, `m2`, … bez mezer v číslování. Dráhy mají názvy podle hierarchie: `A`,
`B`, `out` pro 2+2; `inner`, `mid`, `out` pro 3+1; vaše vlastní názvy pro `custom`.

## Nejvyšší úroveň {#cfg-top}

| Klíč | Typ, výchozí hodnota | Význam |
| --- | --- | --- |
| `topology` | text, povinný | `2+2`, `3+1`, `2+2+1`, `3+1+1`, `2+(2+1)`, `2+2+2`, `3+3` nebo `custom` |
| `hierarchy` | zobrazení dráha → [potomek 1, potomek 2]; jen s `custom` | strom: každá dráha má dva potomky, hvězdu nebo jinou dráhu; právě jedna dráha není ničím potomkem |
| `masses` | zobrazení hvězda → M☉, povinný | jedna kladná hmotnost pro každou hvězdu |
| `orbits` | zobrazení dráha → dráha, povinný | jedna položka pro každou dráhu hierarchie (viz níže) |
| `stars` | zobrazení hvězda → hvězda, výchozí prázdné | vlastnosti nad rámec hmotnosti (viz níže) |
| `discs` | zobrazení název → disk, výchozí prázdné | prachové disky a prstence, které ve světelných křivkách zakrývají hvězdy (viz níže) |
| `stellar_model` | oddíl | odkud pocházejí poloměry, teploty a stavba hvězd (viz níže) |
| `physics` | oddíl | fyzika nad rámec newtonovských hmotných bodů (viz níže) |
| `run` | oddíl, povinný (kromě konfigurace obsahující jen `fit`) | jak se počítá pohyb (viz níže) |
| `ensemble` | oddíl, volitelný | ansámbl Monte Carlo (`swingersim ensemble`); `run` ho ignoruje |
| `fit` | oddíl, volitelný | proložení naměřených dat (`swingersim fit`); `run` ho ignoruje |
| `fit_epoch_bjd_tdb` | číslo, volitelný | zapíše ho proložení v kroku 2 do `stage2_fitted.yaml`: BJD_TDB času $t = 0$ proložených elementů; běhy ho ignorují |

**Tabulka 45.** Klíče na nejvyšší úrovni konfigurace.

Předvolby klíče `topology` jsou zkratkami pro tyto stromy (obrázek 1):

| Předvolba | Strom |
| --- | --- |
| `2+2` | `out: [A, B]`, `A: [m1, m2]`, `B: [m3, m4]` |
| `3+1` | `out: [mid, m4]`, `mid: [inner, m3]`, `inner: [m1, m2]` |
| `2+2+1` | `out5: [out, m5]`, `out: [A, B]`, `A: [m1, m2]`, `B: [m3, m4]` |
| `3+1+1` | `out5: [out, m5]`, `out: [mid, m4]`, `mid: [inner, m3]`, `inner: [m1, m2]` |
| `2+(2+1)` | `out: [A, Bt]`, `A: [m1, m2]`, `Bt: [B, m5]`, `B: [m3, m4]` |
| `2+2+2` | `out: [AB, C]`, `AB: [A, B]`, `A: [m1, m2]`, `B: [m3, m4]`, `C: [m5, m6]` |
| `3+3` | `out: [T1, T2]`, `T1: [A, m3]`, `A: [m1, m2]`, `T2: [B, m6]`, `B: [m4, m5]` |

**Tabulka 46.** Předvolby hierarchií. Každá dráha je zapsána jako `name: [first child, second
child]` (název: [první potomek, druhý potomek]). Trojhvězda se zapisuje jako `custom`, např.
`hierarchy: {out: [inner, m3], inner: [m1, m2]}`. Otestovány jsou stromy až do deseti hvězd.

## `orbits.<name>` {#cfg-orbits}

| Klíč | Typ, výchozí hodnota | Význam |
| --- | --- | --- |
| `P` | číslo > 0 | perioda [d]; zadejte právě jedno z `P` a `a` |
| `a` | číslo > 0 | velká poloosa [au] |
| `e` | 0 ≤ e < 1, výchozí 0 | excentricita |
| `i` | 0–180°, výchozí 0 | sklon k rovině oblohy |
| `Omega` | stupně, výchozí 0 | délka výstupného uzlu |
| `omega` | stupně, výchozí 0 | argument pericentra |
| `M0` | stupně, výchozí 0 | střední anomálie v $t = 0$ |
| `Pdot` | d/yr, výchozí 0 | vnucená sekulární změna periody dráhy dvou hvězd, $P(t) = P + \dot P\,t$ (přenos hmoty nebo ztráta momentu hybnosti polodotykových a kontaktních dvojhvězd); N-částicový backend; v proložení volný jako `orbits.<name>.Pdot` |
| `Pdot_mode` | `angular_momentum` (výchozí) nebo `conservative` | co mění periodu: ztráta momentu hybnosti při stálých hmotnostech, nebo konzervativní přenos mezi oběma hvězdami ($\dot M_1 = (\dot P/P)\,M_1 M_2 / (3(M_1 - M_2))$) |
| `contact_fillout` | 0–1, volitelný | světelné křivky: obě hvězdy této dvojhvězdy tvoří kontaktní dvojhvězdu ve společné obálce s faktorem přeplnění $f$ (Lucy 1968); jejich poloměry plynou z $f$ a mají společnou teplotu |

**Tabulka 47.** Klíče dráhy.

## `stars.<name>` {#cfg-stars}

Vše, co není zadáno, pochází z modelu hvězd.

| Klíč | Typ, výchozí hodnota | Význam |
| --- | --- | --- |
| `R` | R☉ > 0 | poloměr |
| `Teff` | K > 0 | efektivní teplota (světelné křivky; spolu s `R` určuje zářivý výkon) |
| `apsidal_k2` | 0 < k2 < 0,75 | apsidální konstanta (Loveovo číslo = 2 $k_2$) |
| `rg2` | 0 < rg2 ≤ 0,4 | moment setrvačnosti $I/(MR^2)$ |
| `lag_s` | s ≥ 0 | stálé slapové časové zpoždění (Mignard 1979); 0 = bez disipace |
| `P_rot` | dny, `sync` (výchozí) nebo `pseudosync` | rotační perioda; `sync`: střední pohyb její dvojhvězdy; `pseudosync`: pseudosynchronní rotace v excentrické dvojhvězdě (Hut 1981) |
| `obliquity` | 0–180°, výchozí 0 | úhel mezi osou rotace a normálou dráhy |
| `spin_azimuth` | stupně, výchozí 0 | azimut osy rotace kolem normály dráhy, měřený od uzlu |
| `age_offset_myr` | Myr, výchozí 0 | `sse`, `mist`: zdánlivé stáří této hvězdy minus `stellar_model.age_myr` (produkty interakce ve dvojhvězdě) |
| `pulsation` | oddíl, volitelný | pulzace ve světelných křivkách (viz níže) |
| `shape` | `sphere` (výchozí), `rotating`, `roche` | světelné křivky: koule; zploštělý rotátor s gravitačním ztemněním; Rocheova ekvipotenciála v kruhové synchronizované dvojhvězdě |
| `radius_inflation` | −0,5 až 1, výchozí 0 | relativní změna poloměru oproti modelu hvězdy při stálém zářivém výkonu (např. 0,05 pro magneticky aktivní trpaslíky typu M) |

**Tabulka 48.** Klíče hvězdy.

### `stars.<name>.pulsation` {#cfg-pulsation}

| Klíč | Typ, výchozí hodnota | Význam |
| --- | --- | --- |
| `type` | `manual` (výchozí), `lpv`, `delta_scuti`, `gamma_dor`, `beta_cep`, `spb`, `rr_lyrae`, `cepheid` | druh pulzující hvězdy; perioda u `lpv`, `delta_scuti`, `beta_cep`, `cepheid` a `rr_lyrae` plyne z publikovaných vztahů, ostatní potřebují `period_d` nebo `Q` |
| `mode` | `F` (výchozí), `1O`, `2O`, `3O` | radiální mód |
| `period_d` | dny > 0 | perioda; má přednost před vztahem |
| `Q` | dny > 0 | pulzační konstanta, $P = Q/\sqrt{\rho/\rho_\odot}$ |
| `amplitude_mag` | mag ≥ 0, výchozí 0,01 | poloviční amplituda v daném pásmu |
| `phase` | cykly, výchozí 0 | fáze v okamžiku vyzáření $t = 0$ |
| `period_spacing_s`, `n_modes` | s; 1–100, výchozí 1 | gravitační módy (`spb`, `gamma_dor`): řada `n_modes` módů $P_0 + k\,\Delta P$ |
| `extra_modes` | seznam | další módy, každý s `period_d` nebo `frequency_ratio`, `amplitude_mag` (povinný) a `phase` |

**Tabulka 49.** Klíče pulzace.

## `discs.<name>` {#cfg-discs}

Prachový disk nebo prstenec kolem hvězdy nebo dvojhvězdy, který zakrývá hvězdy za sebou: deska
s vnějším poloměrem `r_out_au`, volitelným vnitřním otvorem a celkovou tloušťkou, která zeslabuje
světlo hvězd faktorem $\exp(-\tau L/\mathrm{thickness})$ podél každého zorného paprsku.

| Klíč | Typ, výchozí hodnota | Význam |
| --- | --- | --- |
| `host` | hvězda nebo dráha, povinný | disk má střed ve hvězdě nebo v těžišti hvězd dané dráhy |
| `r_out_au`, `r_in_au` | au; `r_in_au` výchozí 0 | vnější poloměr (povinný) a vnitřní otvor |
| `thickness_au` | au, povinný | celková tloušťka |
| `tau` | > 0, výchozí žádná (neprůhledný) | optická tloušťka přes celou tloušťku, kolmo k disku |
| `i`, `Omega` | stupně, výchozí 90 a 0 | orientace roviny disku, stejně jako u drah |
| `relative_to`, `tilt`, `tilt_node` | dráha; stupně | naklonit disk vůči rovině dráhy (`i`, `Omega` se pak ignorují) |
| `precession_period_yr` | roky | s `relative_to`: precese normály disku kolem normály dráhy (záporná = retrográdní) |
| `dynamics` | `prescribed` (výchozí) nebo `rigid` | precese podle zadání, nebo spočtená pro tuhý disk, na který působí moment sil dvojhvězdy nebo průvodce jeho hostitele (vyžaduje `relative_to`) |
| `sigma_index`, `damping_yr` | 0–2,5, výchozí 1; roky | `rigid`: exponent plošné hustoty; časová škála srovnání viskózní disipací |
| `edge_au` | au, volitelný | měkký vnější okraj (vyžaduje konečné `tau`) |
| `emission` | výchozí false | přidat tepelné světlo plochého disku, který přezařuje světlo hvězd |
| `albedo`, `g` | 0–1, výchozí 0; −1 až 1, výchozí 0,5 | rozptýlené světlo hvězd: albedo jednoho rozptylu a parametr asymetrie |

**Tabulka 50.** Klíče disku. Disky ovlivňují jen světelné křivky a minima, ne dynamiku.

## `stellar_model` {#cfg-stellar-model}

| Klíč | Typ, výchozí hodnota | Význam |
| --- | --- | --- |
| `kind` | `zams` (výchozí), `sse`, `mist`, `table`, `python` | `zams`: hlavní posloupnost nulového stáří (Tout et al. 1996) a polytropy; `sse`: hvězdy ve stáří `age_myr` podle Hurleyho et al. (2000), až do zapálení helia, 0,1–8 M☉; `mist`: vývojové dráhy MIST (Choi et al. 2016) od fáze před hlavní posloupností do konce AGB; `table`: vlastní tabulka CSV; `python`: vlastní funkce |
| `age_myr` | Myr ≥ 0 | `sse`, `mist` (tam povinný): stáří hvězd v $t = 0$; `sse` počítá od ZAMS, `mist` od počátku fáze před hlavní posloupností |
| `feh` | výchozí 0 | `mist`: metalicita [Fe/H] nainstalované mřížky (−1,50 až +0,50 po krocích 0,25) |
| `vvcrit` | 0 nebo 0,4 (výchozí) | `mist`: počáteční rotace modelů |
| `z` | 0,004–0,03, výchozí 0,02 | `sse`: metalicita $Z$ |
| `structure` | `auto` (výchozí), `polytrope`, `claret2023`, `mesa:<grid>` | zdroj $k_2$ a $I/(MR^2)$: modely Clareta (2023) pro hmotnost a $\log g$ hvězdy, polytropy indexu 1,5 nebo 3, nebo vlastní mřížka MESA; `auto` = `claret2023` pro `sse`, polytropy pro `zams` |
| `file` | cesta | `table`: soubor CSV (relativně ke konfiguraci) |
| `entry` | `file.py:function` nebo `module:function` | `python`: funkce |

**Tabulka 51.** Klíče modelu hvězd. `mist` vždy používá Claretovu stavbu (nad 8 M☉ Claret 2004).

## `physics` {#cfg-physics}

| Klíč | Typ, výchozí hodnota | Význam |
| --- | --- | --- |
| `gr` | výchozí true | relativita 1PN (N-částicově), stáčení periastra (sekulárně) |
| `gw` | výchozí false | vyzařování gravitačních vln drahami dvou hvězd (Peters 1964) |
| `tides` | výchozí false | slapová deformace hvězd (apsidální stáčení; se zpožděním i disipace) |
| `rotation` | výchozí false | rotační zploštění a vazba rotace a oběhu (vyžaduje `tides`) |
| `evolution` | výchozí false | hvězdy během běhu sledují své vývojové dráhy (vyžaduje `stellar_model.kind` `sse` nebo `mist`; bez pevně zadaných `R` nebo `Teff`) |
| `tidal_lag` | `fixed` (výchozí) nebo `hurley2002` | zpoždění podle `stars.<m>.lag_s`, nebo z obálky každé hvězdy (Hurley et al. 2002; vyžaduje `tides`) |
| `tidal_efficiency` | > 0, výchozí 1 | `hurley2002`: faktor, kterým se násobí rychlost disipace |
| `tidal_convective_reduction` | `quadratic` (výchozí) nebo `linear` | `hurley2002`: snížení turbulentní viskozity pro slapové periody kratší než doba obratu konvekce |
| `k2_rotation` | výchozí false | opravit $k_2$ na rotaci hvězd (Claret 1999) |
| `inertial_waves` | výchozí false | přidat disipaci inerciálních vln v konvektivních obálkách (vyžaduje `hurley2002`) |
| `gas` | oddíl | parametrická plynná fáze (sekulární backend): `t_end_yr`, `t_decay_yr`, `plane` (`invariable` nebo `{i, Omega}`), `orbits: {<orbit>: {tau_a_yr, tau_e_yr, tau_i_yr}}` |
| `mass_loss` | oddíl | hvězdný vítr v dynamice (viz níže) |
| `mass_transfer` | oddíl | přetečení Rocheova laloku ve dvojhvězdách ze dvou hvězd (viz níže) |
| `mergers` | oddíl | pokračovat po splynutí dvou hvězd (viz níže) |
| `galactic_tide` | oddíl | slapy galaktického disku: `rho0_msun_pc3` (0,185), `pole_i_deg` (0), `pole_Omega_deg` (0) |
| `encounters` | oddíl | míjející hvězdy: `n_star_pc3` (0,1), `sigma_kms` (40), `m_perturber` (0,5), `b_max_au` (výchozí 30 × nejvzdálenější $a$), `seed` (1) |
| `activity` | zobrazení hvězda → cyklus | cyklus magnetické aktivity (Applegate 1992) v okamžicích zákrytů: `p_mod_yr`, `amplitude_s` (povinné), `phase_deg` (0) |
| `magnetic_braking` | oddíl | brzdění rotace větrem (Matt et al. 2015; vyžaduje `evolution` a `rotation`): `chi` (10), `p` (2), `t0_erg` ($9{,}5\times10^{30}$), `stars` (výchozí všechny), `max_step_yr` ($10^6$) |

**Tabulka 52.** Klíče oddílu `physics`. Výchozí hodnoty podklíčů jsou v závorkách.

### `physics.mass_loss` {#cfg-mass-loss}

| Klíč | Typ, výchozí hodnota | Význam |
| --- | --- | --- |
| `source` | `track` (výchozí) nebo `rates` | hmotnosti z vývojových drah (vyžaduje `evolution`), nebo stálé rychlosti |
| `rates` | zobrazení hvězda → M☉/yr | `rates`: $dM/dt$ (záporná = ztráta) |
| `reimers_eta` | výchozí 0,5 | hvězdy SSE: koeficient Reimersova větru |
| `masses` | `initial` (výchozí) nebo `current` | `track`: zadané hmotnosti jsou počáteční, nebo současné ve stáří `age_myr` |
| `branch` | `lowest` (výchozí) nebo `unique` | `current`: kterou počáteční hmotnost vzít, když současné hmotnosti dosahuje více vývojových drah |
| `max_psi` | výchozí 0,05 | sekulární backend: odmítnout běhy, jejichž index ztráty hmoty $\Psi = \lvert\dot M\rvert/(nM)$ tuto hodnotu překračuje |

**Tabulka 53.** Klíče oddílu `physics.mass_loss`.

### `physics.mass_transfer` {#cfg-mass-transfer}

Vyžaduje `physics.evolution: true` se `stellar_model.kind` `sse` nebo `mist`.

| Klíč | Typ, výchozí hodnota | Význam |
| --- | --- | --- |
| `rejuvenation` | `hurley` (výchozí) nebo `fractional_age` | omlazení akreujících hvězd hlavní posloupnosti |
| `stability` | `bse` (výchozí) nebo `zeta` | kritické poměry hmotností podle Hurleyho et al. (2002), nebo adiabatické exponenty dárců-obrů |
| `q_crit_ms` | výchozí 3,0 | kritický poměr hmotností pro dárce na hlavní posloupnosti a dárce s hořením helia v jádře |
| `accretion` | `bse` (výchozí) nebo `fractions` | akrece omezená tepelnou časovou škálou příjemce, nebo pevné podíly |
| `alpha`, `beta` | 0–1, výchozí 0 | `fractions`: hmota ztracená z okolí dárce a příjemce ($\alpha + \beta \leq 1$) |
| `dynamic_onset` | `angular_momentum`, `periastron` nebo žádný (výchozí) | pokračovat, když dynamika dostane hvězdu za její Rocheův lalok: zakulatit dvojhvězdu a zahájit přenos |
| `tertiary` | oddíl (výchozí přítomen) nebo null | přenos z hvězdy na dvojhvězdu, kolem níž obíhá (sekulární backend): `accretion` (`auto`), `split` (`auto`), `inner_orbit` (`energy` nebo `angular_momentum`), `alpha_lambda` (5), `h` (0) |
| `contact` | oddíl nebo žádný (výchozí) | sledovat kontaktní dvojhvězdy: `braking` (true), `q_min` (0,09) |
| `common_envelope` | oddíl nebo žádný (výchozí) | sledovat společnou obálku energetickým formalismem: `alpha` (1), `lam` (0,5), `duration_yr` (1) |

**Tabulka 54.** Klíče oddílu `physics.mass_transfer`. Předpisy jsou v dokumentu *Fyzika a metody*,
kapitola Přenos hmoty, společná obálka a splynutí.

### `physics.mergers` {#cfg-mergers}

| Klíč | Typ, výchozí hodnota | Význam |
| --- | --- | --- |
| `mass_loss` | `glebbeek2013` (výchozí) nebo podíl | hmota odvržená při splynutí |
| `rejuvenation` | `glebbeek` (výchozí), `hurley` nebo číslo | relativní stáří produktu na hlavní posloupnosti vzhledem k hmotnostmi váženým stářím obou hvězd |
| `spin_fraction` | 0–1, výchozí 0,3 | rotace produktu jako podíl kritické rotace (rozpadu) |

**Tabulka 55.** Klíče oddílu `physics.mergers`; `mergers: {}` použije výchozí hodnoty.

## `run` {#cfg-run}

| Klíč | Typ, výchozí hodnota | Význam |
| --- | --- | --- |
| `t_max_yr` | roky, povinný, ≠ 0 | konec běhu; záporný = do minulosti |
| `backend` | `nbody` (výchozí) nebo `secular` | backend ([Jak SWINGERsim funguje](#backends)) |
| `integrator` | `ias15` | N-částicový integrátor |
| `n_output` | ≥ 2, výchozí 1000 | počet ukládaných výstupů |
| `collisions` | výchozí true | zastavit, když se dvě hvězdy dotknou (N-částicově) |
| `stop_on_hierarchy_change` | výchozí false | zastavit při první potvrzené změně hierarchie |
| `hierarchy_persistence` | ≥ 1, výchozí 3 | počet výstupů, po něž musí nová hierarchie trvat, než je potvrzena |
| `stop_on_roche_overflow` | výchozí false | zastavit, když hvězda vyplní svůj Rocheův lalok |
| `stop_when_e_reached` | zobrazení dráha → e | N-částicově: zastavit, jakmile každá uvedená dráha dosáhne své excentricity |
| `max_wall_s` | sekundy | zastavit po této době výpočtu (událost `wall_time`) |
| `max_energy_error` | výchozí $10^{-5}$ | přerušit, překročí-li $\lvert\Delta E/E\rvert$ tuto hodnotu |
| `extremes_every_orbits` | ≥ 0, výchozí 10 | vnitřní vzorkování $e_\mathrm{max}$ a nejbližšího periastra, v periodách nejkratší dráhy; 0 = jen ve výstupech |
| `rate_window_yr` | roky | délka lineárních proložení pro současné rychlosti změny sklonů, na každé straně $t = 0$ (výchozí: dvě periody rušící dráhy, nejméně 6 výstupů) |
| `megno` | výchozí false | počítat MEGNO (newtonovsky; ne se ztrátou nebo přenosem hmoty) |
| `secular_octupole` | výchozí true | sekulární backend: zahrnout oktupólové členy |
| `secular_rtol` | výchozí $10^{-10}$ | sekulární backend: relativní tolerance |
| `evolution_max_dr` | 0–0,2, výchozí 0,01 | `evolution`: největší relativní změna poloměru během jednoho úseku |

**Tabulka 56.** Klíče oddílu `run`.

## `ensemble` {#cfg-ensemble}

| Klíč | Typ, výchozí hodnota | Význam |
| --- | --- | --- |
| `n` | ≥ 1, povinný | počet členů |
| `seed` | výchozí 1 | semínko losování |
| `vary` | zobrazení parametr → hodnota nebo rozdělení | losované parametry: `masses.<m>`, `orbits.<o>.<key>`, `stars.<m>.<key>` (jen číselné klíče) |
| `directions` | seznam z `forward`, `backward`; výchozí `[forward]` | integrovat každého člena do $+\lvert t_\mathrm{max}\rvert$, $-\lvert t_\mathrm{max}\rvert$ nebo obojí |
| `require` | seznam z `stable_ma01`, `stable_v22`, `hierarchical` | znovu vylosovat členy, které kritérium označí za nestabilní (`hierarchical`: každé vnitřní apocentrum leží uvnitř pericentra dráhy kolem něj) |
| `max_tries` | výchozí 1000 | povolený počet losování na člena |
| `ma01_unstable_keep`, `ma01_window` | 0–1, výchozí 1; [dolní, horní] | ponechat jen podíl losování označených jako nestabilní, nebo jen losování, jejichž nejmenší poměr MA01 leží v okně (trénovací množiny pro klasifikátory) |
| `mlp` | výchozí false | přidat ke každému členu verdikt klasifikátorů Vynatheyi et al. (2023) |
| `ghost_test`, `ghost_orbits`, `ghost_max_wall_s` | výchozí false; 100; žádný | přidat ke každému členu N-částicový test duchem |
| `save_runs` | výchozí false | uložit také soubor HDF5 každého člena |
| `observed` | zobrazení sloupec → pozorování | měření, která váží členy (viz níže) |

**Tabulka 57.** Klíče oddílu `ensemble`.

Hodnotou v `vary` je číslo (pevné), `isotropic` (sklon s rovnoměrně rozděleným $\cos i$), `angle`
(rovnoměrně 0–360°) nebo rozdělení s právě jedním z druhů z tabulky 58 a s volitelnými klíči `min`,
`max` (oříznutí: hodnota se losuje znovu) a `mirror` (sklon se s pravděpodobností ½ nahradí hodnotou
$180° - i$, pro zákrytové dvojhvězdy s neznámým smyslem oběhu).

| Druh | Tvar | Význam |
| --- | --- | --- |
| `uniform` | `[low, high]` | rovnoměrné |
| `loguniform` | `[low, high]`, kladné | rovnoměrné v logaritmu |
| `normal` | `[mean, sigma]` | gaussovské |
| `isotropic` | `[i1, i2]` | sklon s rovnoměrně rozděleným $\cos i$ mezi dvěma mezemi |
| `choice` | seznam | jedna z uvedených hodnot |

**Tabulka 58.** Rozdělení pro `ensemble.vary`.

Položka `ensemble.observed` jmenuje sloupec souhrnu ansámblu ([Výstupy](#out-ensemble)) a je jedním
z následujících:

- měření `{value, sigma}` – gaussovská váha;
- meze `{min, max}` (jedna nebo obě) – mimo ně váha 0; `missing_ok: true` nechá chybějící hodnotu
  projít;
- textová hodnota `{equals: intact}` nebo seznam hodnot – pro textové sloupce, jako je `fate`.

Název sloupce může začínat `forward:`, `backward:` nebo `central:` (průměr obou), aby se hodnota
vzala z daného běhu; výchozí je `central` pro současné rychlosti změn (`*_didt0_deg_yr`,
`*_dtiltdt0_deg_yr`) a `forward` pro vše ostatní.

## `fit` {#cfg-fit}

| Klíč | Typ, výchozí hodnota | Význam |
| --- | --- | --- |
| `minima` | cesta, povinný | tabulka minim CSV: `pair`, `type` (`p`/`s`), `time`, `sigma`, `source`, `method` (také `binary`, `kind`, `bjd_tdb`, `sigma_d`, `ref`) |
| `time_system` | `BJD_TDB` (výchozí), `HJD_UTC`, `JD_UTC` | časový systém tabulky; ostatní se převedou (vyžaduje `coordinates`; časy před rokem 1960 se berou jako světový čas, TT = UT + ΔT podle Espenaka a Meeuse 2006) |
| `coordinates` | `{ra_deg, dec_deg}` | souřadnice soustavy v ICRS |
| `ephemerides` | zobrazení pár → `{T0, P}` | lineární efemerida pro epochy a O−C (výchozí: proložená primárními minimy) |
| `exclude` | `{before, after, sources, pairs}` | vyřadit minima podle času, zdroje nebo páru |
| `error_scale` | `none` (výchozí), `global`, `per_source` | přeškálovat chyby tak, aby redukované $\chi^2$ bylo 1 |
| `rv` | oddíl | radiální rychlosti (tabulka 60) |
| `lc` | seznam | segmenty světelných křivek (tabulka 61) |
| `inclination` | seznam | sklony v epochách: `file` (sloupce `mid_hjd` nebo `time`, `i_deg`, `sigma_deg`, volitelně `limit` a `source`), `pair`, `label`, `min_sigma_deg` (0) |
| `stage1` | oddíl | analytický model O−C (tabulka 62) |
| `stage2` | oddíl | společný N-částicový model (tabulka 63) |
| `output` | oddíl | `dir` (`output/fit`), `tables` (`[csv]`, také `latex`), `plots` (`[oc, oc_components]`, také `corner`), `predict` (`years` 5, `pairs`, `start`, `samples` 100) |

**Tabulka 59.** Klíče oddílu `fit`. Cesty jsou relativní vůči konfiguraci.

| Klíč v `fit.rv` | Typ, výchozí hodnota | Význam |
| --- | --- | --- |
| `file` | cesta, povinný | CSV: `star`, `time` (BJD_TDB), `rv`, `sigma` [km/s], `source` |
| `gamma` | km/s, výchozí 0 | systémová rychlost (počáteční hodnota; volná jako `rv.gamma`) |
| `offsets` | zobrazení zdroj → km/s | posuny nulových bodů (počáteční hodnoty; volné jako `rv.offset.<source>`) |
| `stars` | seznam | použité hvězdy (výchozí všechny) |
| `max_sigma` | km/s | vynechat rychlosti s větší chybou |
| `centre` | `mass` (výchozí) nebo `light` | rychlost těžiště, nebo tokem vážená rychlost viditelného povrchu (deformace, rotace, Rossiterův–McLaughlinův efekt) |
| `band` | výchozí `V` | fotometrické pásmo pro tokové váhy při `centre: light` |

**Tabulka 60.** Klíče oddílu `fit.rv`.

| Klíč v `fit.lc[]` | Typ, výchozí hodnota | Význam |
| --- | --- | --- |
| `file` | cesta, povinný | CSV: `time` (BJD_TDB), `flux` (normovaný), `sigma`, volitelně `segment` |
| `label` | text | název ve zprávě |
| `band` | výchozí `TESS` | fotometrické pásmo |
| `exposure_min`, `n_sub` | min ≥ 0, výchozí 0; 1–15, výchozí 1 | doba expozice nebo sdružení bodů a počet modelových bodů na expozici |
| `baseline_degree` | 0–2, výchozí 1 | polynomický trend v každém segmentu |
| `scale` | výchozí true | škálovat model v každém segmentu (zředění); false: hloubky z hvězd tak, jak jsou zadány |
| `gap_d` | výchozí 0,5 | bez sloupce `segment`: rozdělit segmenty v mezerách delších než tato hodnota [d] |
| `segments` | seznam | použít jen tyto segmenty |
| `shapes` | `auto` (výchozí), `as_configured`, `spheres` | tvary hvězd: Rocheovy povrchy tam, kde hvězda vyplňuje svůj lalok nad `roche_threshold`, podle konfigurace, nebo koule |
| `roche_threshold` | 0–1, výchozí 0,5 | poměr objemového poloměru a poloměru Rocheova laloku, nad nímž `auto` použije Rocheův povrch |
| `mesh_resolution`, `coarse_resolution` | výchozí (64, 128) a (24, 48) | počty směrů sítě přesného a hrubého modelu |
| `fidelity` | `mapped` (výchozí) nebo `direct` | prokládat na hrubé síti opravené přesnou, nebo použít přesnou síť při každém vyhodnocení |
| `red_noise` | výchozí false | zvětšit chyby v každém segmentu faktorem červeného šumu β (vyžaduje `error_scale` nebo `clip_sigma`) |

**Tabulka 61.** Klíče položky `fit.lc`.

| Klíč v `fit.stage1` | Typ, výchozí hodnota | Význam |
| --- | --- | --- |
| `terms` | seznam z `lite`, `apsidal`, `quadratic`; výchozí `[lite]` | členy každého páru |
| `pairs` | zobrazení pár → členy | členy pro jednotlivé páry |
| `shared_outer` | výchozí false | jedna vnější dráha pro všechny páry |
| `P3_range_d` | [min, max] | rozsah hledání vnější periody |
| `clip_sigma`, `max_clip_fraction` | žádný; 0,1 | vyřazování odlehlých bodů |
| `sampler` | `none` (výchozí) nebo `emcee` | MCMC po metodě nejmenších čtverců |
| `walkers`, `steps`, `seed` | 48, 4000, 1 | nastavení MCMC |

**Tabulka 62.** Klíče oddílu `fit.stage1`.

| Klíč v `fit.stage2` | Typ, výchozí hodnota | Význam |
| --- | --- | --- |
| `bjd_t0` | BJD_TDB, povinný | epocha času $t = 0$ konfigurace |
| `free` | seznam, povinný | volné parametry: `orbits.<o>.<key>` s klíčem `P`, `M0`, `e`, `omega`, `Omega`, `i`, `ecosw`, `esinw` nebo `Pdot`; `masses.<m>`; `stars.<m>.<key>` s klíčem `R`, `Teff`, `apsidal_k2` nebo `rg2`, `rv.gamma`, `rv.offset.<source>` |
| `pairs` | seznam | páry, jejichž minima se prokládají (výchozí všechny) |
| `max_sigma_d` | dny | vynechat minima s větší chybou |
| `anchor` | výchozí false | nejprve ukotvit $P$ a $M_0$ každého zákrytového páru k minimům |
| `priors` | zobrazení parametr → `[mean, sigma]` nebo `{normal: [mean, sigma]}` | gaussovská apriorní rozdělení volných parametrů (na `e`, `omega` i tehdy, jsou-li volné `ecosw`, `esinw`, a naopak) |
| `bounds` | zobrazení parametr → [dolní, horní] | výchozí: $P$ v rozmezí 0,5–2 × počáteční hodnota, $e$ 0–0,95, $i$ 0–180°, hmotnosti 0,05–20 × |
| `fd_steps` | zobrazení parametr → krok | kroky konečných diferencí (výchozí: $P$ $10^{-6}$ relativně, úhly 0,01°, $e$ $10^{-4}$, hmotnosti $10^{-4}$ relativně) |
| `loss`, `f_scale` | `linear` (výchozí), `soft_l1`, `huber`, `cauchy`; 3 | robustní ztrátová funkce od `f_scale` σ |
| `error_scale` | `none`, `global`, `per_source`; výchozí `fit.error_scale` | korekce chyb po prvním proložení |
| `error_model` | `scale` (výchozí) nebo `jitter` | faktor, kterým se násobí chyby, nebo minimální chyba přičtená kvadraticky |
| `min_source_minima` | výchozí 5 | zdroje s menším počtem minim sdílejí jednu korekci |
| `clip_sigma`, `max_clip_fraction` | žádný; 0,1 | vyřazování odlehlých bodů |
| `refits` | 1–10, výchozí 4 | počet proložení za sebou, dokud se korekce nebo odlehlé body mění |
| `max_iterations` | výchozí 30 | počet iterací optimalizátoru |
| `processes` | výchozí 0 (všechna) | jádra pro výpočet derivací |
| `sampler` | `none` (výchozí) nebo `emcee` | MCMC po metodě nejmenších čtverců |
| `walkers`, `mcmc_steps` | 0 (= 2 × počet volných + 2), 300 | velikost MCMC; první polovina řetězce se zahodí |
| `mcmc_tau_target`, `mcmc_max_steps` | 0, 5000 | prodlužovat řetězec, dokud nemá délku tohoto počtu autokorelačních časů |
| `seed` | výchozí 1 | semínko MCMC |

**Tabulka 63.** Klíče oddílu `fit.stage2`. `bounds`, `fd_steps` a `priors` musejí jmenovat volné
parametry; klíče MCMC nemají se `sampler: none` žádný účinek (upozorní na to varování).

## Konfigurace `scatter` {#cfg-scatter}

Rozptylový experiment (`swingersim scatter`) používá vlastní konfiguraci s klíči `name`,
`description` a `scatter`:

| Klíč v `scatter` | Typ, výchozí hodnota | Význam |
| --- | --- | --- |
| `target`, `projectile` | těleso, povinný | osamocená hvězda nebo dvojhvězda: `masses` (jedna nebo dvě, M☉), `a_au` (dvojhvězdy), `e` (0), `e_thermal` (losovat $e$ z $f(e) = 2e$), `radii_rsun` (zapne srážky), `polytrope_n` (pro slapový záchyt) |
| `v_inf_kms` nebo `v_inf_vc` | jeden povinný | rychlost v nekonečnu [km/s] nebo v jednotkách kritické rychlosti |
| `b_au`, `b_max_au` nebo `r_peri_max_au` | jeden povinný | pevný srážkový parametr, největší srážkový parametr, nebo největší pericentrum (s gravitačním fokusováním) |
| `b_distribution` | `area` (výchozí) nebo `linear` | rovnoměrně v $b^2$ nebo v $b$ |
| `r_start_factor` | výchozí 50 | počáteční vzdálenost v jednotkách největší velké poloosy |
| `t_max_yr`, `max_wall_s`, `check_periods` | žádný, 600, 5 | meze a interval kontrol hierarchie |
| `gr`, `tidal_capture` | výchozí false | relativita 1PN; slapová ztráta energie při těsných průletech |
| `n`, `seed` | 100, 1 | počet setkání a semínko |

**Tabulka 64.** Klíče konfigurace rozptylu.

## Úplný příklad {#cfg-example}

Následující konfigurace popisuje trojhvězdu se zákrytovou vnitřní dvojhvězdou a neviděnou třetí
hvězdou. Běží tak, jak je, s příkazy `run`, `timings`, `lightcurve` a `stability` a její oddíly
`ensemble` a `fit` používají příkazy `ensemble` a `fit` (vedle ní musí existovat tabulka
`minima.csv`).

```yaml
topology: custom
hierarchy:
  out: [inner, m3]          # the inner binary and a third star
  inner: [m1, m2]
masses: {m1: 2.0, m2: 1.8, m3: 1.0}
orbits:
  inner: {P: 3.5, e: 0.0, i: 88.0, Omega: 0.0}
  out:   {P: 150.0, e: 0.3, i: 70.0, Omega: 20.0, omega: 45.0, M0: 10.0}
stars:
  m1: {R: 2.1, Teff: 9000}
  m2: {R: 1.8, Teff: 8200}
stellar_model: {kind: mist, age_myr: 100, feh: 0.0}
physics: {gr: true, tides: true, rotation: true}
run: {backend: nbody, t_max_yr: 100, n_output: 2001}
ensemble:
  n: 500
  seed: 7
  directions: [forward, backward]
  require: [stable_ma01]
  vary:
    masses.m3: {loguniform: [0.3, 3.0]}
    orbits.out.P: {loguniform: [50.0, 2000.0]}
    orbits.out.e: {uniform: [0.0, 0.6]}
    orbits.out.i: isotropic
    orbits.out.Omega: angle
  observed:
    inner_i0: {value: 88.0, sigma: 1.0}
    inner_didt0_deg_yr: {value: -0.05, sigma: 0.02}
    "backward:fate": {equals: intact}
fit:
  minima: minima.csv
  time_system: BJD_TDB
  stage2:
    bjd_t0: 2459000.0
    anchor: true
    free: [orbits.inner.P, orbits.inner.M0, orbits.out.P, orbits.out.M0,
           orbits.out.e, orbits.out.omega, masses.m3]
    error_scale: per_source
    error_model: jitter
    clip_sigma: 5.0
  output:
    dir: output/example_triple_fit
    predict: {years: 5}
```
