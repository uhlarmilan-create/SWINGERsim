# Slovník pojmů {#glossary}

Pojmy používané ve SWINGERsimu a v této příručce, vysvětlené pro čtenáře z jiného oboru astrofyziky
nebo pro studenta, který s vícenásobnými hvězdami začíná. V každé skupině jsou pojmy seřazeny
abecedně; v závorce je anglický termín, pod nímž je najdete v literatuře a v programu. Názvy
v neproporcionálním písmu jsou konfigurační klíče nebo sloupce výstupu.

## Vícenásobné hvězdy a jejich hierarchie {#gl-hierarchy}

| Pojem | Význam |
| --- | --- |
| čtyřhvězda 2+2 (2+2 quadruple) | Dvě dvojhvězdy (A a B) obíhající kolem společného těžiště po široké vnější dráze (`out`). |
| čtyřhvězda 3+1 (3+1 quadruple) | Těsná dvojhvězda (`inner`), kolem níž obíhá třetí hvězda (dráha `mid`), a kolem této trojhvězdy čtvrtá hvězda (`out`). |
| demokratická soustava (democratic system) | Tři a více hvězd bez zřetelné hierarchie, které všechny silně interagují; krátce žijící. |
| dílčí čtyřhvězda (sub-quadruple) | Čtyři hvězdy uvnitř větší hierarchie, které samy tvoří soustavu 2+2 nebo 3+1. |
| dvojitě zákrytová soustava (doubly eclipsing) | Soustava, v níž zakrývají dvě dvojhvězdy. |
| hierarchická soustava (hierarchical system) | Vícenásobná hvězda, jejíž dráhy mají velmi rozdílné velikosti, takže každá je blízká keplerovské elipse pomalu rušené ostatními. |
| pětihvězdy 2+(2+1), 2+2+1, 3+1+1 (quintuples) | Ve 2+(2+1) obíhá pátá hvězda kolem jedné z dvojhvězd; ve 2+2+1 a 3+1+1 obíhá kolem celé čtyřhvězdy. |
| předvolba (preset) | Pojmenovaná hierarchie (`2+2`, `3+1`, …), která zastupuje celý strom (tabulka 46). |
| strom hierarchie (hierarchy tree) | Popis hierarchie jako vnořených dvojic: každá dráha spojuje dva potomky, hvězdu nebo jinou dráhu. |
| šestihvězdy 2+2+2, 3+3 (sextuples) | Tři dvojhvězdy v hierarchii (např. V994 Her), nebo dvě trojhvězdy obíhající kolem sebe. |
| těžiště (centre of mass) | Hmotnostmi vážená střední poloha skupiny hvězd; v hierarchii se každá podsoustava pohybuje, jako by její hmotnost byla soustředěna v těžišti. |
| trojhvězda (triple, 2+1) | Dvojhvězda se třetí hvězdou na širší dráze. |
| vnořená trojice (nested triple) | Dráha spolu s tělesem, které ji ruší (druhou složkou její nadřazené dráhy); kritéria stability trojhvězd se uplatňují na každou vnořenou trojici. |

## Dráhy a geometrie {#gl-orbits}

| Pojem | Význam |
| --- | --- |
| argument periastra (ω, `omega`; argument of periastron) | Úhel v rovině dráhy od výstupného uzlu k periastru. |
| délka výstupného uzlu (Ω, `Omega`; longitude of the ascending node) | Úhel v rovině oblohy k bodu, kde ji dráha protíná při pohybu směrem od pozorovatele; bez astrometrie jsou měřitelné jen rozdíly Ω mezi drahami. |
| excentricita (e; eccentricity) | Tvar dráhy: 0 je kružnice, hodnoty blízké 1 velmi protáhlá elipsa. |
| invariantní rovina (invariable plane) | Rovina kolmá k celkovému momentu hybnosti soustavy; v čase se nemění. |
| Jacobiho souřadnice (Jacobi coordinates) | Souřadnice, v nichž je každá dráha relativním pohybem dvou podsoustav zastoupených svými těžišti. |
| oskulační elementy (osculating elements) | Keplerovská elipsa spočtená z okamžitých poloh a rychlostí; na časových škálách oběhu mírně kmitají. |
| periastrum, apoastrum (periastron, apoastron) | Nejbližší a nejvzdálenější bod dráhy; vzdálenost periastra $q = a(1-e)$. |
| pohled z boku, pohled shora (edge-on, face-on) | Dráha viděná ze své roviny ($i = 90°$, zákryty možné) nebo od svého pólu ($i = 0°$). |
| prográdní, retrográdní (prograde, retrograde) | Dvě dráhy obíhající ve stejném smyslu ($i_\mathrm{mut} < 90°$) nebo v opačném. |
| sklon (i; inclination) | Úhel mezi rovinou dráhy a rovinou oblohy. |
| střední (středované) elementy (mean elements) | Elementy středované přes oběžný pohyb (sekulární backend). |
| střední anomálie (M, `M0`; mean anomaly) | Fáze pohybu po dráze, rovnoměrně rostoucí s časem; `M0` je její hodnota v $t = 0$. |
| velká poloosa (a; semi-major axis) | Polovina nejdelšího průměru elipsy dráhy. |
| vzájemný sklon ($i_\mathrm{mut}$; mutual inclination) | Úhel mezi dvěma rovinami drah; klíčová veličina dynamiky. |

## Dynamika {#gl-dynamics}

| Pojem | Význam |
| --- | --- |
| 1PN, 2.5PN | Řády postnewtonovského rozvoje obecné relativity: 1PN dává stáčení periastra, 2.5PN ztrátu energie gravitačními vlnami. |
| apsidální stáčení (apsidal motion) | Otáčení přímky apsid excentrické dráhy, způsobené deformací hvězd, relativitou a dalšími tělesy; projevuje se rozbíháním primárních a sekundárních minim. |
| backend | Část SWINGERsimu, která integruje pohyb: `nbody` nebo `secular`. |
| cykly ZLK (von Zeipel–Lidov–Kozai; ZLK cycles) | Periodická výměna mezi excentricitou a sklonem vnitřní dráhy rušené vzdáleným skloněným tělesem, při vzájemných sklonech zhruba mezi 39,2° a 140,8°. |
| chaos, MEGNO | Citlivá závislost na počátečních podmínkách; MEGNO (Mean Exponential Growth factor of Nearby Orbits; Cincotta & Simó 2000) směřuje pro pravidelný pohyb k 2 a pro chaos roste. |
| N-částicová integrace (N-body integration) | Přímé numerické řešení pohybových rovnic všech hvězd, bez středování. |
| oktupól, kvadrupól (octupole, quadrupole) | Členy rozvoje interakce podle poměru velikosti vnitřní a vnější dráhy; kvadrupól je hlavní člen. |
| precese uzlů (nodal precession) | Pomalé otáčení roviny dráhy kolem celkového momentu hybnosti; mění sklon k obloze, a tím i zákryty. |
| rezonance středních pohybů (mean-motion resonance) | Téměř celočíselný poměr dvou oběžných period. |
| sekulární dynamika (secular dynamics) | Pomalý vývoj, který zbude po středování přes oběžné pohyby. |
| slapy Galaxie (Galactic tide) | Slapové pole disku Mléčné dráhy působící na velmi široké dráhy. |
| změna hierarchie (hierarchy change) | Změna toho, které hvězdy jsou vázány na které (výměna, vymrštění, rozpad, přestavba, zachycení). |

## Hvězdy, slapy a interakce {#gl-stars}

| Pojem | Význam |
| --- | --- |
| AGB (asymptotická větev obrů; asymptotic giant branch) | Pozdní stadium hvězd malých a středních hmotností s obrovskou chladnou obálkou a silným větrem. |
| Applegateův mechanismus (Applegate mechanism) | Cyklus magnetické aktivity, který mění tvar hvězdy, a tím periodu její dvojhvězdy (Applegate 1992). |
| apsidální konstanta ($k_2$; apsidal-motion constant) | Míra koncentrace hmoty hvězdy ke středu; určuje, jak silně její slapové a rotační výdutě stáčejí periastrum. |
| efektivní teplota ($T_\mathrm{eff}$; effective temperature) | Teplota černého tělesa se stejným povrchovým tokem; $L = 4\pi R^2\sigma T_\mathrm{eff}^4$. |
| fáze před hlavní posloupností (pre-main sequence) | Nejmladší stadium hvězdy, která se stále smršťuje k hlavní posloupnosti. |
| gravitační ztemnění (gravity darkening) | Rotující nebo deformovaná hvězda je jasnější tam, kde má větší povrchové tíhové zrychlení. |
| kontaktní dvojhvězda (contact binary) | Dvojhvězda, jejíž obě hvězdy přeplňují své Rocheovy laloky a sdílejí jednu obálku (hvězdy typu W UMa); faktor přeplnění $f$ od 0 (vnitřní kritická plocha) do 1 (vnější). |
| MIST, SSE, BSE | Modely hvězd: vývojové dráhy MIST (Choi et al. 2016); vzorce vývoje osamocených hvězd a dvojhvězd (Hurley et al. 2000, 2002). |
| moment setrvačnosti ($I/(MR^2)$, `rg2`; moment of inertia) | Druhá mocnina poloměru setrvačnosti hvězdy. |
| okrajové ztemnění (limb darkening) | Pokles povrchové jasnosti hvězdy od středu jejího disku k okraji. |
| omlazení (rejuvenation) | Akreující nebo splynutím vzniklá hvězda hlavní posloupnosti vmíchá do jádra čerstvý vodík a chová se jako mladší hvězda. |
| přenos hmoty (mass transfer) | Proudění plynu z hvězdy, která vyplňuje svůj Rocheův lalok (dárce), na jejího průvodce (příjemce). |
| přenos hmoty z třetí hvězdy (tertiary mass transfer) | Přenos z třetí hvězdy, která vyplní svůj Rocheův lalok, na vnitřní dvojhvězdu, kolem níž obíhá. |
| pulzující hvězda (pulsating star) | Hvězda, jejíž poloměr a jasnost se periodicky mění (δ Scuti, γ Doradus, β Cephei, SPB, RR Lyrae, cefeidy, dlouhoperiodické proměnné). |
| Rocheův lalok (Roche lobe) | Oblast kolem hvězdy ve dvojhvězdě, v níž je hmota vázána k této hvězdě; poloměr podle Eggletona (1983). |
| slapové časové zpoždění (`lag_s`; tidal time lag) | Zpoždění slapové výdutě za průvodcem; určuje slapovou disipaci (Mignard 1979; Hut 1981). |
| společná obálka (common envelope) | Nestabilní přenos hmoty, při němž plyn pohltí obě hvězdy; dráha se zmenšuje, dokud obálka není odvržena nebo hvězdy nesplynou. |
| Wolfova–Rayetova hvězda (Wolf–Rayet star) | Horká hmotná hvězda, která ztratila většinu vodíkové obálky a pohání hustý vítr. |
| ZAMS | Hlavní posloupnost nulového stáří (zero-age main sequence): počátek hoření vodíku; výchozí poloměry SWINGERsimu. |

## Stabilita a osudy {#gl-stability}

| Pojem | Význam |
| --- | --- |
| kritérium Mardlinga–Aarsetha (MA01) | Analytická mez stability hierarchické trojhvězdy (Mardling & Aarseth 2001); poměr nad 1 = stabilní. |
| kritérium Vynatheyi+22 (V22) | Zpřesněné kritérium stability trojhvězd (Vynatheya et al. 2022). |
| osud (fate) | Konečný stav běhu: neporušená soustava, změna hierarchie, srážka, splynutí, přetečení Rocheova laloku, nebo zastavení mezí. |
| poměr stability (stability ratio) | Poměr vzdálenosti periastra vnější dráhy k její kritické hodnotě podle kritéria; nad 1 stabilní. |
| test duchem (ghost test) | N-částicový test stability podle Vynatheyi et al. (2023): soustava se porovná s nepatrně posunutou kopií; rychlé rozbíhání ukazuje nestabilitu. |

## Pozorování {#gl-observations}

| Pojem | Význam |
| --- | --- |
| BJD_TDB, HJD, JD | Juliánská data: barycentrické v barycentrickém dynamickém čase (časový systém naměřených dat SWINGERsimu); heliocentrické a geocentrické v UTC (převádějí se pomocí `fit.time_system`). |
| diagram O−C, ETV | Pozorovaný minus vypočtený okamžik minima (eclipse timing variation, změny okamžiků zákrytů) vůči lineární efemeridě. |
| efekt světelného času (LTTE, LiTE; light-travel-time effect) | Zákryty přicházejí dříve nebo později, jak se dvojhvězda pohybuje po své vnější dráze (Irwin 1952). |
| efemerida (ephemeris) | Vzorec předpovídající okamžiky minim, $T = T_0 + nP$. |
| elipsoidální proměnnost (ellipsoidal variation) | Změna jasnosti dvojhvězdy se slapově deformovanými hvězdami, dvě maxima za oběh. |
| fotodynamický model (photodynamical model) | Světelná křivka spočtená přímo z poloh hvězd v N-částicové integraci. |
| fotometrické pásmo (passband) | Spektrální citlivost filtru a detektoru: UBVRI, TESS, Gaia G, BP, RP. |
| hmotnostní funkce (mass function) | Kombinace hmotností a sklonu měřená z LTTE nebo z radiálních rychlostí. |
| mez zákrytů (eclipse limit) | Nejmenší sklon, při němž pár ještě zakrývá: $\lvert\cos i\rvert < (R_1+R_2)/a$ pro kruhovou dráhu. |
| minimum (primární, sekundární; primary, secondary) | Okamžik zákrytu. Ve SWINGERsimu je primárním zákryt první hvězdy páru; pozorovatelé obvykle nazývají primárním hlubší z nich. |
| Rossiterův–McLaughlinův efekt (Rossiter–McLaughlin effect) | Zkreslení měřené radiální rychlosti během zákrytu, kdy je část rotující hvězdy skryta. |
| střed světla (centre of light) | Tokem vážená rychlost viditelného povrchu hvězdy, kterou měří spektrograf (`fit.rv.centre: light`). |
| světelná křivka (light curve) | Jasnost v závislosti na čase v daném fotometrickém pásmu. |
| třetí světlo (third light) | Světlo hvězd, které se zákrytu neúčastní; zákryt se kvůli němu jeví mělčí. |

## Prokládání a statistika {#gl-statistics}

| Pojem | Význam |
| --- | --- |
| apriorní, aposteriorní rozdělení (prior, posterior) | Rozdělení parametrů před zahrnutím měření a po něm. |
| autokorelační čas (τ; autocorrelation time) | Počet kroků MCMC, po nichž walker „zapomene“ svou polohu; řetězec má být dlouhý mnoho τ. |
| centrální rychlost změny (central rate) | Současná rychlost změny brána jako průměr proložení dopředu a dozadu kolem $t = 0$. |
| červený šum (β; red noise) | Korelovaná rezidua; β je faktor, o který se sdružená rezidua rozptylují víc, než by se rozptyloval bílý šum (Pont et al. 2006; Winn et al. 2008). |
| efektivní počet členů ($n_\mathrm{eff}$; effective number of members) | $(\sum w)^2/\sum w^2$ (Kish 1965): kolika stejně váženým členům odpovídá vážený ansámbl. |
| evidence (Z) | Střední váha ansámblu: pravděpodobnost měření za jeho hypotézy; poměry srovnávají hypotézy se stejnými apriorními rozděleními. |
| $\chi^2$, redukované $\chi^2$ (reduced $\chi^2$) | Součet čtverců normovaných reziduí; vydělené počtem stupňů volnosti je rovno zhruba 1, když model prokládá data v rámci chyb. |
| izotropní (isotropic) | Náhodná orientace s rovnoměrně rozděleným $\cos i$. |
| krok 1, krok 2 (stage 1, stage 2) | Analytický model O−C každého páru a společný N-částicový model celé soustavy v příkazu `swingersim fit`. |
| MCMC, emcee | Vzorkování aposteriorního rozdělení metodou Monte Carlo s Markovovými řetězci; emcee je afinně invariantní ansámblový vzorkovač (Goodman & Weare 2010; Foreman-Mackey et al. 2013). |
| podmínění, váhy (conditioning, weights) | Vážení členů ansámblu podle toho, jak dobře reprodukují měření, $w = \exp(-\chi^2/2)$, mimo meze 0. |
| semínko (seed) | Číslo, které určuje náhodná losování. |
| škálování chyb, jitter (error scaling) | Oprava uvedených chyb pro každý zdroj: faktor, nebo minimální chyba přičtená kvadraticky. |
| ukotvení (anchoring) | Úprava periody a fáze každého zákrytového páru před proložením v kroku 2 tak, aby modelová minima padla na naměřená. |
| vyřazování odlehlých bodů (outlier rejection) | Odstranění měření, která se odchylují o více než daný počet robustních směrodatných odchylek (`clip_sigma`). |
| Wilsonův interval (Wilson interval) | Interval spolehlivosti podílu, který se chová dobře i blízko 0 a 1 (Wilson 1927). |

## Strojové učení {#gl-ml}

| Pojem | Význam |
| --- | --- |
| emulátor (emulator) | Neuronová síť natrénovaná na ansámblu, která předpovídá jeho souhrnné sloupce z vylosovaných parametrů. |
| hluboký soubor sítí (deep ensemble) | Několik sítí natrénovaných z různých náhodných začátků; jejich rozptyl odhaduje nejistotu. |
| neuronová síť, MLP (neural network) | Pružná funkce s mnoha nastavitelnými vahami; vícevrstvý perceptron (multi-layer perceptron) je její nejjednodušší druh. |
| pokrytí, kalibrace (coverage, calibration) | Zda uvedené intervaly obsahují skutečnou hodnotu tak často, jak tvrdí; kontroluje se na odložených členech. |
| R² | Skóre regrese na odložených členech: 1 dokonalé, 0 o nic lepší než průměr. |
| SBI, NPE | Simulační inference s neuronovým odhadem aposteriorního rozdělení (simulation-based inference, neural posterior estimation): síť se z ansámblu naučí aposteriorní rozdělení parametrů. |

## Software {#gl-software}

| Pojem | Význam |
| --- | --- |
| HDF5, CSV, Parquet | Formáty souborů s daty běhů a s tabulkami. |
| IAS15 | Adaptivní integrátor vysokého řádu v REBOUND (Rein & Spiegel 2015). |
| koncept (draft) | Konfigurace napsaná ve webové aplikaci, uchovávaná v `jobs/drafts/`. |
| konfigurace (configuration) | Soubor YAML, který popisuje jednu soustavu a to, co se má počítat. |
| MESA | Kód hvězdného vývoje používaný pro vlastní mřížky (Paxton et al. 2011). |
| Numba | Překladač funkcí v Pythonu; zrychluje sekulární rovnice. |
| operátorové štěpení (operator splitting) | Střídání dvou procesů na různých časových škálách: hvězdy se vyvíjejí mezi úseky, dráhy se integrují uvnitř nich. |
| PHOEBE | Kód pro světelné křivky zákrytových dvojhvězd (Prša et al. 2016), volitelný výpočetní modul SWINGERsimu. |
| REBOUND, REBOUNDx | N-částicový kód a jeho knihovna doplňkových sil (Rein & Liu 2012; Tamayo et al. 2020). |
| úloha, zpracovatel úloh (job, worker) | Soubor úlohy popisuje jeden výpočet; zpracovatel úloh spouští soubory úloh ze složky. |
| YAML | Čitelný textový formát konfigurací a souborů úloh. |
| zpráva HTML (HTML report) | Samostatná webová stránka, kterou každý příkaz zapíše se svými výsledky. |

## Literatura {#references}

- Applegate J. H. (1992). A mechanism for orbital period modulation in close binaries. *ApJ* 385, 621.
- Choi J., Dotter A., Conroy C., Cantiello M., Paxton B., Johnson B. D. (2016). Mesa Isochrones and Stellar Tracks (MIST). I. Solar-scaled models. *ApJ* 823, 102.
- Cincotta P. M., Simó C. (2000). Simple tools to study global dynamics in non-axisymmetric galactic potentials – I. *A&AS* 147, 205.
- Claret A. (1999). *A&A* 350, 56.
- Claret A. (2004). New grids of stellar models including tidal-evolution constants up to carbon burning. I. From 0.8 to 125 M☉ at Z = 0.02. *A&A* 424, 919.
- Claret A. (2023). Theoretical tidal evolution constants for stellar models from the pre-main sequence to the white dwarf stage. *A&A* 674, A67.
- Claret A., Bloemen S. (2011). Gravity and limb-darkening coefficients for the Kepler, CoRoT, Spitzer, uvby, UBVRIJHK, and Sloan photometric systems. *A&A* 529, A75.
- Eggleton P. P. (1983). Approximations to the radii of Roche lobes. *ApJ* 268, 368.
- Foreman-Mackey D., Hogg D. W., Lang D., Goodman J. (2013). emcee: the MCMC hammer. *PASP* 125, 306.
- Goodman J., Weare J. (2010). Ensemble samplers with affine invariance. *Commun. Appl. Math. Comput. Sci.* 5, 65.
- Hamers A. S. et al. (2015). Secular dynamics of hierarchical quadruple systems: the case of a triple system orbited by a fourth body. *MNRAS* 449, 4221.
- Hamers A. S., Lai D. (2017). Secular chaotic dynamics in hierarchical quadruple systems. *MNRAS* 470, 1657.
- Hamers A. S. et al. (2021). Multiple Stellar Evolution (MSE). *MNRAS* 502, 4479.
- Hurley J. R., Pols O. R., Tout C. A. (2000). Comprehensive analytic formulae for stellar evolution as a function of mass and metallicity. *MNRAS* 315, 543.
- Hurley J. R., Tout C. A., Pols O. R. (2002). Evolution of binary stars and the effect of tides on binary populations. *MNRAS* 329, 897.
- Hut P. (1981). Tidal evolution in close binary systems. *A&A* 99, 126.
- Irwin J. B. (1952). The determination of a light-time orbit. *ApJ* 116, 211.
- Jermyn A. S. et al. (2023). Modules for Experiments in Stellar Astrophysics (MESA): time-dependent convection, energy conservation, automatic differentiation, and infrastructure. *ApJS* 265, 15.
- Kish L. (1965). *Survey Sampling*. Wiley, New York.
- Lacy C. H. S., Helt B. E., Vaz L. P. R. (1999). *AJ* 117, 541 (V907 Sco).
- Lu T., Rein H., Tamayo D., Hadden S., Mardling R., Millholland S. C., Laughlin G. (2023). Self-consistent spin, tidal, and dynamical equations of motion in the REBOUNDx framework. *ApJ* 948, 41.
- Lucy L. B. (1968). The structure of contact binaries. *ApJ* 151, 1123.
- Luo L., Katz B., Dong S. (2016). Double-averaging can fail to characterize the long-term evolution of Lidov–Kozai cycles and derivation of an analytical correction. *MNRAS* 458, 3060.
- Mardling R. A., Aarseth S. J. (2001). Tidal interactions in star cluster simulations. *MNRAS* 321, 398.
- Matt S. P. et al. (2015). The mass-dependence of angular momentum evolution in Sun-like stars. *ApJL* 799, L23.
- Mignard F. (1979). The evolution of the lunar orbit revisited, I. *Moon Planets* 20, 301.
- Mikulášek Z. (2015). Phenomenological modelling of eclipsing system light curves. *A&A* 584, A8.
- Naoz S. (2016). The eccentric Kozai–Lidov effect and its applications. *ARA&A* 54, 441.
- Naoz S., Farr W. M., Lithwick Y., Rasio F. A., Teyssandier J. (2013). Secular dynamics in hierarchical three-body systems. *MNRAS* 431, 2155.
- Paxton B. et al. (2011). Modules for Experiments in Stellar Astrophysics (MESA). *ApJS* 192, 3.
- Peters P. C. (1964). Gravitational radiation and the motion of two point masses. *Phys. Rev.* 136, B1224.
- Pont F., Zucker S., Queloz D. (2006). The effect of red noise on planetary transit detection. *MNRAS* 373, 231.
- Prša A. et al. (2016). Physics of eclipsing binaries. II. Toward the increased model fidelity. *ApJS* 227, 29.
- Reeve D. C., Howarth I. D. (2016). Non-LTE limb-darkening coefficients for hot stars. *MNRAS* 456, 1294.
- Rein H., Liu S.-F. (2012). REBOUND: an open-source multi-purpose N-body code for collisional dynamics. *A&A* 537, A128.
- Rein H., Spiegel D. S. (2015). IAS15: a fast, adaptive, high-order integrator for gravitational dynamics. *MNRAS* 446, 1424.
- Sterne T. E. (1939). Apsidal motion in binary stars. *MNRAS* 99, 451.
- Stone N. C., Leigh N. W. C. (2019). A statistical solution to the chaotic, non-hierarchical three-body problem. *Nature* 576, 406.
- Tamayo D., Rein H., Shi P., Hernandez D. M. (2020). REBOUNDx: a library for adding conservative and dissipative forces to otherwise symplectic N-body integrations. *MNRAS* 491, 2885.
- Tian Y. P., Xiang F. Y., Tao X. (2009). *Ap&SS* 319, 119.
- Tout C. A., Pols O. R., Eggleton P. P., Han Z. (1996). Zero-age main-sequence radii and luminosities as analytic functions of mass and metallicity. *MNRAS* 281, 257.
- Vokrouhlický D. (2016). Secular dynamics of gravitationally bound pair of binaries. *MNRAS* 461, 3964.
- Völschow M. et al. (2016). Eclipsing time variations in close binary systems: planetary hypothesis vs. Applegate mechanism. *A&A* 587, A34.
- Vynatheya P., Hamers A. S., Mardling R. A., Bellinger E. P. (2022). Algebraic and machine learning approach to hierarchical triple-star stability. *MNRAS* 516, 4146.
- Vynatheya P., Mardling R. A., Hamers A. S. (2023). Quadruple-star systems are not always nested triples: a machine learning approach to dynamical stability. *MNRAS* 525, 2388.
- Wilson E. B. (1927). Probable inference, the law of succession, and statistical inference. *J. Am. Stat. Assoc.* 22, 209.
- Winn J. N. et al. (2008). The transit light curve project. IX. Evidence for a smaller radius of the exoplanet XO-3b. *ApJ* 683, 1076.
- Zasche P., Vokrouhlický D., Barlow B. N., Mašek M. (2023). *AJ* 165, 81 (V907 Sco).
- Zasche P. et al. (2023b). V994 Herculis: a unique triply eclipsing sextuple star system. *MNRAS* 520, 3127.
