# Stabilita a osudy {#stability-and-fates}

Hierarchická vícenásobná hvězda je *dynamicky stabilní*, když její dráhy zůstávají velmi dlouho
dobře oddělené; nestabilní soustava se rozpadne, vymění partnery nebo vymrští hvězdu, obvykle během
počtu vnějších oběhů, který je ve srovnání s jejím stářím malý. Zda je pozorovaná konfigurace
stabilní, je proto podmínkou pro její neznámé parametry: řešení, které se rozpadne za několik tisíc
let, nemůže popisovat soustavu, která dnes existuje. SWINGERsim nabízí nástroje od levných po
rozhodující – algebraická kritéria, klasifikátory strojového učení, přímý N-částicový test, mapy
chaosu – a klasifikuje osudy běhů.

## Algebraická kritéria pro vnořené trojice {#criteria}

Každé kritérium se uplatňuje na *vnořenou trojici*: vnitřní dráhu a dráhu kolem ní, přičemž zbytek
soustavy je zredukován na hmotné body v těžištích. Čtyřhvězda 2+2 má dvě (A s B jako hmotným bodem,
B s A); 3+1 také dvě (`inner` s m3, `mid` s m4); větší stromy jednu pro každou dvojici drah. Obě
kritéria jsou zapsána jako poměr, který pro stabilní soustavy přesahuje 1.

**Mardling & Aarseth (2001)** (MA01), ve tvaru podle Vynatheyi et al. (2022, rov. 2):

$$ \frac{a_\mathrm{out}(1-e_\mathrm{out})}{a_\mathrm{in}} >
2{,}8\left[\frac{(1+q_\mathrm{out})(1+e_\mathrm{out})}{(1-e_\mathrm{out})^{1/2}}\right]^{2/5}
\left(1 - 0{,}3\,\frac{i_\mathrm{mut}}{180°}\right), \qquad
q_\mathrm{out} = \frac{m_\mathrm{out}}{m_\mathrm{in}} , $$

kde $a_\mathrm{in}$, $a_\mathrm{out}$ jsou velké poloosy vnitřní a vnější dráhy, $e_\mathrm{out}$
vnější excentricita, $m_\mathrm{in}$ hmotnost vnitřní dráhy, $m_\mathrm{out}$ hmotnost vnější složky
a $i_\mathrm{mut}$ vzájemný sklon. Levá strana je vzdálenost vnějšího periastra v jednotkách vnitřní
velké poloosy.

**Vynatheya et al. (2022)** (V22, jejich rov. 3–4) přidávají excentricitu, které vnitřní dráha
dosahuje v cyklech ZLK, a lepší závislost na sklonu:

$$\widetilde Y_\mathrm{crit} = 2{,}4\left[\frac{1+q_\mathrm{out}}{(1+\tilde e_\mathrm{in})(1-e_\mathrm{out})^{1/2}}\right]^{2/5}\left[\frac{1-0{,}2\tilde e_\mathrm{in}+e_\mathrm{out}}{8}(\cos i_\mathrm{mut}-1)+1\right],$$

kde $\widetilde Y = a_\mathrm{out}(1-e_\mathrm{out})/[a_\mathrm{in}(1+\tilde e_\mathrm{in})]$ (vnější
periastrum vůči vnitřnímu apoastru), $\tilde e_\mathrm{in} = \max(e_\mathrm{in}, \tfrac12 e_\mathrm{max}^2)$
a $e_\mathrm{max}$ je kvadrupólové maximum ZLK ([Hierarchie a dráhy](#zlk-cycles)). Vypisuje se
poměr $\widetilde Y/\widetilde Y_\mathrm{crit}$.

![**Obrázek 10.** Mez stability koplanární trojhvězdy podle Mardlinga & Aarsetha (2001). Vodorovná
osa: excentricita vnější dráhy; svislá osa: poměr vnější a vnitřní velké poloosy
$a_\mathrm{out}/a_\mathrm{in}$ (logaritmicky). Černá křivka je mez pro tři stejné hmotnosti
($q_\mathrm{out} = 0{,}5$); nad ní (zeleně) je trojhvězda stabilní po velmi dlouhou dobu, pod ní
(červeně) se rozpadne nebo si hvězdy vymění. Při $e_\mathrm{out} = 0$ je mez
$a_\mathrm{out}/a_\mathrm{in} \approx 3{,}3$; s vnější excentricitou strmě roste, protože rozhoduje
vnější periastrum, nikoli velká poloosa.](../../figures/concepts/stability.png)

Obrázek 10 ukazuje obsah kritéria: záleží na tom, jak blízko se vnější složka v periastru dostane,
takže excentrická vnější dráha musí být několikrát širší než kruhová. U skloněných drah faktor
$(1 - 0{,}3\,i_\mathrm{mut}/180°)$ mez mírně snižuje (retrográdní dráhy jsou stabilnější), kdežto tvar
V22 ji zvyšuje blízko okna ZLK, kde vnitřní excentricita roste.

**Meze.** Obě kritéria jsou aproximacemi simulací trojhvězd. U čtyřhvězd není redukce dvojhvězdy na
hmotný bod vždy spolehlivá, zejména v soustavách 3+1, kde čtvrtá hvězda budí excentricitu prostřední
dráhy (Vynatheya et al. 2023). Souhrn běhu vypisuje oba poměry pro každou dvojici drah jako první
kontrolu.

## Klasifikátory strojového učení {#mlp}

Vynatheya et al. (2023) natrénovali klasifikátory typu vícevrstvého perceptronu (MLP) na
$5\times10^5$ N-částicových integracích čtyřhvězd 2+2 a stejném počtu čtyřhvězd 3+1; správně
klasifikují 94 % a 93 % soustav, oproti 88 % a 66 % pro jejich MLP pro trojhvězdy (Vynatheya et al.
2022) použitý na vnořené trojice a 83 % a 56 % pro kritérium Mardlinga & Aarsetha (2001) (Vynatheya et
al. 2023, tab. 3–4). SWINGERsim používá sítě publikované autory, je-li nainstalován jejich repozitář.
Hierarchie s více než čtyřmi hvězdami se klasifikují po jednotlivých dílčích čtyřhvězdách: každá
dráha, jejíž obě složky jsou dráhy, dává 2+2, každá dráha uvnitř dráhy uvnitř dráhy 3+1, přičemž
hlubší podsoustavy jsou hmotnými body ve svých těžištích (2+2+2 dává čtyři dílčí čtyřhvězdy).
Soustava se počítá za stabilní, jsou-li stabilní všechny; mimo trénovací množinu je to aproximace.
SWINGERsim umí také natrénovat vlastní klasifikátory na ansámblech
([Ansámbly a inference](#emulators)).

## N-částicový test duchem {#ghost-test}

Rozhodující odpovědí pro jednu konfiguraci je přímá integrace s definicí stability podle Vynatheyi et
al. (2023, odd. 2). Soustava a její kopie, *duch*, v níž je velká poloosa méně hmotné vnitřní
dvojhvězdy (2+2) nebo nejvnitřnější dvojhvězdy (3+1) větší o $10^{-6}$ au, se integrují po 100 vnějších
oběhů jako newtonovské hmotné body. Soustava je nestabilní, stane-li se některá dráha nevázanou (alespoň
ve třech po sobě jdoucích kontrolách) nebo překročí-li relativní rozdíl obou velkých poloos $10^{-2}$ –
rozbíhání sousedních trajektorií, které signalizuje chaos. Test trvá sekundy až minutu.

**Meze.** Všechny tři nástroje se týkají dynamické stability po dobu zhruba 100 vnějších oběhů.
Soustava, která testem projde, se může na delších časových škálách dál sekulárně vyvíjet (cykly ZLK,
chaotická difuze); to kontroluje MEGNO a dlouhé běhy. Kritérium s duchem je citlivé na
krátkoperiodické kmity oskulační velké poloosy: u středně hierarchických soustav se obě kopie, jejichž
fáze se rozcházejí, liší asi o 1 %, i když obě zůstanou vázané po 1000 vnějších oběhů. Verdikty
„rozešlo se“ s rozdílem jen mírně nad $10^{-2}$ jsou proto hraniční.

## Mapy chaosu {#chaos-maps}

`swingersim chaosmap` integruje mřížku konfigurací, v nichž se mění jeden nebo dva parametry
(například `orbits.out.a` a `orbits.out.e`), a na konci každého běhu zaznamená MEGNO
([Dynamika](#megno)). Bod se počítá za chaotický, když MEGNO ≥ 4; bod, v němž hvězda unikne (za
dvacetinásobek největší počáteční velké poloosy) nebo se dvě hvězdy srazí, se označí zvlášť. Běhy
počítají s newtonovskými hmotnými body, protože variační rovnice REBOUND pokrývají jen newtonovskou
gravitaci.

![**Obrázek 11.** Mapa chaosu MEGNO pro testovací částici kolem kruhové dvojhvězdy se stejnými
hmotnostmi a vzdáleností složek 1 au, po 3000 letech. Vodorovná osa: velká poloosa dráhy částice
(`orbits.out.a`, au, tj. v jednotkách vzdálenosti složek dvojhvězdy); svislá osa: její excentricita
(`orbits.out.e`). Barva: MEGNO, od modré (2, pravidelný pohyb) přes bílou (4) po tmavě červenou (10,
použito pro únik). Modrá oblast ve velkých vzdálenostech a při malých excentricitách je pravidelná;
červená oblast blízko dvojhvězdy a při vysoké excentricitě je chaotická nebo nestabilní; hranice
s excentricitou stoupá, jak předpovídají kritéria stability.](../../figures/chaosmap_circumbinary.png)

Obrázek 11 se čte jako obrázek 10 s prohozenými osami: stabilní dráhy leží daleko od dvojhvězdy a jsou
téměř kruhové a hranice se s rostoucí excentricitou posouvá ven. Osamocené červené buňky uvnitř modré
oblasti a modré buňky uvnitř červené jsou chaotické ostrovy a rezonance, které jediné kritérium
nedokáže vystihnout. Je-li druhým parametrem počáteční fáze (`orbits.<orbit>.M0`), je nejmenší
vzdálenost, za níž je každá fáze pravidelná, mezí stability: pro testovací částici kolem kruhové
dvojhvězdy se stejnými hmotnostmi nachází SWINGERsim 2,30 vzdálenosti složek oproti 2,27
v publikované mřížce Quarlese et al. (2018) (viz *Testování, část I*).

## Osudy a jejich klasifikace {#fates}

Každý běh i každý člen ansámblu končí *osudem*. Uvádí je tabulka 9.

**Tabulka 9.** Osudy běhu. *Osud* je hodnota uložená ve výstupech (sloupec `fate` ansámblů);
*Význam* co se stalo; *Rozpoznáno pomocí* mechanismus, který jej rozpozná.

| Osud | Význam | Rozpoznáno pomocí |
|---|---|---|
| `intact` | hierarchie přežila celý běh | žádná událost |
| `exchange` | hvězda vyměnila partnera | rozpoznání hierarchie |
| `reconfiguration` | strom změnil tvar, aniž by ztratil hvězdu | rozpoznání hierarchie |
| `ejection` | hvězda opustila soustavu | rozpoznání hierarchie |
| `breakup` | soustava se rozdělila na nevázané části | rozpoznání hierarchie |
| `capture` | vznikla nová vázaná struktura | rozpoznání hierarchie |
| `hierarchy_change` | jiná potvrzená změna stromu | rozpoznání hierarchie |
| `collision` | dotkly se povrchy dvou hvězd | hledání srážek v REBOUND |
| `merger` | dvě hvězdy splynuly (společná obálka, kontaktní dvojhvězda) | historie přenosu hmoty |
| `roche_overflow` | hvězda vyplnila svůj Rocheův lalok a přenos se nesledoval | faktor vyplnění Rocheova laloku |
| `e_reached` | bylo dosaženo požadované excentricity | `run.stop_when_e_reached` |
| `energy_error`, `wall_time`, `gw_unresolved`, `error` | běh se zastavil z numerického nebo výpočetního důvodu | řízení běhu |

Změna hierarchie je potvrzena až poté, co přetrvala `run.hierarchy_persistence` výstupů
([Dynamika](#hierarchy-detection)). Poslední řádek není fyzika: takové členy je třeba počítat zvlášť
a jsou-li časté, změnit nastavení.

## Přežití v obou směrech času {#survival}

Soustava pozorovaná dnes existuje od svého vzniku, takže konfigurace, která se rozpadne v nedávné
minulosti, ji popisovat nemůže. Protože dynamika bez disipace je vratná v čase, platí totéž pro
blízkou budoucnost: nestabilní konfigurace je nestabilní v obou směrech. V N-částicovém ansámblu
CzeV4315 o 10 000 členech integrovaném na ±60 let se 6,0 % členů rozpadne dopředu a stejně mnoho zpět
a dvě třetiny těch, které selžou v jednom směru, selžou během týchž 60 let i ve druhém; ostatní ne jen
proto, že okno je krátké. Jednostranné přežití v krátkém běhu proto není důkazem stability. Přežití se
v ansámblech ukládá jako podmínka na oba běhy, `backward:fate: {equals: intact}`
a `forward:fate: {equals: intact}` ([Ansámbly a inference](#weights)).

Běh v délce desetiletí zachytí jen nejrychlejší nestability. Blízko hranice stability roste doba do
rozpadu zhruba exponenciálně s $(r_\mathrm{p,out}/a_\mathrm{in})^{3/2}$ ($r_\mathrm{p,out}$ je
vzdálenost vnějšího periastra) a mezi sousedními počátečními fázemi se rozptyluje o jeden až dva řády,
což je známka chaosu (Hayashi et al. 2022, 2023; Mushkin & Katz 2020); u silně skloněných drah
(60°–150°) ji oscilace ZLK dále zkracují a nejstabilnější jsou retrográdní koplanární soustavy.
Konfigurace slučitelná s pozorovanou soustavou musí přežít její stáří – obvykle miliony vnějších
oběhů –, takže kritéria a MEGNO jsou doplňujícími testy pro členy blízko hranice.

**Co mění stabilitu.** Bez disipace a při stálých hmotnostech se stabilita hierarchie sama nemění;
přes hranici ji posouvají pomalé procesy (tabulka 10).

**Tabulka 10.** Procesy, které posouvají hierarchii přes hranici stability. *Destabilizující*
a *Stabilizující* uvádějí procesy s literaturou; řádky nejsou párovány.

| Destabilizující | Stabilizující |
|---|---|
| ztráta hmoty vnitřní dvojhvězdy (větry, supernova): hierarchie slábne (Perets & Kratter 2012; Toonen et al. 2022) | slapové zmenšování vnitřní dvojhvězdy, také cykly ZLK se slapovým třením (Fabrycky & Tremaine 2007) |
| rozpínání hvězdy: přetečení Rocheova laloku, přenos hmoty, společná obálka | precese způsobená slapy, rotací a relativitou, která tlumí cykly ZLK |
| cykly ZLK vyhánějící vysokou excentricitu | disipace v plynu během vzniku (utlumené excentricity, srovnané dráhy) |
| setkání s hvězdami a slapy Galaxie na širokých drahách (Kaib et al. 2013; Grishin & Perets 2022) | akrece a migrace v protoplanetárním disku vzniku (Tokovinin & Moe 2020) |
| průchod rezonancemi při změnách period | retrográdní koplanární konfigurace (Hayashi et al. 2022) |

Všechny procesy z tabulky kromě supernov jsou ve SWINGERsimu modelovány; tlumení plynem je
parametrické.

## Rozptylové experimenty {#scattering}

`swingersim scatter` studuje dynamický vznik: dvojhvězda (nebo osamocená hvězda) se setká s jinou
dvojhvězdou (nebo osamocenou hvězdou). Každé setkání se nastaví jako u Huta & Bahcalla (1983) a ve
FEWBODY (Fregeau et al. 2004): izotropní orientace a náhodné fáze dvojhvězd, pevné nebo tepelné
excentricity ($f(e) = 2e$), těžiště na hyperbole s rychlostí v nekonečnu $v_\infty$ a srážkovým
parametrem $b$ ($b^2$ rovnoměrně do $b_\mathrm{max}$; s maximálním periastrem $q$ plyne z gravitační
fokusace $b_\mathrm{max}^2 = q^2 + 2GMq/v_\infty^2$). Kritická rychlost $v_c$ je ta, při níž je celková
energie nulová, $\tfrac12\mu_\infty v_c^2 = \sum G m_1 m_2/(2a)$ ($\mu_\infty$ je redukovaná hmotnost
obou přilétajících objektů, součet se bere přes dvojhvězdy). Integrace (IAS15, hmotné body, volitelně
1PN, srážky, jsou-li zadány poloměry) začíná v násobku největší velké poloosy a končí, když jsou
objekty nejvyšší úrovně nalezené rozpoznáním hierarchie po dvojicích nevázané, vzdalují se a jsou od
sebe dál než počáteční vzdálenost, přičemž každá vázaná podsoustava tří a více hvězd je stabilní podle
kritéria MA01. Soustava, která zůstane vázaná, se počítá za zachycenou, jakmile je stabilní po
10 svých vnějších oběhů; jinak se běh zastaví na svém časovém nebo výpočetním limitu jako
*nerozhodnutý*.

Výsledky se řídí Antogninim & Thompsonem (2016): *průlet* (flyby), *výměna* (exchange), *jednoduchá
ionizace*, *výměna + ionizace*, *trojhvězda + osamocená hvězda*, *úplná ionizace*, *vázaná 2+2* nebo
*vázaná 3+1* (zachycená čtyřhvězda) a *srážka*. Účinné průřezy jsou $\sigma = \pi b_\mathrm{max}^2
\times$ podíl daného výsledku. Při slapovém zachycení odebere každý těsný průchod periastrem dvou
hvězd energii slapů, které vybudí (Press & Teukolsky 1977),
$E = \sum G M_\mathrm{o}^2/R\,[(R/R_p)^6 T_2(\eta) + (R/R_p)^8 T_3(\eta)]$ s
$\eta = (M/(M+M_\mathrm{o}))^{1/2}(R_p/R)^{3/2}$, kde $M$, $R$ jsou hmotnost a poloměr rušené hvězdy,
$M_\mathrm{o}$ hmotnost druhé hvězdy, $R_p$ vzdálenost periastra a $T_2$, $T_3$ aproximace pro
polytropy podle Portegiese Zwarta & Meinena (1993). Relativní rychlost se v periastru sníží, takže
nevázané hvězdy mohou být zachyceny do dvojhvězd: dvě hvězdy 0,6 M☉ ($n = 3$, 0,65 R☉) setkávající se
rychlostí 10 km s⁻¹ jsou zachyceny do 2,37 poloměru hvězdy, hvězdy s $n = 1{,}5$ do 3,33 poloměru, blízko
hodnotě zhruba 3,4 podle McMillana et al. (1987). Odpor plynu během setkání se nemodeluje.

## Statistická teorie úniku {#escape-theory}

Když tři hvězdy interagují rezonančně (bez hierarchie, v „demokratické“ fázi), je výsledek prakticky
náhodný a Stone & Leigh (2019) ukázali, že jeho rozdělení je dáno objemem fázového prostoru, který
připouštějí zachovávající se energie a moment hybnosti (jejich rov. 2, s poloměrem silné interakce
$R = 2a_B(1+e_B)$, kde $a_B$ a $e_B$ jsou elementy výsledné dvojhvězdy). `swingersim escape-theory`
vyhodnotí tento integrál numericky pro každou ze tří hvězd, které mohou uniknout, a udá
pravděpodobnosti úniku a rozdělení energie, excentricity a orientace dvojhvězdy, za sekundu místo
tisíců N-částicových běhů. Reprodukuje limitní případy z článku (rozdělení energie
$\propto |E_B|^{-7/2}$ pro nulový moment hybnosti a $|E_B|^{-4}$ blízko maxima; jejich rozdělení
excentricity) a souhlasí s vlastními rozptylovými experimenty SWINGERsimu: při rezonančních setkáních
dvojhvězdy 1 + 0,6 M☉ s hvězdou 0,8 M☉ uniká nejlehčí hvězda v 75 % N-částicových běhů a v teorii
v 77 %. Teorie pokrývá jen rezonanční setkání tří těles; okamžité výměny a průlety leží mimo ni.

## Implementace {#implementation-stability}

**Externí kód a data.** Klasifikátory MLP jsou natrénované sítě Vynatheyi et al. (2023), použité tak,
jak byly publikovány; N-částicové integrace testu duchem, map chaosu a rozptylových experimentů
používají REBOUND/IAS15 a jeho MEGNO (Rein & Liu 2012; Rein & Spiegel 2015). **Vlastní kód.** Kritéria
MA01 a V22 pro každou vnořenou trojici libovolného stromu; rozklad velkých hierarchií na dílčí
čtyřhvězdy a výpočet vstupních parametrů klasifikátorů; test duchem; mapy chaosu a jejich
klasifikace; klasifikace osudů; nastavení, ukončení a klasifikace výsledků rozptylových experimentů;
ztráta energie při slapovém zachycení s aproximacemi Portegiese Zwarta & Meinena (1993); a numerické
vyhodnocení integrálu fázového prostoru podle Stonea & Leighe (2019).
