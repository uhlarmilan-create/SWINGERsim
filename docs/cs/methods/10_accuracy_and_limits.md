# Přesnost a meze platnosti {#accuracy-and-limits}

Tato kapitola shrnuje, co předchozí kapitoly říkají o přesnosti a rozsahu platnosti: jak přesně
numerické metody řeší své rovnice, která fyzika se nemodeluje, v jakých režimech přestávají modely
platit a jak byla ověřena každá fyzikální složka.

## Numerická přesnost {#numerical-accuracy}

Tabulka 15 uvádí numerickou přesnost hlavních metod – chybu, s jakou SWINGERsim řeší své vlastní
rovnice, nikoli chybu fyziky.

**Tabulka 15.** Numerická přesnost. *Metoda* je krok výpočtu; *Míra* porovnávaná veličina; *Přesnost*
dosažená úroveň; *Kapitola* místo, kde je metoda popsána. Zkratky: HP hlavní posloupnost,
HM Hertzsprungova mezera, VO větev obrů.

| Metoda | Míra | Přesnost | Kapitola |
|---|---|---|---|
| elementy → stav → elementy | elementy po zpětném převodu | $10^{-12}$ | Hierarchie a dráhy |
| N-částicová integrace (IAS15) | relativní chyba energie, newtonovská | ~$10^{-14}$ | Dynamika |
| N-částicová integrace s 1PN | 1PN hamiltonián | $10^{-13}$–$10^{-11}$ | Dynamika |
| běh dopředu a zpět | obnovené počáteční elementy | $10^{-12}$ | Dynamika |
| sekulární backendy | relativní tolerance DOP853; kontrola zachování $H$ | $10^{-10}$ (výchozí) | Dynamika |
| interpolace SSE | poloměr oproti přímému výpočtu SSE (95. percentil) | 0,8 % HP, 1,7 % HM, 4 % VO | Hvězdy |
| interpolace MIST | poloměr na konci HP, věky (dvojnásobný krok mřížky) | 1 %, 1–2 % | Hvězdy |
| zákryty koulí | tok oproti Mandelovi & Agolovi (2002) | $10^{-8}$ | Pozorovatelné veličiny |
| povrchové sítě (64 × 128) | tok oproti přesným zákrytům koulí | $3\times10^{-4}$ | Pozorovatelné veličiny |
| zákryt diskem (4000 bodů) | tok na ostrém okraji | ~$10^{-3}$ | Pozorovatelné veličiny |
| syntetická minima v proložení | okamžiky oproti hustě vzorkovaným okamžikům | 0,06 s | Prokládání pozorovaných dat |
| předepsaná změna periody | odchylka od paraboly | < 0,05 s | Přenos hmoty a vývoj |
| radiální rychlosti | nesoulad okamžiku vyslání | < 1 ms | Pozorovatelné veličiny |

Numerická přesnost je zřídka omezujícím faktorem. Meze plynou z fyziky a z dat.

## Meze režimů {#regime-limits}

- **Sekulární backendy** vyžadují silné hierarchie ($\varepsilon_\mathrm{SA}$ hluboko pod 0,01), žádné
  rezonance středních pohybů a žádná těsná přiblížení; v chaotických oblastech se s N-částicovým řešením
  shodují jen omezenou dobu ([Dynamika](#secular-validity)).
- **Chaos** omezuje každý jednotlivý běh na několik Ljapunovových časů, dopředu i zpět; za touto hranicí
  mají smysl jen ansámbly ([Dynamika](#megno)).
- **Modely hvězd**: SSE končí zapálením helia a pokrývá 0,1–8 M☉; MIST pokrývá 0,1–300 M☉ až po opuštění
  AGB (hmotné hvězdy: hoření uhlíku); pozůstatky se nevyvíjejí; hmotnosti a poloměry vyvinutých hmotných
  hvězd závisejí na modelu na úrovni desítek procent ([Hvězdy](#mist-model)).
- **Strukturní konstanty** z polytrop jsou hrubé (pro hvězdy 2–3 M☉ asi trojnásobně); Claretovy mřížky
  pokrývají 0,2–8 M☉ (Claret 2023) a až 125 M☉ při jedné metalicitě (Claret 2004).
- **Slapy**: disipace rovnovážného slapu je nejistá o velké faktory; cirkularizace na hlavní posloupnosti
  ve starých hvězdokupách není reprodukována ([Slapy a rotace](#dissipation-from-structure)).
- **Přenos hmoty** se řídí předpisy BSE pro kruhové dvojhvězdy; heliové hvězdy a akrece větru se
  nemodelují ([Přenos hmoty a vývoj](#mass-transfer-limits)).
- **Světelné křivky**: Rocheovy tvary jen pro kruhové synchronizované dvojhvězdy; členy malých deformací
  v modelu koulí se mýlí o 3–17 % v amplitudě elipsoidální proměnnosti při $R/a$ = 0,10–0,25 a pro
  těsné páry je třeba je nahradit sítěmi ([Pozorovatelné veličiny](#distorted-stars)).
- **Nástroje stability** se týkají zhruba 100 vnějších oběhů; klasifikátory MLP jsou mimo svou
  trénovací množinu aproximacemi ([Stabilita a osudy](#ghost-test)).
- **Proložení** nacházejí lokální minima; jejich lineární chyby platí v přímých údolích $\chi^2$ a chyby
  z MCMC jen pro zkonvergované řetězce ([Prokládání](#fitting-limits)).

## Co se nemodeluje {#not-modelled}

Tabulka 16 uvádí efekty, které leží mimo rozsah SWINGERsimu, a čím jsou nahrazeny, pokud vůbec.

**Tabulka 16.** Efekty mimo rozsah SWINGERsimu. *Efekt* je fyzika; *Důsledek* co jeho absence znamená
pro výsledky a co jej nahrazuje, pokud něco.

| Efekt | Důsledek |
|---|---|
| hydrodynamika (proudění plynu, hydrodynamika společné obálky, srážky jako děje v tekutině) | plyn vstupuje jen parametricky (tlumení drah během vzniku); společné obálky a splynutí přes energetické a hmotnostní předpisy |
| magnetohydrodynamika (dynama, rozlišená magnetická pole) | magnetismus vstupuje přes parametrizované zákony brzdění a aktivity kalibrované na MHD a pozorováních |
| odpor plynu v rozptylových experimentech | slapové zachycení je zahrnuto, odpor nikoli |
| excentrické nebo asynchronní deformované dvojhvězdy ve světelných křivkách ze sítí | Rocheova geometrie předpokládá kruhové synchronizované dráhy; excentrické deformované dvojhvězdy zvládá hybrid s PHOEBE |
| hvězdy po AGB, bílí trpaslíci, neutronové hvězdy, černé díry jako vyvíjející se hvězdy | mohou být hmotnými body; jádro po společné obálce je hmotným bodem stálé hmotnosti |
| heliové hvězdy, akrece větru | dárce obnažený na jádro ukončí historii; větry jsou izotropní a opouštějí soustavu |
| slapová disipace na hlavní posloupnosti nad rámec rovnovážného slapu a dvouvrstvého odhadu inerciálních vln | periody cirkularizace na hlavní posloupnosti ve starých hvězdokupách vycházejí příliš krátké (M67: 7,5 d místo 12,1 d; Barker 2022 vysvětluje M67 inerciálními vlnami spočtenými z úplných profilů hvězd) |
| hvězdné skvrny ve světelných křivkách | cykly aktivity vstupují do okamžiků zákrytů, skvrny nikoli do světla |
| tvary pulzací, amplitudy ze stavby hvězdy, rozestupy gravitačních módů ze stavby | módy jsou sinusoidy pevné amplitudy; hvězdy SPB a γ Dor potřebují $P_0$ a $\Delta P$ od uživatele |
| hmotnost a gravitace disku, zborcení, dopředný rozptyl v tenkých částech | disky jsou nehmotné desky nebo tuze precedující roviny |
| astrometrie (pohyb fotocentra, astrometrické dráhy z Gaie) | polohy jsou ve výstupu, žádná astrometrická pozorovatelná veličina se však nepočítá ani neprokládá |
| supernovy a impulzy při vzniku pozůstatku (natal kicks) | nemodelují se |

## Princip ověřování {#validation-principle}

Každý fyzikální modul SWINGERsimu se před použitím porovná s publikovanými výsledky: s analytickými
vzorci, publikovanými numerickými výsledky, jinými kódy (příklady REBOUND/REBOUNDx, MSE, PHOEBE,
BSE/COSMIC, MESA, batman, ellc) a pozorovanými soustavami. Každé srovnání je jedno číslo s tolerancí
stanovenou před srovnáním a jeho verdikt se zaznamená: vyhověl, vysvětleno (mimo toleranci
z dokumentovaného důvodu, který není chybou SWINGERsimu – například publikovaná čísla, která si
navzájem odporují, nebo uvedené omezení fyzikálního modelu), nebo nevyhověl. Tolerance se nerozšiřuje,
aby zakryla nesoulad. Jednovětná shrnutí v předchozích kapitolách („reprodukuje … na …“) jsou výsledky
těchto srovnání; samotná srovnání s jejich zdroji, vstupy a výsledky dokumentuje svazek *Testování
SWINGERsimu, část I: Ověření proti literatuře*. Části II–IV testovacích dokumentů uplatňují celý
řetězec na publikované vícenásobné hvězdy, na zákrytové dvojhvězdy s proměnným sklonem a na zvláštní
případy.

## Implementace {#implementation-accuracy}

Údaje o přesnosti v tabulce 15 měří vlastní sada testů a sbírka srovnávacích úloh SWINGERsimu
(reprodukce publikovaných výsledků), které jsou součástí kódu. Integrátory, jejichž zaokrouhlovací
chování určuje dolní mez přesnosti, jsou externí: IAS15 z REBOUND (Rein & Spiegel 2015) a DOP853
ze SciPy.
