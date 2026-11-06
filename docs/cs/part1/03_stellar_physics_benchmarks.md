# Hvězdy: stavba, slapy, vývoj a výměna hmoty {#stellar}

Hvězdy vícenásobné soustavy nejsou body. Mají velikost a vnitřní stavbu, deformují je slapy
a rotace, vyvíjejí se, ztrácejí hmotu hvězdným větrem a vyměňují si hmotu se svými průvodci. Tato
kapitola ověřuje každý modul SWINGERsimu, který to popisuje, od strukturních konstant jednotlivé
hvězdy až po dvojhvězdu, která se vyvine v kontaktní a splyne.

## Strukturní konstanty, Rocheovy laloky a rotující hvězdy {#structure}

**Definice.** **Apsidální konstanta** $k_2$ měří, jak silně se hvězda deformuje slapovou nebo
odstředivou silou; spolu s poloměrem $R$ určuje dodatečné apsidální stáčení excentrické
dvojhvězdy, které škáluje jako $k_2R^5$. **Poloměr setrvačnosti** $\beta$ (neboli
$r_g^2 = I/MR^2$) určuje moment setrvačnosti $I$. **Rocheův lalok** je oblast kolem hvězdy ve
dvojhvězdě, v níž převládá její vlastní gravitace; jeho poloměr ekvivalentního objemu $R_L$ udává
s přesností lepší než 1 % vzorec Eggletona (1983). Hvězda rotující tuze s úhlovou rychlostí
$\Omega$ má v Rocheově modelu rovníkový poloměr, který s $\Omega$ roste a při kritické rotaci
dosáhne 1,5násobku polárního poloměru.

**Reference.** Klasické apsidální konstanty polytrop indexu $n$ = 1,5 a 3 (Brooker & Olle 1955);
publikované tabulky $\log k_2$ a $\beta$ podél vývojových drah hvězd podle Clareta (2004, do
125 M☉) a Clareta (2023); přesný Rocheův potenciál, integrovaný přímo testem, pro Eggletonův vzorec
a pro povrchovou síť SWINGERsimu; uzavřené tvary Rocheova modelu rotující hvězdy; a rotační oprava
$k_2$, $\delta\log k_2 = -0{,}87\,\lambda_s$ s $\lambda_s = 2V^2/(3gR)$, tedy faktor
$10^{-0{,}87\cdot 2V^2/3gR}$ (Claret 1999).

**Tabulka 14.** Strukturní konstanty, Rocheovy laloky a rotující hvězdy. Sloupce jsou definovány v oddílu [Jak číst tuto zprávu](#how-to-read).

| Test | Třída | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|---|
| `polytrope-k2-n1.5` | own | 0,1433 | 0,143279 | −0,01 % | 0,2 % | PASS |
| `polytrope-k2-n3` | own | 0,0144 | 0,014443 | +0,30 % | 0,5 % | PASS |
| `claret2023-k2-2msun` | own/tab | −2,478 | −2,47794 | +6·10⁻⁵ | 0,002 | PASS |
| `claret2023-beta-1msun` | own/tab | 0,2265 | 0,226499 | −6,3·10⁻⁷ | 0,0005 | PASS |
| `claret2004-log-k2` | own/tab | −2,81523 | −2,81569 | −0,000463 | 0,002 | PASS |
| `claret2004-beta` | own/tab | 0,161706 | 0,161699 | −6,6·10⁻⁶ | 0,001 | PASS |
| `k2-rotation` | own | 1 | 1 | +0,00 % | 10⁻⁷ % | PASS |
| `roche-eggleton-q1` | own | 0,379861 | 0,378921 | −0,25 % | 1 % | PASS |
| `roche-eggleton-q0.1` | own | 0,2054 | 0,206773 | +0,67 % | 1 % | PASS |
| `roche-lobe-volume` | own | 1 | 1,00509 | +0,51 % | 1 % | PASS |
| `rotating-star-critical` | own | 1,5 | 1,5 | −0,00 % | 10⁻⁴ % | PASS |
| `rotating-star-equator` | own | 0 | −2,2·10⁻¹⁶ | −2,2·10⁻¹⁶ | 10⁻⁹ | PASS |

**Výsledek (tabulka 14).** SWINGERsim pro svůj výchozí model hvězdy integruje Laneovu–Emdenovu
a Radauovu rovnici; apsidální konstanty polytrop souhlasí s Brookerem & Ollem (1955) na 0,01 %
a 0,3 % (automatické testy ověřují také analytické hodnoty pro $n$ = 1 na 10⁻⁸). Interpolace
Claretových tabulek v $\log g$ podél vývojových drah reprodukuje publikované řádky do poslední
číslice (6·10⁻⁵ v $\log k_2$, 6·10⁻⁷ v $\beta$). Eggletonův vzorec, který SWINGERsim používá pro
všechny faktory vyplnění Rocheova laloku, se od přesného poloměru ekvivalentního objemu liší
o −0,25 % při $q$ = 1 a +0,67 % při $q$ = 0,1, v rámci 1 %, které Eggleton uvádí pro všechny
poměry hmotností. Povrchová síť hvězdy vyplňující Rocheův lalok (používaná pro světelné křivky
deformovaných hvězd, [kapitola 4](#meshes)) uzavírá objem, jehož ekvivalentní poloměr je o 0,51 %
větší než Eggletonův; zbytek odpovídá rozlišení sítě.

## Apsidální stáčení a slapová disipace {#tides}

**Definice.** Slapy mají konzervativní část – statickou výduť, která přidává apsidální stáčení –
a disipativní část – zpoždění výdutě za spojnicí hvězd, které přenáší moment hybnosti mezi drahou
a rotacemi a disipuje energii. N-částicový backend SWINGERsimu používá efekt `tides_spin` knihovny
REBOUNDx (rovnovážný slap s konstantním časovým zpožděním podle Eggletona et al. 1998; Lu et al.
2023) s parametry hvězd nastavenými SWINGERsimem; sekulární backendy používají vlastní rovnice
SWINGERsimu středované přes oběh (Hut 1981).

**Reference.** Apsidální perioda ze slapové a rotační deformace obou hvězd (Sterne 1939); rychlosti
$\dot a$, $\dot e$ a $\dot\Omega$ modelu s konstantním časovým zpožděním středované přes oběh
(Hut 1981; Leconte et al. 2010); rychlost zakulacování podle Huta (1981, rov. 10) pro vnořený
sekulární backend; a tři příspěvky k apsidálnímu stáčení zákrytového páru A pozorované čtyřhvězdy
CzeV343 (slapový, rotační a relativistický) spočtené Pejchou et al. (2022, tabulka 3).

**Tabulka 15.** Apsidální stáčení a slapová disipace. Sloupce jsou definovány v oddílu [Jak číst tuto zprávu](#how-to-read).

| Test | Třída | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|---|
| `apsidal-sterne-static` | R | 318,227 roku | 318,201 roku | −0,01 % | 0,5 % | PASS |
| `apsidal-sterne-rotating` | R | 291,213 roku | 291,185 roku | −0,01 % | 0,5 % | PASS |
| `tides-leconte-a` | R | 1 | 0,999285 | −0,07 % | 1 % | PASS |
| `tides-leconte-e` | R | 1 | 0,999356 | −0,06 % | 1 % | PASS |
| `tides-leconte-spin` | R | 1 | 1,00113 | +0,11 % | 1 % | PASS |
| `tides-leconte-e-eccentric` | R | 1 | 0,996923 | −0,31 % | 1 % | PASS |
| `secular-nested-tides-hut81` | own | 1 | 1 | −0,00 % | 10⁻⁴ % | PASS |
| `czev343-apsidal-tidal` | own | 4,14·10⁻⁴ rad/d | 4,22272·10⁻⁴ rad/d | +2,00 % | 3 % | PASS |
| `czev343-apsidal-rotational` | own | 6,4·10⁻⁵ rad/d | 6,32533·10⁻⁵ rad/d | −1,17 % | 5 % | PASS |
| `czev343-apsidal-gr` | own | 1,47·10⁻⁵ rad/d | 1,47874·10⁻⁵ rad/d | +0,59 % | 2 % | PASS |
| `czev343-apsidal-period` | R | 33,4 roku | 34,491 roku | +3,27 % | 6 % | PASS |

**Výsledek (tabulka 15).** Osamocená excentrická dvojhvězda rozměrných hvězd ($R/a\approx0{,}14$,
$k_2$ = 0,015, $e$ = 0,3), nejprve bez rotace a pak se synchronní rotací, se v N-částicovém běhu
dlouhém 60 let stáčí s apsidální periodou podle Sterneova vzorce s přesností 10⁻⁴. Se zahrnutým
časovým zpožděním souhlasí N-částicové rychlosti změny $a$, $e$ a rotace s teorií středovanou přes
oběh na 0,06–0,31 %, a to jak při $e$ = 0,3, tak při $e$ = 0,6 s hvězdou rotující 1,2násobkem
středního pohybu. Rychlosti se berou z rozdílu běhů s disipací a bez ní, čímž se odstraní
krátkoperiodické členy konzervativního slapu. Tyto testy také ustalují konvenci: parametr `tau`
v REBOUNDx je dvojnásobkem časového zpoždění SWINGERsimu; bez faktoru ½ ve SWINGERsimu by byly
rychlosti dvakrát větší. Vnořený sekulární backend reprodukuje Hutovu rychlost zakulacování přesně
(jde o tutéž rovnici). Pro CzeV343 souhlasí tři příspěvky s Pejchou et al. (2022) na 0,6–2 %
a úplná N-částicová čtyřhvězda se slapy a rotací dává apsidální periodu páru A 34,5 roku proti
pozorovaným 33,4 ± 0,2 roku (Pejcha et al. 2022, oddíl 4.3) (+3,3 %; tolerance 6 % odpovídá
nejistotě parametrů hvězd). CzeV343 je podrobně studována v části IV.

## Vývojové dráhy hvězd malých a středních hmotností {#tracks}

**Definice.** SWINGERsim sleduje vývoj hvězd pomocí tabelovaných vývojových drah. **Dráhy SSE** jsou
analytické vzorce Hurleyho, Polse & Touta (2000), jednou vyčíslené publikovanými fortranovými
rutinami SSE na mřížce 320 hmotností × 5 metalicit a interpolované SWINGERsimem podél
„ekvivalentních vývojových bodů“ (EEP, podílů každé vývojové fáze); končí zapálením helia. **Dráhy
MIST** (Choi et al. 2016; Dotter 2016) jsou modely MESA převedené z publikovaných archivů EEP, od
0,1 do 300 M☉, s rotací i bez ní, až do konce asymptotické větve obrů (AGB). Veličiny obálky
(hmotnost a poloměr konvektivní obálky) se řídí BSE (Hurley et al. 2002).

**Reference.** Poloměry ZAMS podle Touta et al. (1996); současné Slunce; helioseismologicky určená
spodní hranice sluneční konvektivní zóny, 0,713 $R_\odot$ (Christensen-Dalsgaard et al. 1991);
publikovaná stáří EEP z MIST; MIST proti SSE (dvě rodiny modelů z různých kódů); Reimersův vítr
na větvi rudých obrů; a poloměr plně konvektivního červeného trpaslíka CM Dra A, který standardní
modely podhodnocují (Feiden & Chaboyer 2014).

**Tabulka 16.** Vývojové dráhy hvězd malých a středních hmotností (třída own/tab: interpolace
publikovaných drah ve SWINGERsimu). Sloupce jsou definovány v oddílu [Jak číst tuto zprávu](#how-to-read).

| Test | Třída | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|---|
| `sse-zams-radius-1` | own/tab | 0,88825 R☉ | 0,888312 R☉ | +0,01 % | 0,1 % | PASS |
| `sse-zams-radius-2` | own/tab | 1,79434 R☉ | 1,79435 R☉ | +0,00 % | 0,1 % | PASS |
| `sse-sun-radius` | own/tab | 1 R☉ | 0,984843 R☉ | −1,52 % | 5 % | PASS |
| `sse-sun-luminosity` | own/tab | 1 L☉ | 0,95786 L☉ | −4,21 % | 5 % | PASS |
| `sse-sun-convection-zone` | own/tab | 0,713 | 0,733903 | +0,0209 | 0,03 | PASS |
| `mist-grid-tams-age-1` | own/tab | 9919,42 Myr | 9919,42 Myr | +0,00 % | 0,001 % | PASS |
| `mist-interpolation` | own/tab | 1 | 1,01068 | +1,07 % | 2 % | PASS |
| `mist-ms-lifetime-1` | own/tab | 9878 Myr | 11 003,5 Myr | +11,39 % | 15 % | PASS |
| `mist-ms-lifetime-2` | own/tab | 1054 Myr | 1164,17 Myr | +10,45 % | 15 % | PASS |
| `mist-tams-radius-2` | own/tab | 3,913 R☉ | 3,56583 R☉ | −8,87 % | 15 % | PASS |
| `sse-wind-rgb-1msun` | own/tab | 0,04158 M☉ | 0,0421789 M☉ | +1,44 % | 10 % | PASS |
| `cmdra-a-radius-excess-solar` | own/tab | 0,035 | 0,0335989 | −0,0014 | 0,015 | PASS |
| `cmdra-a-radius-excess-lowz` | own/tab | 0,06 | 0,0641756 | +0,00418 | 0,015 | PASS |

**Výsledek (tabulka 16).** Mřížka SSE reprodukuje poloměry ZAMS podle Touta et al. (1996)
s přesností 10⁻⁴, i pro 2,5 M☉, což leží mezi hmotnostmi mřížky. Hvězda o hmotnosti 1 M☉
a $Z$ = 0,02 má ve stáří 4,57 Gyr $R$ = 0,985 R☉ a $L$ = 0,958 L☉, v rámci přesnosti 5 %, kterou
pro své vzorce uvádějí Hurley et al. (2000); spodní hranice její konvektivní obálky leží na
0,734 $R$, o 0,021 výše než helioseismologická hodnota. Mřížka MIST reprodukuje publikovaná stáří
EEP přesně a vývojová dráha interpolovaná ze sousedů ve dvojnásobném kroku mřížky (1,14 M☉
z 1,12 a 1,16 M☉) má na konci hlavní posloupnosti poloměr v rozmezí 1,1 % od tabelovaného.
Doby života na hlavní posloupnosti pro 1 a 2 M☉ jsou v MIST o 11 % a 10 % delší než v SSE a poloměr
hvězdy 2 M☉ na konci hlavní posloupnosti je o 9 % menší. To jsou známé rozdíly mezi oběma rodinami
modelů (konvektivní přestřelování, škála slunečního složení $Z$ = 0,02 proti 0,0142); určují
systematickou nejistotu stáří a počátku přetečení Rocheova laloku ve SWINGERsimu a toleranci
15 %. Reimersův vítr SSE s $\eta$ = 0,1 odebere hvězdě 1 M☉ do zapálení helia 0,0422 M☉, proti
0,0416 M☉ v modelu MIST se stejným zákonem. CM Dra A je o 3,4 % větší než poloměr podle MIST pro její
hmotnost při sluneční metalicitě a o 6,4 % při [Fe/H] = −0,3, jak pro standardní modely zjišťují
Feiden & Chaboyer (2014) (3–4 % a 5–7 %): SWINGERsim reprodukuje známý problém nafouknutých
poloměrů červených trpaslíků, místo aby jej zakrýval.

## Hmotné hvězdy a obři do konce AGB {#massive-giants}

**Reference.** Pro hmotné hvězdy: interpolace mřížky MIST s vynecháním jednoho bodu (vývojová dráha
spočtená ze svých sousedů a srovnaná s tabelovanou); předpis pro vítr podle Vinka et al. (2001,
rov. 24); doby hoření vodíku v ženevských modelech (Ekström et al. 2012), nezávislém kódu s rotací
i bez ní. Pro obry: vztah mezi počáteční a konečnou hmotností bílých trpaslíků podle Cummingse
et al. (2018), změřený pro 80 bílých trpaslíků ve 13 hvězdokupách se stářími předchůdců podle MIST
(rozptyl 0,06 M☉). Hmotnost heliovo-uhlíkového jádra na konci vývojové dráhy MIST je téměř rovna
hmotnosti budoucího bílého trpaslíka.

**Tabulka 17.** Hmotné hvězdy a konec AGB. Sloupce jsou definovány v oddílu [Jak číst tuto zprávu](#how-to-read).

| Test | Třída | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|---|
| `mist-interpolation-massive-r` | own/tab | 1 | 0,999431 | −0,06 % | 1 % | PASS |
| `mist-interpolation-massive-m` | own/tab | 1 | 1,00107 | +0,11 % | 3 % | PASS |
| `mist-vink2001` | own/tab | 1 | 1,05511 | +5,51 % | 10 % | PASS |
| `mist-geneva-tauh-15-0` | own/tab | 11,015 Myr | 12,3538 Myr | +12,15 % | 15 % | PASS |
| `mist-geneva-tauh-15-4` | own/tab | 13,447 Myr | 12,6825 Myr | −5,69 % | 15 % | PASS |
| `mist-geneva-tauh-60-0` | own/tab | 3,53 Myr | 3,71791 Myr | +5,32 % | 15 % | PASS |
| `mist-geneva-tauh-60-4` | own/tab | 4,465 Myr | 3,94657 Myr | −11,61 % | 15 % | PASS |
| `mist-geneva-tauh-120-0` | own/tab | 2,672 Myr | 2,78283 Myr | +4,15 % | 15 % | PASS |
| `mist-geneva-tauh-120-4` | own/tab | 3,178 Myr | 3,33128 Myr | +4,82 % | 15 % | PASS |
| `agb-final-mass-1p5` | own/tab | 0,609 M☉ | 0,565426 M☉ | −0,0436 M☉ | 0,1 M☉ | PASS |
| `agb-final-mass-2p5` | own/tab | 0,689 M☉ | 0,610396 M☉ | −0,0786 M☉ | 0,1 M☉ | PASS |
| `agb-final-mass-3p0` | own/tab | 0,745 M☉ | 0,653775 M☉ | −0,0912 M☉ | 0,1 M☉ | PASS |
| `agb-final-mass-4p0` | own/tab | 0,899 M☉ | 0,835342 M☉ | −0,0637 M☉ | 0,1 M☉ | PASS |

**Výsledek (tabulka 17).** Poloměr hvězdy 60 M☉ uprostřed její hlavní posloupnosti, interpolovaný
z 55 a 65 M☉, leží v rozmezí 0,06 % od tabelovaného a současná hmotnost hvězdy 100 M☉ na konci
hlavní posloupnosti v rozmezí 0,11 %. Meze interpolace leží v pozdějším vývoji: od pozdní hlavní
posloupnosti dále mřížka nerozlišuje nafouknuté obálky hvězd 40–70 M☉, s chybami poloměru až
20–40 % a chybami hmotnosti až 7 % (například −24 % v poloměru pro 60 M☉ na konci hlavní
posloupnosti bez rotace). Rychlost ztráty hmoty hvězdy 40 M☉ v MIST se řídí Vinkem et al. (2001)
s posunem +5,5 %, který odpovídá metalicitě předané předpisu. Doby hoření vodíku pro 15, 60
a 120 M☉ souhlasí s ženevskými modely na 4–12 %, s rotací i bez ní. Hmotnosti na konci hlavní
posloupnosti se liší více (60 M☉ bez rotace: 46,6 proti 36,3 M☉; 120 M☉ s rotací: 25,3 proti
34,6 M☉): právě tento rozptyl mezi rodinami modelů, nikoli interpolace, je nejistotou historie
hmotnosti velmi hmotných hvězd.

Hmotnosti bílých trpaslíků z hvězd MIST o hmotnostech 1,5–4 M☉ leží o 0,04–0,09 M☉ pod vztahem
Cummingse et al. (2018), což je posun, který pro tyto modely uvádějí sami autoři; tolerance
0,1 M☉ je tento posun plus rozptyl vztahu. Test srovnává jádro, protože vývojové dráhy MIST končí
s obálkou až 0,15–0,2 M☉, kterou by odstranila fáze po AGB. Odvozené veličiny ověřují automatické
testy: konvektivní obálka Slunce má 0,034 M☉ a 0,29 $R$, s dobou obratu 21 d.

## Inerciální vlny a slapové zakulacování drah {#circularisation}

**Definice.** V konvektivních obálkách je slapové proudění tlumeno turbulentní viskozitou
(**rovnovážný slap**). V rotujících hvězdách slap budí také **inerciální vlny**, vlny s vratnou
silou Coriolisovy síly, které mohou disipovat mnohem více energie (Ogilvie 2013; Barker 2020,
2022). Účinnost se vyjadřuje modifikovaným slapovým faktorem kvality $Q'$ (větší $Q'$ = slabší
disipace). **Perioda zakulacení** populace dvojhvězd stejného stáří je perioda, pod níž jsou dráhy
kruhové; pozorovatelé ji definují z konečných drah (Meibom & Mathieu 2005).

**Reference.** Uzavřený vzorec pro disipaci inerciálních vln ve dvouvrstvé hvězdě (Ogilvie 2013,
v podobě Barkera 2020, rov. 39); $Q' \approx 10^7$ pro hvězdy hlavní posloupnosti o hmotnostech 0,3
a 1 M☉ rotující s periodou 10 d (Barker 2022, z profilů MESA); zakulacení ve fázi před hlavní
posloupností, 7,2–8,5 d (Zahn & Bouchet 1989); perioda zakulacení dvojhvězd před hlavní
posloupností, 7,1 ± 1,2 d, a perioda otevřené hvězdokupy M67 ve stáří 4 Gyr, 12,1 +1,0/−1,5 d (obojí
Meibom & Mathieu 2005, tabulka 3; citováno Mazehem 2008).

**Tabulka 18.** Inerciální vlny a slapové zakulacování drah. Sloupce jsou definovány v oddílu [Jak číst tuto zprávu](#how-to-read).

| Test | Třída | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|---|
| `iw-two-layer` | own | 1 | 1,00031 | +0,03 % | 0,5 % | PASS |
| `iw-q-solar` | own | 7 | 7,0122 | +0,0122 | 0,3 | PASS |
| `iw-q-m-dwarf` | own | 7 | 6,94697 | −0,053 | 0,3 | PASS |
| `tides-circularization-pms` | own | 7,85 d | 7,9307 d | +1,03 % | 10 % | PASS |
| `tides-circularization-pms-population` | own | 7,1 d | 6,5831 d | −0,517 d | 1,2 d | PASS |
| `tides-circularization-m67` | own | 12,1 d | 3,7847 d | −68,72 % | 30 % | **EXPLAINED** |
| `tides-circularization-m67-iw` | own | 12,1 d | 7,4572 d | −38,37 % | 30 % | **EXPLAINED** |
| `tides-circularization-m67-final` | own | 12,1 d | 7,5337 d | −4,57 d | 1,2 d | **EXPLAINED** |

**Výsledek (tabulka 18).** Řešič disipace inerciálních vln s úplným profilem hvězdy ve SWINGERsimu
reprodukuje uzavřený dvouvrstvý vzorec na 0,03 % pro vyhlazenou dvouvrstvou hvězdu; tím se ověřuje
řešič i jeho normování. Pro hvězdy hlavní posloupnosti 0,3 a 1 M☉ rotující s periodou 10 d dává
$\log Q'$ = 6,95 a 7,01, stejně jako Barker (2022). Dvojhvězdy dvou stejných hvězd 1 M☉ sledované
po vývojových drahách MIST před hlavní posloupností od 1 Myr s rovnovážným slapem a inerciálními
vlnami jsou při dosažení hlavní posloupnosti kruhové až do 7,9 d (Zahn & Bouchet: 7,2–8,5 d)
a populace takových dvojhvězd ve stáří 3 Myr má v pozorovatelské definici (konečná perioda,
$e$ < 0,01) periodu zakulacení 6,6 d, proti pozorovaným 7,1 ± 1,2 d u dvojhvězd před hlavní
posloupností.

Tři testy M67 jsou vysvětleny (EXPLAINED) jako otevřená otázka slapového modelu SWINGERsimu. Se
samotným rovnovážným slapem podle Hurleyho et al. (2002) jsou dvojhvězdy dvou stejných hvězd
slunečního typu, začínající na hlavní posloupnosti s $e$ = 0,3, po 4 Gyr kruhové jen pod 3,8 d;
s inerciálními vlnami pod 7,5 d; a s celou historií od fáze před hlavní posloupností 7,5 d
v konečné periodě. Fáze před hlavní posloupností je reprodukována, hlavní posloupnost však přidá
jen asi 1 d, kdežto M67 vyžaduje růst na 12 d. Barker (2022) zjišťuje, že inerciální vlny mohou
vysvětlit periody zakulacení dvojhvězd slunečního typu s periodou zakulacení, která s věkem na
hlavní posloupnosti roste tak, jak se pozoruje. Odhad disipace inerciálních vln ve SWINGERsimu je
rovněž středován přes frekvence, používá však dvouvrstvé přiblížení (Ogilvie 2013) stavby hvězd
SSE a MIST a sleduje rotace hvězd, kdežto Barker počítá disipaci z úplných profilů MESA
a předpokládá synchronizované rotace. Rozdíl tedy spočívá ve stavbě hvězd a v historii rotace;
která z obou příčin jej způsobuje, zde rozhodnuto není. (Rezonanční zamykání, které SWINGERsim
nemodeluje, působí hlavně ve fázi před hlavní posloupností, se zanedbatelným příspěvkem na hlavní
posloupnosti; Zanazzi & Wu 2021.) Tolerance se nerozšiřuje; rozdíl je uveden jako mez slapového
modelu SWINGERsimu.

## Ztráta hmoty hvězdným větrem {#mass-loss}

**Reference.** Dvojhvězda, jejíž hvězda ztrácí hmotu pomalu (ve srovnání s oběžnou periodou),
zachovává adiabatický invariant $aM$ (Jeans 1924; Hadjidemetriou 1963). Během jednoho oběhu se
excentricita řídí vztahem $e = e_0 + \Psi(1-e_0^2)^{3/2}\sin f/(1+e_0\cos f)$, získaným integrací
rov. 4 Verase et al. (2011) přes pravou anomálii $f$, kde $\Psi$ je hmotnost ztracená za jeden
oběh v jednotkách celkové hmotnosti; počátečně kruhová dráha dosáhne $e = 2\Psi$. Okamžitá ztráta
podílu $\delta$ hmotnosti dává $a'/a = (1-\delta)/(1-2\delta)$ a $e' = \delta/(1-\delta)$ (Blaauw
1961; Hills 1983). V trojhvězdě, jejíž vnitřní primární složka ztrácí hmotu, mohou rozpínající se
dráhy překročit mez stability („triple evolution dynamical instability“ podle Peretse & Kratterové
2012).

V N-částicovém backendu se hmotnosti mění operátorem `modify_mass` knihovny REBOUNDx podle
časového plánu ztráty hmoty SWINGERsimu (třída R); v sekulárním a vývojovém backendu se dráhy
přeškálují vlastními rovnicemi SWINGERsimu (třída own).

**Tabulka 19.** Ztráta hmoty. Sloupce jsou definovány v oddílu [Jak číst tuto zprávu](#how-to-read).

| Test | Třída | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|---|
| `ml-adiabatic-invariant` | R | 0 | 3,57·10⁻⁹ | +3,57·10⁻⁹ | 10⁻⁶ | PASS |
| `ml-osculating-e-03` | R | 0 | 0,0126471 | +0,0126 | 0,05 | PASS |
| `ml-osculating-e-07` | R | 0 | 0,0123538 | +0,0124 | 0,05 | PASS |
| `ml-circular-emax` | R | 2 | 2,00265 | +0,13 % | 2 % | PASS |
| `ml-impulsive-a` | R | 1,21429 | 1,21491 | +0,05 % | 1 % | PASS |
| `ml-impulsive-e` | R | 0,176471 | 0,176892 | +0,24 % | 1 % | PASS |
| `ml-zlk-nbody-secular` | own | 1 | 1,0021 | +0,21 % | 1 % | PASS |
| `ml-tedi-onset` | R | 1 | 1,06475 | +6,47 % | 15 % | PASS |
| `wind-transfer-jeans` | own | 0 | 1,69·10⁻¹¹ | +1,69·10⁻¹¹ | 10⁻⁸ | PASS |

**Výsledek (tabulka 19).** Dvojhvězda, jejíž hvězda pomalu ztratí 0,3 M☉, zachovává $aM$
v N-částicovém backendu na 4·10⁻⁹; sekulární backendy jej z konstrukce zachovávají přesně. Během
dvou oběhů se excentricita řídí výše uvedeným vzorcem na 1,3 % z $\Psi$ pro $e_0$ = 0,3 a 0,7
a kruhová dráha dosáhne $e/\Psi$ = 2,003. Veras et al. tisknou jmenovatel své rov. 17 jako
$1-e_0\cos f$; tento tvar z jejich rov. 4 neplyne a s N-částicovými běhy nesouhlasí asi o $\Psi$,
takže jej považujeme za tiskovou chybu. Ztráta 15 % hmotnosti během 10⁻³ oběhu dává okamžité
výsledky na 0,05 % a 0,24 %; zbytek odpovídá konečné době trvání ztráty. Zatímco rušící těleso ZLK
ztrácí polovinu své hmotnosti, N-částicový a sekulární backend se shodují na okamžiku čtvrtého
maxima ZLK na 0,2 %. Pro hmotnosti Peretse & Kratterové (2012) (vnitřní dvojhvězda 5 + 2 M☉, třetí
hvězda 2,25 M☉) a vnitřní primární složku, která ztratí 3 M☉ za 8000 let, začne vnitřní
excentricita růst o 6,5 % později, než adiabaticky se rozpínající dráhy překročí mez MA01; soustava
se nerozpadla ani po dalších 2000 vnějších obězích, v souladu s dlouhým přežíváním na hraně, které
zjistili Toonen et al. (2022). S větrem a přenosem hmoty současně zachovává dvojhvězda
25 + 18 M☉ veličinu $a(M_1+M_2)$ na 2·10⁻¹¹, než přenos začne.

## Přenos hmoty a omlazení {#mass-transfer}

**Definice.** Hvězda, která vyplní svůj Rocheův lalok, přenáší hmotu na svého průvodce
(**přetečení Rocheova laloku**). Jak dráha odpoví, závisí na tom, kolik přenesené hmoty si příjemce
ponechá a kolik momentu hybnosti odnese ztracená hmota. Modul přenosu hmoty je vlastním kódem
SWINGERsimu: sleduje dárce po jeho vývojové dráze, počítá rychlost přenosu, odezvu dráhy
a Rocheova laloku a stabilitu přenosu. Příjemce na hlavní posloupnosti, který získává hmotu, je
**omlazen**: čerstvý vodík způsobí, že se jeví mladší (Hurley et al. 2000).

**Reference.** Uzavřené tvary odezvy dráhy při konstantních podílech hmoty ztracené z dárce
a z příjemce (Soberman, Phinney & van den Heuvel 1997, rov. 20–28) a pro konzervativní přenos
a izotropní reemisi (Postnov & Yungelson 2014, rov. 37 a 43); odezva Rocheova laloku
$\zeta_L = d\ln R_L/d\ln M_d$, přiblížená výrazem $2{,}13q - 1{,}67$ (Tout et al. 1997, citováno
Hurleym et al. 2002); **příklad Algolu z BSE** podle Hurleyho et al. (2002, oddíl 3.1) –
2,9 + 0,9 M☉, $Z$ = 0,02, se začátkem při $P$ = 8 d a $e$ = 0,7, kterou slapy zakulatí na
$P$ = 2,7 d v okamžiku, kdy primární složka vyplní Rocheův lalok; testy začínají z tohoto stavu –
s BSE spuštěným nezávisle v COSMIC 4.2.1 (Breivik et al. 2020); a oprava stáří omlazeného příjemce
(Hurley et al. 2000).

**Tabulka 20.** Přenos hmoty a omlazení (vše třída own). Sloupce jsou definovány v oddílu [Jak číst tuto zprávu](#how-to-read).

| Test | Třída | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|---|
| `mt-soberman-general` | own | 1,29573 $a/a_0$ | 1,29573 $a/a_0$ | −0,00 % | 10⁻⁵ % | PASS |
| `mt-conservative` | own | 2,25 $a/a_0$ | 2,25 $a/a_0$ | −0,00 % | 10⁻⁵ % | PASS |
| `mt-isotropic-reemission` | own | 2,03003 $a/a_0$ | 2,03003 $a/a_0$ | −0,00 % | 10⁻⁵ % | PASS |
| `mt-zeta-lobe-q2` | own | 2,59 | 2,60707 | +0,0171 | 0,05 | PASS |
| `mt-bse-algol-donor` | own | 0,42 M☉ | 0,404336 M☉ | −0,0157 M☉ | 0,05 M☉ | PASS |
| `mt-cosmic-algol-donor` | own | 0,401 M☉ | 0,404336 M☉ | +0,00334 M☉ | 0,03 M☉ | PASS |
| `mt-bse-algol-accretor` | own | 2,85 M☉ | 2,51105 M☉ | −0,339 M☉ | 0,15 M☉ | **EXPLAINED** |
| `mt-bse-algol-period` | own | 30 d | 21,9301 d | −26,90 % | 15 % | **EXPLAINED** |
| `mt-bse-algol-giant-phase` | own | 30 d | 30,3561 d | +1,19 % | 3 % | PASS |
| `rejuvenation-hurley` | own | 1 | 0,999174 | −0,08 % | 0,2 % | PASS |

**Výsledek (tabulka 20).** Numericky integrovaná odezva dráhy souhlasí s uzavřenými tvary
na 3·10⁻¹³ (obecné podíly: 30 % přenesené hmoty ztraceno z dárce, 40 % reemitováno z příjemce,
30 % akreováno) a 2·10⁻¹² (konzervativní přenos, izotropní reemise). $\zeta_L$ při $q$ = 2 je
2,607 proti 2,59 z aproximace; aproximace je přesná na 0,03 pro $q$ = 0,5–3. V příkladu Algolu
z BSE ukončí SWINGERsim, který hvězdy přemodelovává stejně jako BSE, první přenos s dárcem
0,404 M☉ (článek 0,42 M☉, BSE v COSMIC 0,401 M☉; obojí PASS), s příjemcem 2,51 M☉ a s $P$ = 21,9 d.
Článek uvádí 2,85 M☉ a 30 d. Tyto dva testy jsou vysvětleny (EXPLAINED), protože publikovaná
historie není vnitřně konzistentní: z 2,7 d na 17 d na konci Hertzsprungovy mezery by musel
orbitální moment hybnosti *vzrůst* nejméně o 8 % (horní mez bez jakékoli ztráty dává 13,5 d),
kdežto následující fáze obra, ze 17 d na 30 d, je přesně konzervativní. SWINGERsim, spuštěný
z publikovaného stavu na konci Hertzsprungovy mezery (0,53 + 2,74 M☉, 17 d), tuto fázi obra
reprodukuje (30,4 d proti 30 d). Samotný kód BSE, spuštěný v COSMIC 4.2.1 se stejným vstupem
a s akreční mezí BSE, končí při 3,2 d s příjemcem 1,07 M☉: kolik hmoty si příjemce na tepelné
časové škále ponechá, závisí na jeho tepelné časové škále (COSMIC v odpovídající rovnici používá
3,13·10⁷ místo 10⁷ let) a na implementaci BSE. Jde o nejméně jistou veličinu vývoje Algolů
a SWINGERsim ji nechává nastavit uživateli (`accretion: fractions`). Stáří příjemce 2 M☉ poté,
co získá 1 % své hmotnosti, se řídí vzorcem omlazení podle Hurleyho et al. (2000) na 0,08 %.

## Společná obálka, magnetické brzdění, Applegateův mechanismus a slapový záchyt {#ce-magnetism}

**Definice.** Stane-li se přenos hmoty dynamicky nestabilním, obálka dárce pohltí obě hvězdy
(**společná obálka**); v energetickém formalismu orbitální energie uvolněná při spirálování
dovnitř uvolní obálku s účinností $\alpha_\mathrm{CE}$. **Magnetické brzdění** je ztráta rotačního
momentu hybnosti chladné hvězdy jejím magnetizovaným větrem. **Applegateův mechanismus**
vysvětluje cyklické změny periody těsných dvojhvězd magnetickým cyklem, který mění kvadrupólový
moment jedné z hvězd; otázkou je, zda má hvězda pro to dostatek energie. **Slapový záchyt**: dvě
hvězdy, které proletí dostatečně blízko u sebe, ztratí slapy dost orbitální energie, aby se staly
vázanými.

**Reference.** Dráha s periodou 35 h po společné obálce předchůdce kataklyzmické proměnné podle
Touta et al. (1997) – 7 + 1 M☉, $P$ = 1200 d, jejíž primární složka vyplní Rocheův lalok na AGB při
5,36 M☉ s jádrem 1,014 M☉ (nejde o příklad kataklyzmické proměnné podle Hurleyho et al. 2002,
oddíl 3.2, který začíná z 6,0 + 1,3 M☉ na kruhové dráze s periodou 630 d) – a dráha s periodou
0,04 d po druhé společné obálce v příkladu Algolu z BSE; impulzní odezva vnější dráhy na náhlou
ztrátu obálky (Hadjidemetriou 1963; Veras et al. 2011); moment síly větru podle Matta et al. (2015)
a Skumanichův zákon $\Omega\propto t^{-1/2}$ (Skumanich 1972); energetické poměry v modelu tenké
slupky pro HU Aqr, NN Ser a V471 Tau (Völschow et al. 2016, tabulka 4, sloupec modelu tenké slupky
podle Tiana et al. 2009); a poloměry záchytu podle Raye et al. (1987), citované Portegiesem Zwartem
& Meinenem (1993), a podle McMillana et al. (1987).

**Tabulka 21.** Společná obálka, magnetismus a slapový záchyt. Sloupce jsou definovány v oddílu [Jak číst tuto zprávu](#how-to-read).

| Test | Třída | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|---|
| `ce-tout97-period` | own | 35 h | 34,84 h | −0,46 % | 10 % | PASS |
| `ce-bse-algol-second-period` | own | 0,04 d | 0,0365638 d | −8,59 % | 15 % | PASS |
| `ce-impulsive-outer-a` | R | 1 | 1,01274 | +1,27 % | 3 % | PASS |
| `ce-impulsive-outer-e` | R | 0 | 0,00442155 | +0,00442 | 0,02 | PASS |
| `matt15-sun-period` | own | 27,97 d | 27,7635 d | −0,74 % | 2 % | PASS |
| `matt15-unsaturated-powerlaw` | own | 0,5 | 0,491949 | −0,00805 | 0,02 | PASS |
| `applegate-huaqr` | own | 0,1 | 0,105659 | +5,66 % | 10 % | PASS |
| `applegate-nnser` | own | 3,3 | 3,30338 | +0,10 % | 10 % | PASS |
| `applegate-v471tau` | own | 0,014 | 0,0143688 | +2,63 % | 10 % | PASS |
| `tidal-capture-radius-n15-mcmillan` | own | 3,4 | 3,33454 | −1,93 % | 5 % | PASS |
| `tidal-capture-radius-n3` | own | 2,72 | 2,36599 | −13,01 % | 5 % | **EXPLAINED** |
| `tidal-capture-radius-n15` | own | 3,84 | 3,33454 | −13,16 % | 5 % | **EXPLAINED** |

**Výsledek (tabulka 21).** Energetický formalismus reprodukuje dráhu s periodou 35 h (34,8 h; dárce
vyplní svůj lalok na dráze rozšířené jeho větrem) a dráhu s periodou 0,04 d (0,037 d
s $\alpha_\mathrm{CE}$ = 3). Dvojhvězda 1,8 + 0,5 M☉, jejíž obr projde společnou obálkou uvnitř
kruhové vnější dráhy s periodou 10 let (N-částicově), rozšíří vnější dráhu o 1,3 % více, než
předpovídá impulzní přiblížení, a učiní ji excentrickou s odchylkou do 0,004 od předpovědi; zbytek
odpovídá konečné době trvání odvržení. Model Slunce začínající při $P$ = 1 d rotuje ve slunečním
stáří s periodou 27,76 d (sluneční hodnota $2\pi/\Omega_\odot$ je 27,97 d) a rotace se řídí
Skumanichovým zákonem (poměr 0,492 proti 0,5 mezi 1,14 a 4,55 Gyr). Applegateovy energetické poměry
souhlasí s Völschowem et al. na 0,1–5,7 %.

Dva testy slapového záchytu jsou vysvětleny (EXPLAINED) jako rozdíl v definici. SWINGERsim používá
aproximace Portegiese Zwarta & Meinena (1993) (implementované v NBODY6); pro dvě hvězdy 0,6 M☉
setkávající se rychlostí 10 km/s dávají poloměry záchytu 2,37 (polytropa $n$ = 3) a 3,33
(n = 1,5) poloměru hvězdy, o 13 % méně než hodnoty Raye et al. (1987) citované v tomto článku.
Při citovaných poloměrech je však slapová energie 5krát a 13krát menší než orbitální energie, kterou
je třeba odebrat, takže citované hodnoty nemohou patřit ke stejné rychlosti a stavbě hvězdy; patří
k jiným konvencím Raye et al. Tytéž aproximace souhlasí s McMillanem et al. (1987), asi 3,4 $R$ pro
$n$ = 1,5, na 1,9 % (PASS).

## Kontaktní dvojhvězdy a dynamický počátek přenosu hmoty {#contact}

**Definice.** **Kontaktní dvojhvězda** je pár, jehož obě hvězdy přeplňují své Rocheovy laloky
a sdílejí společnou obálku (hvězdy typu W UMa). Ztrácí moment hybnosti (Gazeas & Stępień 2008)
a splyne, když dosáhne Darwinovy nestability nebo ztrácí hmotu vnějším Lagrangeovým bodem L2
(Rasio 1995). K **dynamickému počátku** přenosu hmoty dochází, když cykly ZLK zvýší excentricitu
vnitřní dvojhvězdy tak, že se hvězda v periastru dotkne svého Rocheova laloku; zakulatí-li pak
slapy dráhu se zachováním momentu hybnosti, je vzdálenost na kruhové dráze $a(1-e_c^2)$ a dárce
vyplňuje svůj lalok $(1+e_c)$krát méně než v okamžiku dotyku (Hut 1981).

**Tabulka 22.** Kontaktní dvojhvězdy a dynamický počátek. Sloupce jsou definovány v oddílu [Jak číst tuto zprávu](#how-to-read).

| Test | Třída | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|---|
| `contact-jdot` | own | 2,18893 Gyr | 2,18898 Gyr | +0,00 % | 0,1 % | PASS |
| `dynamic-onset-circularisation` | own | 1 | 1 | −0,00 % | 0,01 % | PASS |
| `dynamic-onset-nbody-vs-secular` | own | 1 | 1,0173 | +1,73 % | 3 % | PASS |

**Výsledek (tabulka 22).** Doba ztráty momentu hybnosti kontaktní dvojhvězdy 1,0 + 0,5 M☉
s $P$ = 0,35 d, spočtená v jednotkách SWINGERsimu, se rovná hodnotě spočtené v jednotkách cgs podle
Gazease & Stępieńa (2008) na 2·10⁻⁵. Vyplnění laloku dárce po zakulacení při $e_c$ = 0,89, vynásobené
$(1+e_c)$, je rovno 1 na 2·10⁻⁹. Okamžik dotyku vyvolaného cykly ZLK (asi 61 let) se mezi
N-částicovým a sekulárním backendem liší o 1,7 %; rozdíl tvoří krátkoperiodické členy vnější dráhy
($a_\mathrm{out}/a_\mathrm{in}$ = 15), které dvojitě středované rovnice vynechávají. Oba backendy
kladou počátek do excentricity, při níž $a(1-e_c) = R/r_L$, kde $r_L$ je Eggletonův poloměr laloku
v jednotkách $a$. Úplný vývoj dvojhvězdy do kontaktu a její splynutí ukazuje
[kapitola 8](#life-cycles).
