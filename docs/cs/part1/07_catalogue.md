# Katalog publikovaných vícenásobných hvězd {#catalogue}

Srovnávací testy kapitol 2–4 ověřují vždy jeden jev. Skutečná vícenásobná hvězda je ověřuje
společně, s hvězdami ze všech částí Hertzsprungova–Russellova diagramu. Distribuce SWINGERsimu
obsahuje 29 konfigurací 23 publikovaných soustav (tabulka 34): Wolfovy–Rayetovy hvězdy a hvězdy
třídy O, pulzující hvězdy typu β Cep, δ Sct a cefeidy, obry, mladé hvězdy typu T Tauri, bílé
trpaslíky, milisekundový pulzar a hvězdy skryté prachovými disky. Každá konfigurace uvádí své
zdroje (článek a tabulku), každý převod (například relativní dráha používá $\omega$ = publikované
$\omega$ primární složky + 180°; střední anomálie plyne z okamžiku průchodu periastrem) a každý
předpoklad, označený ASSUMED nebo UNKNOWN; kde je veličina neznámá, losuje ji ansámbl z apriorního
rozdělení. Níže je srovnána každá veličina, kterou SWINGERsim počítá a která byla také publikována.
Rozdíl označujeme *souhlasí*, leží-li v rámci publikované nejistoty nebo asi 10 %, *vysvětleno*,
byl-li vysledován ke vstupu, definici nebo publikovaným číslům, a *předpověď*, není-li s čím
srovnat.

**Tabulka 34.** Katalog. Sloupec 2 uvádí hierarchii v zápisu SWINGERsimu (2+1 = dvojhvězda se
třetí hvězdou; 1+2 = osamocená hvězda, kolem níž obíhá dvojhvězda), sloupec 3 hvězdy, sloupec 4
hlavní zdroj konfigurace.

| Soustava | Hierarchie | Hvězdy | Zdroj |
|---|---|---|---|
| Algol | 2+1 | B8 V + K0 IV (vyplňuje svůj lalok) + Am | Baron et al. 2012 |
| BU CMi | 2+2 | čtyři hvězdy třídy A | Pribulla et al. 2023 |
| VW LMi | 2+2 | kontaktní pár typu W UMa + pár hvězd G | Pribulla et al. 2020 |
| ξ Tau | 3+1 | B9 V + B9 V, B5 V, F5 V | Nemravová et al. 2016 |
| TIC 120362137 | 3+1 | F + F + F + G | Borkovits et al. 2026 |
| TIC 454140642 | 2+2 | čtyři hvězdy třídy F | Kostov et al. 2021 |
| TIC 470710327 | 2+1 | B + B + O9.5 V | Eisner et al. 2022 |
| δ Ori A | 2+1 | O9.5 II + B2 V + B0 IV | Oplištilová et al. 2023 |
| σ Ori | 2+1 | O9.5 V + B0.5 V + B | Schaefer et al. 2016 |
| λ Sco | 2+1 | β Cep + hvězda před hlavní posloupností + B2 V | Tango et al. 2006 |
| Polárka | 2+1 | cefeida pulzující v prvním svrchním tónu + F6 V + F3 V | Evans et al. 2024 |
| HD 181068 | 1+2 | rudý obr + G + K | Borkovits et al. 2013 |
| PSR J0337+1715 | 2+1 | milisekundový pulzar + dva bílí trpaslíci | Ransom et al. 2014 |
| KIC 4150611 | 2+(2+1) | δ Sct/γ Dor + M + M, G + G | Hełminiak et al. 2017 |
| HD 98800 | 2+2 + disk | čtyři hvězdy T Tauri, 10 Myr, polární disk | Zúñiga-Fernández et al. 2026; Kennedy et al. 2019 |
| HD 5980 | 2+2 | WR/LBV + WN4, veleobr třídy O + neviditelná složka | Koenigsberger et al. 2014 |
| HD 150136 | 2+1 | O3 V + O5.5 V + O6.5 V | Mahy et al. 2018 |
| α Cen | 2+1 | G2 V + K1 V + M5.5 V | Akeson et al. 2021; Kervella et al. 2017 |
| GG Tau | 2+(2+1) + prstenec | pět hvězd T Tauri / hnědý trpaslík | Köhler 2011; Keppler et al. 2020 |
| TIC 168789840 | 2+2+2 | tři zákrytové páry F + K/M | Powell et al. 2021 |
| ε Aur | 1+1 + disk | veleobr třídy F + hvězda B v prachovém disku | Stefanik et al. 2010; Hoard et al. 2010 |
| EE Cep | 1+1 + disk | B5 III + tmavý disk (ilustrativní) | Pieńkowski et al. 2020 |
| KH 15D | 2 + prstenec | dvojhvězda T Tauri za prstencem (ilustrativní) | Winn et al. 2006 |

Soustavy s prachovými disky jsou rozebrány v [kapitole 4](#discs); ostatní následují, seskupené
podle fyziky, kterou ověřují. Třídy ve smyslu tabulky 1: běhy precese a časování používají
N-částicový backend (R) nebo sekulární backendy (own), vývojové běhy a časování pulzací používají
vlastní kód SWINGERsimu.

## Apsidální stáčení a precese uzlů {#catalogue-precession}

**Metoda.** Otáčení elips (apsidální) a rovin drah (uzlové) se měří přímo ze zintegrovaných
dráhových elementů – se všemi hvězdami, relativitou, slapy a rotací – a srovnává se
s publikovanými rychlostmi, které autoři odvodili z okamžiků zákrytů, fotodynamických modelů nebo
analytických vzorců. Běhy pokrývají 30–1000 let.

**Tabulka 35.** Apsidální stáčení a precese uzlů v katalogu. Sloupec 4 je hodnota SWINGERsimu,
sloupec 5 publikovaná hodnota, sloupec 6 jejich rozdíl.

| Soustava | Veličina | SWINGERsim | Publikováno | Rozdíl | Stav |
|---|---|---|---|---|---|
| BU CMi | apsidální perioda A | 27,8 roku | 28,7 roku (Pribulla et al. 2023) | −3 % | souhlasí |
| BU CMi | apsidální perioda B | 24,8 roku | 25,1 roku | −1 % | souhlasí |
| ξ Tau | rychlost apsidálního stáčení dráhy s periodou 145 d | 2,14 °/rok | 2,02 ± 0,31 °/rok | +6 % (0,4σ) | souhlasí |
| TIC 120362137 | apsidální perioda dráhy s periodou 51 d | 32,4 roku | 34,79 roku (Borkovits et al. 2026) | −7 % | souhlasí |
| TIC 120362137 | perioda precese uzlů vnějších rovin | 282 let | 265,8 roku | +6 % | souhlasí |
| TIC 454140642 | apsidální periody A, B | 69,7, 90,9 roku | 70, 97 let (Kostov et al. 2021) | −0,4 %, −6 % | souhlasí |
| HD 181068 | perioda precese uzlů | 12,8 roku | 13–14 let (Borkovits et al. 2013) | −5 % | souhlasí |
| VW LMi | apsidální perioda odděleného páru | 95,1 roku | 95,5 roku (vlastní N-částicový model článku) | −0,4 % | souhlasí |
| δ Ori A | rychlost apsidálního stáčení dráhy s periodou 5,73 d | 2,42 °/rok | 1,45 °/rok | +67 % | vysvětleno: vstup |
| TIC 470710327 | perioda precese uzlů zákrytového páru | 14,1 roku | – | – | předpověď |

![**Obrázek 13.** Precese změřená ve čtyřech bězích katalogu; v každém panelu úhel s časem rovnoměrně roste, sklon je rychlostí precese a 360° dělených rychlostí je perioda. *Vlevo nahoře:* délka periastra obou párů BU CMi (A a B). *Vpravo nahoře:* totéž pro oddělený pár VW LMi (modře) proti čárkované přímce periody 78,6 roku, kterou Pribulla et al. (2020) proložili třemi roky dat. *Vlevo dole:* uzly tří drah TIC 120362137 v invariantní rovině (rovině kolmé k celkovému momentu hybnosti) během 1000 let. *Vpravo dole:* sklony obou vnitřních drah ξ Tau k rovině oblohy během 300 let.](../../figures/studies/catalogue_precession.png)

**Výsledek (tabulka 35).** Obrázek 13 ukazuje čtyři z běhů. V BU CMi se 360° dosáhne asi po 28
letech (A) a 25 letech (B), o 3 % a 1 % rychleji, než je publikováno. Ve VW LMi se oddělený pár
otočí za 95 let, pomaleji než čárkovaná přímka 78,6 roku: hodnota 78,6 roku podle Pribully et al.
(2020) je proložením jen tří let radiálních rychlostí, kdežto vlastní analytický a N-částicový
model článku dávají 96,2 a 95,5 roku (a data DDO z let 1998–2008 dávají 91 ± 11 let); SWINGERsim
souhlasí s dynamickým modelem článku na 0,4 %. Uzly TIC 120362137 se otáčejí všechny společně,
zpětně, o 360° asi za 280 let: roviny se stáčejí jako jeden blok. Pro tuto soustavu udává tabulka 1
Borkovitse et al. (2026) pro excentrickou dráhu s periodou 51 d hodnotu 34,79 roku, proti
32,4 roku ze SWINGERsimu (−7 %). Téměř kruhový pár s periodou 3,3 d ($e$ = 0,007) je třetí hvězdou
vychylován na $e$ mezi 0,0001 a 0,016; jeho oskulační $\omega$ skáče a měřená apsidální perioda
není definována (článek pro ni uvádí jen nominální teoretickou hodnotu). Sklon zákrytového
vnitřního páru ξ Tau kolísá o ±0,6° s periodou asi 20 let, takže hloubka jeho zákrytů by se měla
s touto periodou měnit – předpověď.

Jediný velký rozdíl, δ Ori A, pochází ze vstupu. Jejímu apsidálnímu stáčení dominují slapy
a rotace (samotná dvojhvězda dává apsidální periodu 155 let) a slapová část škáluje jako $k_2R^5$
primární složky. S poloměrem podle Oplištilové et al. (2023) a $k_2$ podle Clareta (2004) je
SWINGERsim 1,7krát rychlejší než rychlost podle Pabla et al. (2015); publikovaná hodnota
1,45 °/rok vyžaduje $R_1\approx$ 11,5 R☉ (o 11 % méně) nebo $\log k_2\approx$ −2,85. Publikovaná
rychlost je proto spíše měřením $k_2R^5$ složky Aa1 než testem dynamiky. U TIC 470710327 se při
vzájemném sklonu 17° stáčí pár s periodou 1,1 d za 14 let a jeho zákryty by se měly změnit během
let (předpověď).

## Změny okamžiků zákrytů {#catalogue-timings}

**Metoda.** `swingersim timings` počítá z N-částicové integrace přesné okamžiky konjunkcí včetně
doby šíření světla pro TIC 470710327 (5 let), TIC 454140642 (5 let), BU CMi (30 let) a ξ Tau
(10 let). Výsledné diagramy O−C jsou tím, co by získal pozorovatel.

![**Obrázek 14.** Syntetické diagramy O−C čtyř soustav katalogu; vodorovné osy: čas; svislé osy: o kolik dříve nebo později přichází každý zákryt, než předpovídají přísně periodické hodiny. *Vlevo nahoře* (TIC 470710327, minuty): pár s periodou 1,1 d s dráhou třetí hvězdy s periodou 52 d. *Vpravo nahoře* (ξ Tau, minuty): pár s periodou 7,1 d s dráhou s periodou 145 d. *Vlevo dole* (TIC 454140642, minuty): primární minima obou párů, A (modře) a B (oranžově), s vnější dráhou s periodou 432 d. *Vpravo dole* (BU CMi, hodiny): primární (modře) a sekundární (oranžově) minima páru A během 30 let.](../../figures/part1/catalogue_timings.png)

**Tabulka 36.** Změny okamžiků zákrytů. „Poloviční amplituda“ je polovina rozkmitu O−C proložené
periodické křivky.

| Soustava | Veličina | SWINGERsim | Publikováno | Rozdíl | Stav |
|---|---|---|---|---|---|
| TIC 470710327 | poloviční amplituda ETV páru s periodou 1,1 d (dráha 52 d) | 4,97 min | ≈ 6 min (odečteno z obrázku) | −17 % | souhlasí |
| ξ Tau | poloviční amplituda ETV páru s periodou 7,1 d (dráha 145 d) | 19,7 min | – | – | předpověď |
| TIC 454140642 | poloviční amplitudy ETV A a B (dráha 432 d) | 51,4 a 40,2 min, ve fázi | převládá dynamická ETV (Kostov et al. 2021) | – | souhlasí |
| BU CMi | apsidální perioda A z průsečíků O−C | 27,78 roku | 28,7 roku | −3 % | souhlasí |
| BU CMi | apsidální perioda B z průsečíků O−C | 24,79 roku | 25,1 roku | −1 % | souhlasí |

**Výsledek (tabulka 36).** Na obrázku 14 přicházejí zákryty páru s periodou 1,1 d v TIC 470710327
asi ±5 min dříve a později s periodou třetí hvězdy 52 d; publikovaná hodnota, odečtená z obrázku,
je asi 6 min. ξ Tau ukazuje ±20 min s dráhou s periodou 145 d; čistý efekt světelného času této
dráhy by byl nejvýše 8,6 min, takže převládá dynamická ETV – třetí hvězda je dost blízko, aby
během každého vnějšího oběhu měnila periodu páru. Poučný je panel vlevo dole. Kdyby záleželo jen
na době šíření světla, pohybovaly by se oba páry čtyřhvězdy 2+2 v protifázi (A pozdě, když je
B brzy), zde s asi ±8 min (7,6 a 7,9 min). Spočtené ETV jsou pětkrát až sedmkrát větší a *ve
fázi*: v této kompaktní čtyřhvězdě ($P_\mathrm{out}/P_\mathrm{in}\approx$ 32–42) převládá
gravitační působení každého páru na druhý, jak popisují Kostov et al. (2021). Protifáze proto není
obecným znakem čtyřhvězdy 2+2 a srovnání s minimy skutečné kompaktní čtyřhvězdy 2+2 musí používat
úplný dynamický model, jako to dělá `swingersim timings`. V BU CMi se primární a sekundární
minima A rozcházejí a znovu protínají, jak se elipsa otáčí; doba mezi dvěma průsečíky je polovina
apsidální periody, což dává 27,78 a 24,79 roku – totéž jako měření z $\omega$ (27,76 a 24,82 roku)
nezávislou cestou.

## Cykly ZLK a dlouhodobé přežití {#catalogue-zlk}

**Metoda.** Tři soustavy se silně skloněnými drahami: σ Ori (sekulárně, 3 Myr), HD 98800
(sekulárně, 10 Myr) a α Cen s Proximou (sekulárně se slapy Galaxie, 1 Gyr).

**Tabulka 37.** Cykly ZLK a přežití. Sloupec 3 uvádí výsledek běhu, sloupec 4 publikované
pozorování.

| Soustava | Běh | SWINGERsim | Publikováno | Stav |
|---|---|---|---|---|
| σ Ori | sekulárně, 3 Myr | $e$ dráhy Aa–Ab 0,69–0,90, $i_\mathrm{mut}$ 116°–134°, perioda ZLK 39,5 kyr; nejbližší periastrum 3,4 $(R_1+R_2)$ | pozorováno $e$ = 0,78 (Schaefer et al. 2016) | souhlasí; přežije |
| HD 98800 | sekulárně, 10 Myr | $e_A$ 0,27–0,63, $e_B$ 0,42–0,80, periody ZLK 46 a 29 kyr | pozorováno 0,48 a 0,78 (Zúñiga-Fernández et al. 2026) | souhlasí; přežije |
| α Cen | sekulárně + slapy Galaxie, 1 Gyr | Proxima zůstává vázaná; $e_\mathrm{out}$ 0,50 → 0,92, periastrum až 724 au; $i_\mathrm{mut}$ 64°–83° | vázaná (Kervella et al. 2017) | souhlasí; růst $e$ je předpovědí |

![**Obrázek 15.** Cykly ZLK a slapy Galaxie: excentricita (svislé osy) v závislosti na čase. *Vlevo:* σ Ori Aa–Ab, první 0,15 běhu dlouhého 3 Myr. *Uprostřed:* dva páry HD 98800, A (modře) a B (oranžově), první 0,5 běhu dlouhého 10 Myr. *Vpravo:* dráha Proximy kolem α Cen AB během 1 Gyr.](../../figures/studies/zlk_systems.png)

**Výsledek (tabulka 37).** Na obrázku 15 (vlevo) kmitá excentricita σ Ori Aa–Ab pravidelně mezi
0,69 a 0,90 každých 40 kyr a pozorovaná hodnota 0,78 leží uprostřed kmitu; nejbližší přiblížení
zůstává 3,4násobkem součtu poloměrů, takže soustava přežije svůj věk. V HD 98800 kmitají oba páry
(0,27–0,63 a 0,42–0,80) a pozorované hodnoty 0,48 a 0,78 leží v těchto rozsazích: velké pozorované
excentricity přirozeně vysvětlují cykly ZLK, jak naznačují Zúñiga-Fernández et al. (2026). V α Cen
slapy Galaxie pomalu zvyšují excentricitu dráhy Proximy z 0,50 na 0,92 během miliardy let, aniž by
ji uvolnily; ansámbl v oddílu [Ansámbly publikovaných soustav](#catalogue-ensembles) to potvrzuje
i s procházejícími hvězdami. Samotné slapy Galaxie jsou ověřeny proti analytické maximální
excentricitě široké dvojhvězdy ([kapitola 2](#environment)).

## Rostoucí hvězdy a hvězdný vítr {#catalogue-growing}

**Metoda.** Soustavy, v nichž se hvězdy během běhu mění: HD 181068 (N-částicově během 30 let;
vývojový backend s MIST od 362,2 Myr), HD 150136 (vývojově, MIST), HD 5980 (sekulárně s větrem,
100 kyr) a PSR J0337+1715 (N-částicově během 30 let) jako kontrola přesné koplanarity.

**Tabulka 38.** Rostoucí hvězdy a hvězdný vítr.

| Soustava | Běh | SWINGERsim | Publikováno | Stav |
|---|---|---|---|---|
| HD 181068 | N-částicově, 30 let | koplanární ($i_\mathrm{mut}$ 0,8°), perioda precese uzlů 12,8 roku | 0,8 ± 1,4°; 13–14 let (Borkovits et al. 2013) | souhlasí |
| HD 181068 | vývojově, MIST | obr roste z 12,4 na 38,3 R☉ a vyplní svůj Rocheův lalok po 2,1 Myr; přenos hmoty z třetí hvězdy končí trojnásobnou společnou obálkou v čase 2,1055 Myr (stáří 364,3 Myr) | – | předpověď |
| HD 150136 | vývojově, MIST | primární složka O3 V vyplní svůj lalok po 2,45 Myr (stáří 3,45 Myr) | – | předpověď |
| HD 5980 | sekulárně s větrem, 100 kyr | $a_\mathrm{out}$ 41,3 → 42,9 au = přesně $a_0M_0/M$; hmotnosti 61 → 56,5 a 66 → 64 M☉ | Jeans (1924) | souhlasí |
| PSR J0337+1715 | N-částicově, 30 let | $i_\mathrm{mut}$ zůstává 0,00° | koplanární na 0,01° (Ransom et al. 2014) | souhlasí |

![**Obrázek 16.** Rostoucí hvězdy. *Vlevo a uprostřed:* poloměr rostoucí hvězdy (modře) a velikost jejího Rocheova laloku (čárkovaně) v závislosti na čase pro obra v HD 181068 (jehož lalok je kolem vnitřního páru) a pro primární složku O3 V v HD 150136; kde modrá čára dosáhne čárkované, začne hvězda přenášet hmotu. *Vpravo:* velká poloosa vnější dráhy HD 5980 vydělená svou počáteční hodnotou (modře) a poměr počáteční a současné celkové hmotnosti (čárkovaně) v závislosti na čase.](../../figures/studies/growing_stars.png)

**Výsledek (tabulka 38).** Na obrázku 16 (vlevo) roste obr v HD 181068 za 2,1 Myr z 12 na 38 R☉,
až se setká se svým Rocheovým lalokem o poloměru 38,3 R☉. Protože obr (3,0 M☉) je těžší než
vnitřní pár (0,92 + 0,87 M☉), je přenos na pár nestabilní (Hamers et al. 2022) a SWINGERsim ukončí
běh trojnásobnou společnou obálkou v okamžiku dotyku. HD 181068 je učebnicovým případem přenosu
hmoty z třetí hvězdy na dvojhvězdu (de Vries et al. 2014; Tokovinin & Moe 2020). Primární složka
O3 V v HD 150136 roste asi za 2,5 Myr z 9,5 na 14,6 R☉, zatímco její lalok trochu roste, protože
její vítr rozšiřuje dráhu; stane se velmi mladou hmotnou soustavou podobnou Algolu. V pravém
panelu splývají obě čáry pro HD 5980: vnější dráha se rozšiřuje přesně tak rychle, jak soustava
ztrácí hmotu, $a\propto 1/M$, což je Jeansův zákon pro pomalou ztrátu hmoty. PSR J0337+1715,
milisekundový pulzar se dvěma bílými trpaslíky, zůstává přesně koplanární, jak se pozoruje.

## Časování pulzací {#catalogue-pulsations}

**Metoda.** Světelné křivky SWINGERsimu obsahují pulzace a zpoždění se měří tak jako pozorovatel:
v krátkých oknech se proloží sinusoida s frekvencí pulzací a sleduje se její fáze
([kapitola 4](#pulsations)). Tři soustavy: λ Sco (pulzující hvězda typu β Cep, 1100 d, pásmo V),
Polárka (cefeida pulzující v prvním svrchním tónu s pozorovanou periodou 3,97 d, 30 let, V)
a KIC 4150611 (δ Sct, 200 d, pásmo TESS).

**Tabulka 39.** Časování pulzací. Sloupec 3 je poloviční amplituda zpoždění změřená ze syntetické
světelné křivky, sloupec 4 poloviční amplituda efektu světelného času dráhy použité v běhu (Irwin
1952), sloupec 5 publikovaná hodnota.

| Soustava | Frekvence pulzací | Změřeno | Dráha běhu | Publikováno | Rozdíl | Stav |
|---|---|---|---|---|---|---|
| λ Sco | 4,679 c/d | 1067 s | 1069 s | 1193 s | −11 % | vysvětleno: vstup |
| Polárka | 0,252 c/d | 1309 s | 1336 s | 1335 s | −2 % | souhlasí |
| KIC 4150611 | 20,243 c/d | 97,1 s | 97,9 s | 98 ± 9,5 s | −1 % | souhlasí |

![**Obrázek 17.** Časování pulzací ve třech soustavách katalogu. *Nahoře:* několik dní každé syntetické světelné křivky (relativní tok v závislosti na čase); pravidelné vlny jsou pulzace, u λ Sco také dva zákryty. *Dole:* zpoždění fáze pulzací v sekundách (body), změřené v krátkých oknech podél celé světelné křivky, s křivkou světelného času dráhy, která nese pulzující hvězdu (čára); čárkované čáry vyznačují ± publikovanou poloviční amplitudu. Text nad každým panelem uvádí změřenou, modelovou a publikovanou hodnotu.](../../figures/part1/pulsator_delays.png)

**Výsledek (tabulka 39).** Na obrázku 17 sledují změřená zpoždění (body) křivku světelného času
vlastní dráhy běhu (čára) ve všech třech soustavách s přesností 2 %: měření je přesné. Zbývající
rozdíly vůči publikovaným číslům pocházejí ze vstupu. Pro λ Sco používá konfigurace
interferometrickou dráhu a hmotnosti podle Tanga et al. (2006), které dávají $a_1\sin i$ = 2,15 au;
publikovaných 1193 s patří ke starší dráze podle Uytterhoevenové et al. (2004). U Polárky je často
citovaných 1430 s hodnotou $a_1\sin i/c$; poloviční amplituda křivky světelného času je menší
o faktor $\sqrt{1-e^2\cos^2\omega}$, 1335 s (Evans et al. 2024), což SWINGERsim reprodukuje na 2 %
(spočten je jen jeden oběh dlouhý 30 let, takže křivka není vystředěna na nulu; její výška od
nejnižšího k nejvyššímu bodu je dvojnásobkem poloviční amplitudy). Perioda prvního svrchního tónu,
kterou model hvězdy Polárky předpovídá ze vztahu perioda–svítivost–hmotnost, je 3,62 d proti
pozorovaným 3,97 d (−9 %), v rámci rozptylu takových vztahů; běh používá pozorovanou periodu. Pro
KIC 4150611 bere konfigurace hmotnost průvodce $M_\mathrm{Ab}$ = 0,86 M☉, v souladu s poměry
velikostí a period podle Kempa et al. (2024); jejich dynamický poměr hmotností 3,61 těmto poměrům
odporuje (ty implikují 1,70) a dal by jen 54 s. S konzistentními hmotnostmi dává dráha publikovaných
98 s (Hełminiak et al. 2017 z radiálních rychlostí; Shibahashi & Kurtz 2012 z frekvenční modulace).
Časování pulzací je tedy citlivým testem hmotností v konfiguraci.

## Ansámbly publikovaných soustav {#catalogue-ensembles}

**Metoda.** Kde je parametr skutečné soustavy neznámý (vnější dráha, uzel, vzájemný sklon), losuje
jej ansámbl z apriorního rozdělení a ptá se, co pozorovaná soustava implikuje.

**Tabulka 40.** Ansámbly publikovaných soustav. Sloupec 3 je počet členů, sloupec 4 otázka,
sloupec 5 odpověď.

| Soustava | Backend | Členy | Otázka | Odpověď |
|---|---|---|---|---|
| HD 150136 | N-částicový | 500 | záleží na neznámém vzájemném sklonu během století? | $i_\mathrm{mut}$ se změní o < 0,04°; pár s periodou 2,7 d zůstává zákrytový; jeho sklon se posune o 0–0,57° za století podle $i_\mathrm{mut}$ |
| HD 5980 | N-částicový | 1000 | je soustava 2+2 vázaná a bezpečná? | 29,7 % se srazí (medián největší $e_B$ 0,91), většinou pro vnější periody pod 20 let |
| GG Tau | N-částicový | 500 | osud mladé pětihvězdy během 20 kyr | neporušena 80,8 %, přestavba 9,2 %, srážka 7,6 %, vymrštění 2,4 % |
| TIC 168789840 | sekulární | 500 | jak často zakrývá každý pár? | A 31 %, B 91 %, C 35 % času; všechny tři společně 11 % |
| KIC 4150611 | sekulární | 500 | které vnější dráhy udrží pár Ab kruhový a je soustava stabilní? | s konzistentními hmotnostmi zůstávají všechny členy neporušené a stabilní podle MA01 a Vynatheyi et al. (2022); pro neznámé vnější dráhy si 42,6 % udrží $e$ < 0,05 během 1 Myr, ostatní jsou vybuzeny až na 0,7 |
| α Cen | sekulární + slapy Galaxie + procházející hvězdy | 300 | zůstane Proxima vázaná 1 Gyr? | 100 % vázaných; medián největší $e$ 0,90; medián nejbližšího periastra 849 au, minimum 70 au |

![**Obrázek 18.** Ansámbly publikovaných soustav, každý v jednom panelu. *HD 5980:* podíl členů, které skončí srážkou, v závislosti na vnější periodě. *GG Tau:* osudy (modře neporušena, oranžově srážka, zeleně přestavba, fialově vymrštění) pro rozsahy velké poloosy Aa–Ab. *TIC 168789840:* střední podíl času, kdy zakrývá každý ze tří párů. *KIC 4150611:* každý bod je člen, největší excentricita páru Ab v závislosti na velikosti vnější dráhy; čárkovaně: $e$ = 0,05. *α Cen:* histogram největší excentricity dráhy Proximy během 1 Gyr. *HD 150136:* každý bod je člen, rychlost změny sklonu zákrytového páru v závislosti na neznámém vzájemném sklonu.](../../figures/part1/catalogue_ensembles.png)

**Výsledek (tabulka 40).** Obrázek 18 ukazuje šest ansámblů. U HD 5980 se většina členů s vnějšími
periodami pod 20 let srazí, dlouhé periody jsou bezpečné: tvoří-li obě dvojhvězdy vázanou
čtyřhvězdu 2+2, je vnější perioda pravděpodobně delší než 20 let – dynamické omezení pro
spektroskopiky. GG Tau zůstává neporušená asi v 80 % členů všude; 74 % členů, které MA01 označuje
za nestabilní, zůstane 20 kyr neporušených, což jen ukazuje, že 20 kyr je krátká doba ve srovnání
s časovou škálou MA01 (MA01 označuje za stabilní 34,8 %). V TIC 168789840 zakrývá pár B téměř
vždy, A a C třetinu času: trojnásobné zákryty pozorované TESS jsou shodou okolností s pravděpodobností
1 : 9. Členy KIC 4150611 buď udrží pár Ab kruhový, nebo jej při jakékoli velikosti vnější dráhy
vybudí na 0,6–0,7. Proxima zůstává vázaná ve všech 300 členech, většinou s největšími
excentricitami 0,85–0,97 (Kervella et al. 2017 ji shledávají vázanou). U HD 150136 posun sklonu
mizí při vzájemných sklonech 0°, 90° a 180°, takže změřená změna hloubky zákrytu by určila vzájemný
sklon; na pozorovatelných časových škálách jsou efekty zanedbatelné, jak uvádějí Mahy et al.
(2018).

## Shrnutí katalogu {#catalogue-summary}

Všech 23 soustav proběhlo bez selhání. Každá veličina, kterou lze srovnat s publikací – apsidální
a uzlové periody, amplitudy změn okamžiků zákrytů, rozsahy ZLK, koplanarity, amplitudy efektu
světelného času, velikosti Rocheových laloků a geometrie disků – souhlasí v rámci publikované
nejistoty nebo asi 10 %, nebo byl rozdíl vysledován ke vstupu konfigurace (δ Ori A, λ Sco) nebo
k odlišné definici publikované veličiny (VW LMi, TIC 120362137, Polárka). Obrázek se všemi
srovnáními na jedné stupnici je v [kapitole 11](#conclusions).
