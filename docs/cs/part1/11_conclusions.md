# Závěry {#conclusions}

## Co je ověřeno a s jakou přesností {#validated}

SWINGERsim obsahuje 219 srovnávacích testů, každý se svým zdrojem a tolerancí stanovenou před
srovnáním. S konečnou verzí kódu jich 205 vyhovuje (PASS) a 14 leží mimo svou toleranci
z doloženého důvodu, který není chybou SWINGERsimu; žádný neselhává. Tři čtvrtiny z nich (165)
ověřují kód napsaný přímo pro SWINGERsim, nejen externí N-částicový integrátor. Tabulka 47 shrnuje
přesnost dosaženou pro každou skupinu fyziky; jednotlivá čísla jsou v kapitolách 2–4
a v [příloze A](#appendix-benchmarks).

**Tabulka 47.** Ověřená fyzika. Sloupec 2 jmenuje reference, sloupec 3 dosaženou přesnost (největší
rozdíl mezi vyhovujícími testy skupiny, nebo typický rozdíl, je-li tak uvedeno), sloupec 4 třídu
ověřovaného kódu (tabulka 1).

| Fyzika | Reference | Přesnost | Třída |
|---|---|---|---|
| newtonovská a 1PN N-částicová dynamika | Merkur, limity ZLK, Söderhjelm, MA01, MEGNO | 10⁻⁴ (relativita), 0,3 % (ZLK, s očekávaným krátkoperiodickým přebytkem) | R |
| sekulární backendy (trojhvězda, 2+2, vnořený) | analytické limity, Naoz et al. (2013), Vokrouhlický (2016), MSE, N-částicový běh šestihvězdy | 10⁻⁴ v analytických limitách; 0,04° na publikovaných křivkách; 10⁻⁴ v $e_\mathrm{max}$ proti MSE; 1,2 % v periodě ZLK proti N-částicovému běhu | own |
| ansámbly chaotických soustav | Pejcha et al. (2013), Antognini & Thompson (2016) | podíly pro trojhvězdy a zesílení souhlasí; podíly pro čtyřhvězdy závisejí na délce běhu (viz níže); účinné průřezy 1–2 % (výměna 23 %, 1,8σ) | R, own |
| slapy a apsidální stáčení | Sterne (1939), Hut (1981), Leconte et al. (2010), CzeV343 | 10⁻⁴ (apsidální), 0,3 % (disipativní rychlosti), 3 % (pozorovaná apsidální perioda) | R, own |
| strukturní konstanty hvězd a vývojové dráhy | Brooker & Olle, Claret, Tout et al., MIST, ženevské modely, Slunce, bílí trpaslíci | poslední číslice tabulek; 5 % pro Slunce; 4–12 % mezi rodinami modelů | own/tab |
| inerciální vlny a zakulacování drah před hlavní posloupností | Ogilvie, Barker, Zahn & Bouchet, Meibom & Mathieu | 0,03 % (řešič), 1 % (fáze před hlavní posloupností) | own |
| ztráta hmoty | Jeans, Hadjidemetriou, Veras et al., Blaauw/Hills, Perets & Kratter | 4·10⁻⁹ (adiabatický invariant), 0,24 % (impulzní), 6,5 % (počátek nestability) | R, own |
| přenos hmoty a společná obálka | Soberman et al., Postnov & Yungelson, BSE, COSMIC, Tout et al. | 3·10⁻¹³ (uzavřené tvary); 0,02 M☉ u obnaženého dárce; 0,5–9 % v periodách po společné obálce | own |
| magnetismus, Applegateův mechanismus | Matt et al. (2015), Völschow et al. (2016) | 0,7 % (rotace Slunce), 6 % (energetické poměry) | own |
| okolí a gravitační vlny | Grishin & Perets, Hamers & Samsing, PSR B1913+16 | 3·10⁻⁵ (slapy Galaxie), 0,7 % (procházející hvězdy), 0,02 % (zkracování dráhy) | own, R |
| světelné křivky na koulích a sítích, kontaktní dvojhvězdy | batman, ellc, PHOEBE 2.5 | 10⁻⁸ (zákryt), 4·10⁻⁴ (PHOEBE, koule), 0,003 mag (kontaktní dvojhvězda) | own |
| prachové disky | přesná geometrie, Farago & Laskar, ε Aur, HD 98800 | přesná geometrie; 14 % v trvání zákrytu ε Aur (ostrý okraj); 2,5° u disku HD 98800 | own |
| okamžiky zákrytů a proložení v kroku 1 | Irwin, Kaplan, Giménez & Bastero, OCFit | 10⁻⁴ (LTTE), 0,3–0,6 % (Rømerovo zpoždění), 0,2 % ($a\sin i$ z N-částicových minim) | own |
| časování pulzací | Murphy et al. (2014), vztahy pro δ Sct a cefeidy | 3·10⁻⁷ (analyticky), 1,3 % (KIC 9651065), 2 % (pulzující hvězdy katalogu) | own |
| klasifikátory stability | Vynatheya et al. (2023), Quarles et al. (2018) | skóre MLP na 0,01; MA01 na 0,05 se vzorkováním článku | own, ext |

Kromě srovnávacích testů SWINGERsim reprodukuje publikovaný model QZ Car (velikosti drah na
0,07 %, podíly světla na 0,001), amplitudy precese, změn okamžiků zákrytů a časování pulzací soustav
katalogu (v rámci publikované nejistoty nebo asi 10 % všude, kde je vstup konzistentní)
a pozorované fáze ZLK soustav σ Ori a HD 98800. Přinesl nové, ověřitelné výsledky: hierarchie
2+(2+1) soustavy CzeV4315 je preferována před 2+2 asi faktorem 70 a slapy snižují podíl jejích
selhání na polovinu; čtyřhvězda 2+2 vysvětluje ASASSN-V J075203 jen v úzkém ostrově pravidelného
pohybu se dvěma předpovězenými projevy; změny okamžiků zákrytů kompaktní čtyřhvězdy 2+2
TIC 454140642 probíhají ve fázi; HD 181068 končí trojnásobnou společnou obálkou; a precese disků
EE Cep a KH 15D omezuje jejich velikost.

![**Obrázek 24.** Každé kvantitativní srovnání SWINGERsimu s publikovanou hodnotou v této zprávě na jedné stupnici. Každý řádek je jedno číslo; vodorovná poloha bodu je hodnota SWINGERsimu dělená publikovanou hodnotou (logaritmická osa; 1 = stejné číslo, 0,5 = polovina, 2 = dvojnásobek). Šedý pás vyznačuje ±10 %; vodorovné úsečky jsou publikované nejistoty 1σ dělené publikovanou hodnotou. Modře: souhlasí; oranžově: rozdíl byl vysledován k jeho příčině. Sloupec vpravo uvádí obě čísla, hodnotu SWINGERsimu jako první.](../../figures/part1/literature_overview.png)

Obrázek 24 dává všechna srovnání celých soustav a publikovaných simulací na jednu stupnici. Čtěte
jej řádek po řádku: modrý bod uvnitř šedého pásu je shoda lepší než 10 %; modrý bod mimo pás
s úsečkou, která dosahuje svislé čáry, je shoda v rámci publikované chyby. Z 38 řádků je 28 modrých.
Každý z 10 oranžových řádků má vysledovanou příčinu: vstup konfigurace nebo článku (δ Ori A, λ Sco,
hvězdy QZ Car, které získaly hmotu), publikovaná čísla, která si odporují (Algol z BSE, V994 Her C),
délku publikovaných běhů (Pejcha et al.) nebo otevřenou otázku slapového modelu SWINGERsimu (M67).

## Čtrnáct vysvětlených rozdílů {#explained}

**Tabulka 48.** 14 srovnávacích testů s verdiktem EXPLAINED. Sloupec 2 uvádí hodnotu SWINGERsimu,
sloupec 3 publikovanou hodnotu, sloupec 4 důvod, sloupec 5 druh důvodu.

| Test | SWINGERsim | Publikováno | Důvod | Druh |
|---|---|---|---|---|
| `v994s-gr-C` | 0,081 | 0,06 | vlastní tabulka 3 článku dává 8,1 %; A a B, spočtené stejně, souhlasí | přehlédnutí v článku |
| `mt-bse-algol-accretor` | 2,51 M☉ | 2,85 M☉ | publikovaná historie by vyžadovala, aby dráha v Hertzsprungově mezeře získala moment hybnosti; konzervativní fáze obra je reprodukována na 1,2 % | publikovaná čísla si odporují |
| `mt-bse-algol-period` | 21,9 d | 30 d | jako výše | publikovaná čísla si odporují |
| `qzcar-evolutionary-mass-ac2` | 26,0 M☉ | 33,2 M☉ | Ac2 pravděpodobně získala hmotu (její proložená rotace je 410 km/s); vývojové dráhy osamocených hvězd ji popsat nemohou | vstup / definice |
| `qzcar-evolutionary-mass-aa1` | 42,5 M☉ | 69,8 M☉ | Brož et al. naznačují dřívější přenos Aa2 → Aa1; omlazená hvězda, která získala hmotu, se jeví mladší a lehčí | vstup / definice |
| `vok16-pejcha-fA` | 0,21 | 0,31 | podíly rostou s délkou běhu; při 50 kyr jsou 0,36; tři nezávislé integrace souhlasí | délka publikovaných běhů |
| `vok16-pejcha-fB` | 0,40 | 0,481 | jako výše (0,57 při 50 kyr) | délka publikovaných běhů |
| `vyn23-3p1-ma01-score` | 0,79 | 0,56 | skóre závisí na vzorkování; vlastní tabulka 4 článku implikuje 0,69; při losování blízko hranice jako v článku: 0,51 | vzorkování publikované množiny |
| `tidal-capture-radius-n3` | 2,37 | 2,72 | citované poloměry patří k jiným konvencím; tytéž aproximace souhlasí s McMillanem et al. (1987) | definice |
| `tidal-capture-radius-n15` | 3,33 | 3,84 | jako výše (3,33 proti asi 3,4) | definice |
| `tides-circularization-m67` | 3,8 d | 12,1 d | samotný rovnovážný slap je na hlavní posloupnosti příliš slabý; Barker (2022) vysvětluje M67 inerciálními vlnami, které dvouvrstvý odhad SWINGERsimu se stavbou SSE/MIST a vyvíjejícími se rotacemi nereprodukuje (stavba hvězd a historie rotace) | otevřená otázka slapového modelu |
| `tides-circularization-m67-iw` | 7,5 d | 12,1 d | jako výše, s inerciálními vlnami | otevřená otázka slapového modelu |
| `tides-circularization-m67-final` | 7,5 d | 12,1 d | jako výše, s fází před hlavní posloupností (ta je reprodukována: 6,6 d proti 7,1 ± 1,2 d) | otevřená otázka slapového modelu |
| `eecep-disc-precession` | 183 let | 65 let | tuhý disk ilustrativní velikosti se stáčí příliš pomalu; 65 let vyžaduje vnější poloměr asi 1,6 au | omezení disku |

Žádný ze čtrnácti rozdílů není chybou SWINGERsimu; tři testy M67 zůstávají otevřenou otázkou jeho
slapového modelu. Čtyři jsou problémy publikovaných čísel (přehlédnutí, vnitřně nekonzistentní
historie, která se týká dvou testů, a skóre, které vlastní tabulka článku nepodporuje), dva se
týkají hvězd, které nejsou modely osamocených hvězd, dva odlišné definice, dva délky publikovaných
běhů, tři otevřené otázky slapového modelu SWINGERsimu a jeden mění rozdíl v omezení objektu.

## Co ověřeno není {#not-validated}

Následující části SWINGERsimu nebo režimy jeho fyziky nemají v této zprávě nezávislou referenci.
Je třeba je s tímto vědomím používat.

- **Slapové zakulacování drah na hlavní posloupnosti:** zakulacení M67 do 12 d SWINGERsim
  nereprodukuje. Barker (2022) je vysvětluje inerciálními vlnami spočtenými z úplných profilů hvězd;
  dvouvrstvý odhad SWINGERsimu středovaný přes frekvence se stavbou SSE a MIST dává 7,5 d. Rozdíl
  spočívá ve stavbě hvězd a v historii rotace a je otevřenou otázkou.
- **Proudy při přenosu hmoty** jsou předpisy (rychlosti BSE, akreční podíly podle Hamerse et al.
  2022), nikoli hydrodynamika; akreovaný podíl během přenosu na tepelné časové škále je nejméně
  jistým číslem vývoje Algolů a je volbou uživatele. Přenos hmoty z třetí hvězdy a trojnásobné
  společné obálky nemají kvantitativní srovnávací test; ověřují se zákony zachování a očekávaným
  výsledkem.
- **Pozdní fáze velmi hmotných hvězd:** od pozdní hlavní posloupnosti dále mřížka MIST nerozlišuje
  nafouknuté obálky hvězd 40–70 M☉ (chyby poloměru až 20–40 %) a hmotnosti na konci hlavní
  posloupnosti se mezi rodinami modelů liší až o třetinu.
- **Kontaktní dvojhvězdy** mají jednu společnou teplotu a žádné oscilace tepelné relaxace; splynutí
  je ověřeno jen proti očekávanému chování.
- **Prachové disky** gravitačně nepůsobí na dráhy a hvězda uvnitř disku nepřidává své světlo.
- **Atmosféry pod 3500 K** a s $\log g$ < 0 používají nejbližší uzel tabulek okrajového ztemnění.
- **Klasifikátory stability pro více než čtyři hvězdy** používají sítě pro 2+2 a 3+1 na dílčí
  čtyřhvězdy, což je přiblížení mimo jejich trénovací množiny; referencí jsou algebraická kritéria
  pro každou vnořenou trojici a test duchem.
- **Dosud neprovedená srovnání:** rezonance 1 : 1 dvojhvězd čtyřhvězdy 2+2 (Breiter & Vokrouhlický
  2018), sekulární chaos v soustavách 2+2 a 3+1 (Hamers & Lai 2017), osudy hierarchických soustav
  s vývojem hvězd proti MSE, podíly trojhvězd destabilizovaných větrem a přenosem hmoty (Toonen et al.
  2022) a širší vzorek pozorovaných soustav, jako je Multiple Star Catalog (Tokovinin 2018).

Záměrně mimo rozsah SWINGERsimu jsou relativistické časování pulzarů, kompaktní pozůstatky
a impulzy při výbuších supernov, Wolfovy–Rayetovy hvězdy až do kolapsu jádra, srážející se
hvězdné větry, skvrny a erupce ve světelných křivkách, neradiální pulzační módy a výběrové efekty
přehlídek.

## Shrnutí {#summary}

Každý fyzikální modul SWINGERsimu reprodukuje nezávislé publikované výsledky: analytické limity
s přesností 10⁻⁴ nebo lepší, publikované simulace a jiné programy na zlomek stupně nebo několik
procent a pozorované soustavy v rámci jejich nejistot. Kde rozdíl zůstává, byl vysledován
k publikovanému číslu, jeho definici nebo délce publikovaných běhů, nebo je uveden jako otevřená
otázka (zakulacování drah v M67 na hlavní posloupnosti). Většina těchto důkazů se týká kódu
napsaného pro SWINGERsim – sekulárních backendů, vazby vývoje, slapů a přenosu hmoty na dráhy,
světelných křivek, okamžiků zákrytů, prokládání a ansámblů. Části II–IV ověřují tentýž kód na
skutečných pozorováních soustav SQUADRA.
