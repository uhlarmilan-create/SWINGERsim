# Soustavy SQUADRA {#squadra-systems}

SQUADRA je přehlídka kandidátů na čtyřhvězdy 2+2, kterou vedou J. Kolář, P. Zasche
a spolupracovníci: jde o soustavy, v nichž TESS vidí dva zákrytové páry v obrazu jedné hvězdy
(Kolář et al. 2025: 781 kandidátů). Pro SWINGERsim je otázkou, zda lze pozorované chování – hloubky
zákrytů, které se mění, a okamžiky zákrytů, které se odchylují od konstantní periody – vysvětlit
dynamickým modelem a která hierarchie je správná. Tato kapitola nereprodukuje publikovaný model,
protože pro tyto soustavy žádný neexistuje; ověřuje aparát ansámblů SWINGERsimu na skutečných
objektech: podmínění pozorováními, evidenci soupeřících modelů, shodu sekulárního
a N-částicového backendu a symetrii dynamiky v čase. Podrobné analýzy soustav SQUADRA se
skutečnými minimy a světelnými křivkami jsou předmětem částí II–IV.

**Definice.** V ansámblu dostane každý člen **váhu** – věrohodnost pozorování při jeho
parametrech, například $\exp(-\chi^2/2)$ pro změřenou rychlost změny s její chybou. Střední váha přes
všechny členy losované z apriorního rozdělení je **evidence** $Z$ modelu; poměr evidencí dvou modelů
(Bayesův faktor) říká, který z nich vysvětluje data lépe. **Efektivní počet členů**
$n_\mathrm{eff} = (\sum w)^2/\sum w^2$ (Kish 1965) udává, kolik členů skutečně nese váhu;
$n_\mathrm{eff}$ blízké 1 znamená, že vše nese jediný člen a evidence není spolehlivá. **Osud** člena
říká, jak jeho běh skončil: neporušen, srážka, vymrštění, přestavba (změna hierarchie), výměna nebo
rozpad.

## CzeV4315: 2+2, nebo 2+(2+1)? {#czev4315}

**Soustava.** CzeV4315 ukazuje dva zákrytové páry: A (6,74 d, excentrický, $e$ = 0,64), jehož
zákryty jsou léta stabilní, a B (0,92 d), jehož sklon se mění o −1,4 ± 0,2 °/rok – jeho zákryty
rychle mělčí. Pozoruje se změna okamžiků zákrytů s periodou asi 119 d. Jde o prostou čtyřhvězdu
2+2, nebo o pětihvězdu 2+(2+1), v níž má B vlastního blízkého průvodce na dráze s periodou 119 d?

**Vstup.** Ansámbl 2+2 s 20 000 členy se sekulárním backendem (±300 let; vnější perioda
logaritmicky rovnoměrně mezi 40 a 3000 d, izotropní orientace, odhadnuté hmotnosti) a ansámbl
2+(2+1) s 10 000 členy s N-částicovým backendem a relativitou (±60 let; průvodce B na dráze
s periodou 119 ± 3 d). Pozorování použitá jako podmínky: rychlost změny sklonu B, žádné přetečení
Rocheova laloku a volitelně to, že A zakrývala nejméně 5 let v minulosti a bude zakrývat nejméně
0,25 roku dopředu.

**Tabulka 29.** Evidence dvou modelů CzeV4315. Sloupec 3 uvádí pozorování, jimiž je ansámbl
podmíněn; sloupce 4–5 evidenci a efektivní počet členů.

| Model | Členy | Podmínky | Evidence $Z$ | $n_\mathrm{eff}$ |
|---|---|---|---|---|
| 2+2, široké apriorní rozdělení vnější periody | 20 000 | rychlost změny sklonu B | 0,0063 | 177 |
| 2+2, široké apriorní rozdělení vnější periody | 20 000 | + A zakrývá po léta | 0,00013 | 4,8 |
| 2+(2+1) | 10 000 | sklon B + A zakrývá + oddělená | 0,0104 | 165 |
| 2+(2+1) | 10 000 | totéž + neporušena v obou směrech času | 0,0095 | 153 |

**Výsledek.** Tabulka 29 ukazuje, že pětihvězda 2+(2+1) je za stejných podmínek preferována před
čtyřhvězdou 2+2 asi faktorem 70. Nejvíce rozlišuje pozorování, že A zakrývá po léta: odstraní
téměř všechny členy 2+2. Členy 2+2, které data prokládají, leží na hraně stability (vážený medián
poměru MA01 1,43, vnější perioda asi 116 d). Zde sekulárnímu backendu, který předpokládá stabilní
hierarchii, nelze věřit samotnému. N-částicová kontrola 40 nejlepších členů 2+2 během 100 let
s MEGNO nachází 4 srážky, 14 chaotických členů (MEGNO 4–35, Ljapunovovy doby asi 1,5–25 let) a 22
pravidelných: asi polovina řešení 2+2 nemůže popisovat soustavu, která dnes existuje, takže
evidence 2+2 je nadhodnocená a preference 2+(2+1) je konzervativní.

![**Obrázek 9.** Nejlepší konfigurace 2+2 soustavy CzeV4315 v sekulárním ansámblu: srážkový parametr obou párů během 600 let. Vodorovná osa: čas od 300 let v minulosti do 300 let dopředu (čárkovaná čára: dnešek). Svislá osa: nejmenší promítnutá vzdálenost obou hvězd páru v konjunkci v jednotkách součtu jejich poloměrů; pod 1 (vystínovaný pás) pár zakrývá. Modře: pár A; oranžově: pár B.](../../figures/showcase/czev4315.png)

Obrázek 9 ukazuje, proč to má čtyřhvězda 2+2 těžké. Pár A (modře) zůstává hluboko pod 1 – zakrývá
neustále. Pár B (oranžově) kmitá nahoru a dolů, jak se rovina jeho dráhy stáčí; každých několik
desetiletí překročí 1, takže jeho zákryty na léta mizí a vracejí se (B zakrývá 80 % času, v oknech
asi 20 let). Dnes je B na cestě dolů – to je pozorovaná rychlá změna hloubky jeho zákrytů.
Čtyřhvězda 2+2 reprodukuje tuto rychlost změny sklonu jen s vnějšími periodami 50–170 d (medián
116 d), blízko meze stability; v modelu 2+(2+1) budí precesi B přímo jeho průvodce s periodou 119 d.

**Symetrie v čase.** Konzervativní dynamika je symetrická v čase: konfigurace, která je nestabilní,
to ukáže v minulosti stejně často jako v budoucnosti. Ve 20 000 bězích ansámblu 2+(2+1) (10 000
členů v obou směrech) selže 5,99 % členů dopředu a 5,97 % dozadu. V obou směrech v rámci okna
selže 3,93 %, kde by nezávislá selhání dala 0,36 %: nestabilita je vlastností konfigurace a počáteční
fáze rozhoduje jen o tom, kdy se projeví. Kumulativní podíl selhání roste bez nasycení (0,4 % při
20 letech, 2,4 % při 40 letech, 6,0 % při 60 letech), jak se očekává blízko hranice stability, kde
doba rozpadu roste s hierarchií zhruba exponenciálně (Hayashi et al. 2022; Mushkin & Katz 2020).
4,2 % běhů, které končí srážkou, k ní dospěje třemi cestami: cykly ZLK dvojhvězdy B buzenými jejím
průvodcem (vzájemné sklony 60°–120°; B dosáhne $e\approx0{,}45$, kdy se její hvězdy dotknou),
excentricitou průvodce zvýšenou párem A na vnější dráze asi na 0,9, takže jeho periastrum se
přiblíží na 1–3 velké poloosy B a vnitřní trojhvězda se stane nestabilní, a – zřídka – růstem
excentricity A z 0,64 na asi 0,86.

## CzeV4315 se slapy a rotací {#czev4315-tides}

**Otázka.** Ansámbl 2+(2+1) používal hmotné body s relativitou. Skutečné hvězdy jsou deformovány
slapy a rotací; deformace přidává rychlé apsidální stáčení a disipuje energii. S apsidálními
konstantami těchto hvězd (Claret 2023) má slapové a rotační stáčení B periodu asi 9 let, kratší
než periody ZLK 14–118 let. Změní to, které členy přežijí?

**Vstup.** Stejné apriorní rozdělení a stejné náhodné semínko jako předtím, takže členy 0–4999 jsou
přesně tatáž losování, plus slapy a rotace; N-částicově, 100 let dopředu a dozadu (10 000 běhů).

**Tabulka 30.** CzeV4315 se slapy a rotací: osudy 5000 členů po 100 letech v každém směru.

| Osud (100 let) | Dopředu | Dozadu |
|---|---|---|
| neporušena | 94,68 % | 94,34 % |
| srážka | 2,50 % | 2,52 % |
| vymrštění | 1,26 % | 1,38 % |
| přestavba | 0,92 % | 1,12 % |
| výměna | 0,36 % | 0,38 % |
| rozpad | 0,28 % | 0,26 % |
| neporušena ne v obou směrech | 4,20 % (0,30 % při nezávislosti) | |

![**Obrázek 10.** CzeV4315 se slapy a rotací. *Vlevo:* pro každý osud kromě neporušené soustavy procento z 5000 členů, které tak skončily během 100 let; modré sloupce: běhy do budoucnosti, oranžové sloupce: běhy do minulosti. *Uprostřed:* procento týchž členů, které nezůstaly neporušené během prvních 60 let, s hmotnými body a relativitou (tmavě šedě) a s přidanými slapy a rotací (zeleně), dopředu a dozadu. *Vpravo:* histogram (obě osy logaritmické) největší excentricity, které pár B s periodou 0,92 d dosáhl v každém běhu, s hmotnými body (šedě vyplněno) a se slapy a rotací (zelený obrys).](../../figures/part1/job010_fates.png)

**Výsledek.** Obrázek 10 a tabulka 30 ukazují tři věci. V levém panelu mají modré a oranžové sloupce
pro každý osud téměř stejnou délku (srážky 2,50 % dopředu, 2,52 % dozadu): symetrie v čase platí
i se slapy. V prostředním panelu klesne podíl selhání během prvních 60 let z 6,2 % na 3,3–3,6 %,
když se přidají slapy a rotace; 223 běhů, které se s hmotnými body srazily, zůstane neporušených.
Pravý panel ukazuje proč: s hmotnými body je excentricita B zvyšována na jakoukoli hodnotu až do 1
(medián 0,45 v bězích, které se srazily), kdežto se slapy zůstávají téměř všechny běhy blízko
$10^{-3}$ (medián 0,0015) – slapové a rotační výdutě způsobují tak rychlé stáčení dráhy B, že její
průvodce už nemůže zvyšovat její excentricitu. Toto potlačení cyklů ZLK silami krátkého dosahu je
známé (Fabrycky & Tremaine 2007; Liu, Muñoz & Lai 2015); zde se objevuje u skutečného kandidáta.
Při opakování se slapy během 100 let přežijí příklady všech tří cest ke srážce, kdežto hluboký pád
průvodce (periastrum ve vzdálenosti 1,6 $a_B$) ke srážce stále vede: síly krátkého dosahu odstraní
sekulární srážky a část dynamických. Praktický závěr: hmotné body nadhodnocují selhání kompaktních
soustav 2+(2+1) dvojnásobně a ansámbly takových soustav potřebují slapy a rotaci. CzeV4315 je se
svými skutečnými minimy analyzována v části IV.

## ASASSN-V J075203 {#j075203}

**Soustava.** U tohoto kandidáta na čtyřhvězdu 2+2 se B naklání o 1,03 ± 0,2 °/rok, A si udržuje
konstantní sklon (0 ± 0,3 °/rok) a A vykazuje rychlé apsidální stáčení (apsidální perioda
$U\approx3{,}5\pm1{,}5$ roku). Dokáže čtyřhvězda 2+2 reprodukovat všechny tři současně?

**Tabulka 31.** ASASSN-V J075203: čtyři kroky analýzy. Sloupec 3 uvádí, co každý krok zjistil.

| Krok | Členy | Výsledek |
|---|---|---|
| izotropní apriorní rozdělení, sekulárně | 20 000 | žádný člen neprokládá všechna tři pozorování (nejlepší $\chi^2$ 6,2; $n_\mathrm{eff}$ = 1,0) |
| vnější dráha blízko roviny A ($i\sim N(80{,}7°, 3°)$, $P_\mathrm{out}$ 40–600 d), sekulárně ±100 let | 40 000 | $n_\mathrm{eff}$ = 53; 32 členů s $\chi^2$ < 3; $P_\mathrm{out}$ 50–80 d, A téměř koplanární s vnější dráhou |
| N-částicově ±30 let pro 32 členů | 32 | 2 se srazí; 7 prokládá všechna tři ($\chi^2$ < 3,5, $P_\mathrm{out}$ 55–83 d) |
| N-částicově s MEGNO ±1000 let pro 7 členů | 7 | 4 pravidelné, 1 chaotický (Ljapunovova doba asi 13,5 roku), 2 se srazí během 50–150 let |

**Tabulka 32.** Sedm členů, které N-částicově prokládají všechna tři pozorování. Sloupce: číslo
člena, vnější perioda, poměr stability MA01 (nad 1 = stabilní podle kritéria), MEGNO
(2 = pravidelný) a verdikt.

| Člen | $P_\mathrm{out}$ [d] | Poměr MA01 | MEGNO | Verdikt |
|---|---|---|---|---|
| 642 | 55,3 | 1,06 | 2 | pravidelný |
| 12744 | 62,3 | 1,16 | 2 | pravidelný |
| 39269 | 66,3 | 1,22 | 2 | pravidelný |
| 38234 | 77,5 | 1,21 | 2 | pravidelný |
| 1043 | – | 1,08 | 40–55 | chaotický, neporušený 1000 let |
| 2654, 30610 | – | 0,87–0,89 | – | srazí se, vyloučeny |

**Výsledek.** Tabulky 31 a 32 sledují hledání. S izotropním apriorním rozdělením nese celou váhu
jeden člen ($n_\mathrm{eff}$ = 1,0): evidence takového ansámblu nemá smysl a apriorní rozdělení je
třeba zúžit, což dělá druhý krok. Členy, které data prokládají, leží na mezi stability (MA01
0,9–1,2, poměr period $P_\mathrm{out}/P_A\approx$ 6–9). Tam středování přes oběh v sekulárním
backendu selhává: N-částicová apsidální perioda A je asi poloviční oproti sekulární, jak se
očekává při malém poměru period (Luo, Katz & Dong 2016). Všechny testy přežijí jen 4 pravidelné
členy s vnějšími periodami 55–78 d. Čtyřhvězda 2+2 je proto možná jen v tomto úzkém ostrově
pravidelného pohybu na mezi stability a dává dvě předpovědi, které mohou ověřit existující minima:
excentricita A kmitá během let o 0,13–0,17 a minima A vykazují změnu okamžiků zákrytů s vnější
periodou 55–78 d. Pokud se neprojeví, je čtyřhvězda 2+2 vyloučena. Ponaučení pro program je obecné:
blízko meze stability je sekulární výsledek předvýběrem a odpověď dává N-částicový backend.

## Statistika zákrytů populací 2+2 {#populations}

**Otázka.** Jak často zakrývají oba páry čtyřhvězdy 2+2 a jak rychle se zákryty objevují a mizí?
To kalibruje, co může SQUADRA očekávat, že uvidí.

**Vstup.** 50 000 syntetických soustav 2+2 se sekulárním backendem během 3000 let. Předpokládané
apriorní rozdělení: hmotnosti rovnoměrně v rozmezí 0,8–3 M☉, vnitřní periody logaritmicky
rovnoměrně v rozmezí 0,3–30 d, vnitřní $e$ rovnoměrně v rozmezí 0–0,2, vnější periody
logaritmicky rovnoměrně od 100 d do 3000 let, vnější $e$ rovnoměrně v rozmezí 0–0,8, izotropní
orientace, stabilní podle kritéria Vynatheyi et al. (2022).

**Tabulka 33.** Statistika zákrytů syntetické populace 2+2. Řádky: podíl párů zakrývajících dnes
a podíly těch, jejichž zákryty skončí (mezi dnes zakrývajícími) nebo začnou (mezi dnes
nezakrývajícími) během 10, 100 a 1000 let. Sloupce: pár A a pár B.

| Veličina | A | B |
|---|---|---|
| zakrývá nyní | 31,7 % | 32,3 % |
| zákryty skončí během 10 / 100 / 1000 let | 3,4 / 10,2 / 19,0 % | 3,0 / 10,2 / 19,4 % |
| zákryty začnou během 10 / 100 / 1000 let | 2,6 / 12,9 / 24,7 % | 2,6 / 13,0 / 24,3 % |

![**Obrázek 11.** Syntetická populace 2+2 (50 000 soustav): podíl zakrývajících párů, jejichž zákryty skončí (vlevo), a nezakrývajících párů, jejichž zákryty začnou (vpravo), během dané doby od dneška (vodorovná osa, logaritmická, v letech). Modře: pár A; oranžově: pár B.](../../figures/showcase/population.png)

![**Obrázek 12.** Tatáž populace seřazená podle vnější periody (vodorovná osa, logaritmická). Každý bod je podíl soustav v daném intervalu period; vystínovaný pás je jeho statistická nejistota. Červeně: pár A se naklání rychleji než 0,1°/rok. Modře: pár A začne nebo přestane zakrývat během 100 let. Zeleně: oba páry zakrývají v některém okamžiku během 3000 let; tečkovaná čára je podíl očekávaný v náhodném okamžiku pro nezávislé páry.](../../figures/population_eclipses.png)

**Výsledek (tabulka 33).** Na obrázku 11 odečtěte například hodnotu při 100 letech: asi 10 %
zakrývajících párů přestane zakrývat a asi 13 % skrytých začne. Křivky A a B leží na sobě, jak
musí u párů losovaných ze stejného apriorního rozdělení, a P(B zakrývá | A zakrývá) se rovná
P(B zakrývá) – kontrola, že statistika zákrytů zachází s izotropními drahami správně (srov.
Vokrouhlický 2016, oddíl 3.3). Oba páry zakrývají společně v průměru 10 % času. Obrázek 12
ukazuje, že rychlé změny vyžadují krátké vnější periody: mezi zakrývajícími páry se náklon rychlejší
než 0,1 °/rok vyskytuje v 92 % soustav s vnějšími periodami pod 0,6 roku, ve 27 % při asi 2 letech,
v 5 % při 4 letech a nad 10 let nikdy. Rychlosti nad 1 °/rok – rychlosti změřené pro CzeV4315
a ASASSN-V J075203 – vyžadují vnější periody 105–460 d, což souhlasí s pozorovanými rychle se
naklánějícími soustavami (vnější periody 0,3–3 roku; Kolář et al. 2025; V994 Her; CzeV4315).
Pravděpodobnost, že pár zakrývá, 32 % pro toto apriorní rozdělení, souhlasí s geometrickým odhadem
$(R_1+R_2)/a$ středovaným přes apriorní rozdělení. Populace je předpokládaná, nikoli proložená;
srovnání se 781 kandidáty vyžaduje výběrovou funkci TESS, která je mimo rozsah SWINGERsimu.

## CzeV343 {#czev343}

Čtvrtou soustavou SQUADRA v části I je CzeV343, čtyřhvězda 2+2 v rezonanci period, kterou dynamicky
studovali Pejcha et al. (2022). Její příspěvky k apsidálnímu stáčení jsou čtyřmi ze srovnávacích
testů SWINGERsimu ([kapitola 3](#tides)): slapová, rotační a relativistická rychlost páru A souhlasí
s článkem na 0,6–2 % a úplná N-částicová čtyřhvězda dává pozorovanou apsidální periodu s přesností
3,3 %. CzeV343 je podrobně analyzována v části IV.
