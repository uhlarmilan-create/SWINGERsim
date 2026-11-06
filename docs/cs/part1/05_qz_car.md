# QZ Car proti Brožovi et al. (2022) {#qz-car}

QZ Car (HD 93206) v mlhovině Carina je čtyřhvězda 2+2 hvězd spektrální třídy O, patřících
k nejhmotnějším známým hvězdám: zákrytová dvojhvězda Ac (Ac1 25,5 + Ac2 33,2 M☉, $P$ = 6,0 d),
excentrická dvojhvězda Aa (Aa1 69,8 + Aa2 8,85 M☉, $P$ = 20,7 d) a vnější dráha s periodou
14 728 d (asi 40 let). Brož et al. (2022) proložili 40 let zákrytů, spekter a interferometrie svým
N-částicovým kódem Xitau a publikovali úplný model. Jde o nejhmotnější soustavu v našich testech
a o přísný test geometrie a fotometrie SWINGERsimu, protože velikosti drah a podíly světla čtyř
hvězd jsou známy přesně. Jejich tabulka 2 je vstupem konfigurací QZ Car ve SWINGERsimu; poloměry
a efektivní teploty jsou převzaty z jejich tabulek 2 a 3. V konfiguraci jsou čtyři hvězdy
m1 = Ac1, m2 = Ac2, m3 = Aa1 a m4 = Aa2.

## Geometrie, světlo a vývojové hmotnosti {#qzcar-benchmarks}

**Tabulka 28.** Srovnávací testy QZ Car. První tři řádky ověřují hierarchickou inicializaci
SWINGERsimu (velké poloosy z period a hmotností), další čtyři jeho fotometrii (podíly světla
v pásmu V z poloměrů a teplot), poslední tři vývojové dráhy MIST (hmotnost modelu osamocené hvězdy
s pozorovaným poloměrem a teplotou). Sloupce jsou definovány v oddílu [Jak číst tuto zprávu](#how-to-read).

| Test | Třída | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|---|
| `qzcar-broz-a1` | own | 54,032 R☉ | 53,9946 R☉ | −0,07 % | 0,2 % | PASS |
| `qzcar-broz-a2` | own | 136,132 R☉ | 136,081 R☉ | −0,04 % | 0,2 % | PASS |
| `qzcar-broz-a3` | own | 13 053 R☉ | 13 046 R☉ | −0,05 % | 0,2 % | PASS |
| `qzcar-broz-lv-ac1` | own | 0,22 | 0,220775 | +0,000775 | 0,005 | PASS |
| `qzcar-broz-lv-ac2` | own | 0,133 | 0,133856 | +0,000856 | 0,005 | PASS |
| `qzcar-broz-lv-aa1` | own | 0,635 | 0,635024 | +2,4·10⁻⁵ | 0,005 | PASS |
| `qzcar-broz-lv-aa2` | own | 0,01 | 0,0103458 | +0,000346 | 0,005 | PASS |
| `qzcar-evolutionary-mass-ac1` | own/tab | 25,5 M☉ | 27,0258 M☉ | +5,98 % | 15 % | PASS |
| `qzcar-evolutionary-mass-ac2` | own/tab | 33,2 M☉ | 26,0243 M☉ | −21,61 % | 15 % | **EXPLAINED** |
| `qzcar-evolutionary-mass-aa1` | own/tab | 69,8 M☉ | 42,5172 M☉ | −39,09 % | 15 % | **EXPLAINED** |

Tabulka 28 uvádí deset srovnávacích testů; následující odstavce ji procházejí po skupinách.

**Geometrie (úroveň B).** Tři velké poloosy vyplývají z publikovaných period a hmotností
s přesností 0,04–0,07 %; tolerance 0,2 % pokrývá zaokrouhlení hmotností v tabulce 2 na tři
platné číslice.

**Světlo (úroveň C).** Černá tělesa s okrajovým ztemněním podle Clareta/TLUSTY a s publikovanými
poloměry a teplotami dávají podíly čtyř hvězd na světle v pásmu V – 0,221, 0,134, 0,635 a 0,010 –
s přesností 0,001 vůči hodnotám, které Brož et al. odvodili ze syntetických spekter. Tolerance
0,005 je rozdíl mezi černým tělesem a modelovou atmosférou při těchto teplotách.

**Vývojové hmotnosti (úroveň D).** **Vývojová hmotnost** hvězdy je hmotnost vývojové dráhy MIST,
která dosáhne jejího pozorovaného poloměru a teploty; **dynamická hmotnost** plyne z drah.
U normální osamocené hvězdy se shodují. Ac1 je na vývojových drahách MIST hvězda 27,0 M☉ ve stáří
5,6 Myr, v rozmezí 6 % od své dynamické hmotnosti. Ac2 (26,0 M☉, 5,1 Myr) a především Aa1
(42,5 M☉, 3,9 Myr) svým dynamickým hmotnostem 33,2 a 69,8 M☉ neodpovídají. Tyto dva testy jsou
vysvětleny (EXPLAINED): obě hvězdy pravděpodobně hmotu získaly. Nominální model Brože et al.
(2022, tabulka 2) prokládá teplotu Ac2 (32 979 K) a nechává její rotaci volnou; proložení dává
410 km/s, znak hvězdy, která akreovala moment hybnosti při přenosu Ac1 → Ac2, jejž Brož et al.
popisují. U Aa1 Brož et al. poznamenávají, že je ve srovnání s Aa2 extrémně hmotná, a naznačují
dřívější přenos Aa2 → Aa1. Hvězda, která hmotu získala, je omlazena (do jádra je vmíchán čerstvý
vodík) a jeví se mladší a méně hmotná než osamocená hvězda stejné hmotnosti; vývojové dráhy
osamocených hvězd ji popsat nemohou. Umístěny na své dynamické hmotnosti dosahují vývojové dráhy
MIST pozorovaných poloměrů ve stáří 6,0 Myr (Ac1), 3,5 Myr (Ac2) a 2,5 Myr (Aa1) a rotující dráhy
nikdy nedají Aa1 její pozorovaný poloměr 28 R☉ při 70 M☉: vývojové dráhy osamocených hvězd nedokážou
učinit soustavu stejně starou. SWINGERsim proto dovoluje, aby hvězda začínala s efektivním stářím
(`stars.<name>.age_offset_myr`), takže každá složka začíná se svým pozorovaným poloměrem, a během
běhu omlazuje hvězdy, které hmotu získávají ([kapitola 3](#mass-transfer)).

## Syntetická minima {#qzcar-minima}

**Vstup a metoda.** Minima zákrytového páru Ac byla spočtena přes jednu vnější periodu (40,5 roku)
kulovým modelem ze světelné křivky v pásmu V s kadencí 10 min a s ročními pozorovacími sezónami
(4933 minim) a přes jeden rok hybridem PHOEBE, který deformuje hvězdy do jejich Rocheových tvarů
(4000 trojúhelníků na hvězdu).

![**Obrázek 7.** Diagram O−C syntetických minim zákrytového páru Ac soustavy QZ Car přes jednu vnější periodu asi 40 let. Vodorovná osa: čas; svislá osa: O−C ve dnech, okamžik každého minima minus okamžik předpovězený s konstantní periodou. Modře: primární minima; oranžově: sekundární minima. Vlna je efektem světelného času vnější dráhy.](../../figures/showcase/qzcar_oc.png)

**Výsledek.** Obrázek 7 ukazuje O−C 4933 minim. Kdyby byla Ac osamocená, ležely by všechny body na
vodorovné přímce. Místo toho sledují vlnu asi ±0,2 d (rozsah O−C 0,38 d): když se Ac na vnější
dráze od nás vzdaluje, přicházejí její minima pozdě, když se přibližuje, přicházejí brzy. Vlna není
sinusoida, protože vnější dráha je excentrická. Její amplituda se rovná $a_{12}\sin i/c$
publikované vnější dráhy. Primární i sekundární minima leží na téže křivce, takže obojí je
časováno stejně dobře. Měřeno vůči barycentrickým konjunkcím přicházejí primární minima o 16 s
pozdě a sekundární o 10 s brzy: Rømerovo zpoždění uvnitř Ac je ±16 s (Kaplan 2010) a zbytek přidává
mírně excentrická dráha Ac. Minima z PHOEBE během jednoho roku souhlasí s minimy kulového modelu
na 1 s pro sekundární minimum a leží během tří týdnů v rozmezí +6 s a −5 s od konjunkcí. Časování
hmotné, slapově deformované zákrytové dvojhvězdy je tedy dvěma nezávislými modely světelné
křivky reprodukováno na několik sekund.

## Současná dynamika {#qzcar-dynamics}

**Vstup.** N-částicový ansámbl 2000 členů s hmotnostmi a dráhovými elementy losovanými v rámci
nejistot Brože et al. (2022) a s neznámými uzly dvojhvězd losovanými izotropně, zintegrovaný
100 let dopředu a 100 let dozadu, s podmínkou, že Ac zakrývá od roku 1955.

**Výsledek.** Všech 2000 členů zůstává v obou směrech neporušených a všechny jsou stabilní podle
kritérií Mardlinga–Aarsetha a Vynatheyi et al. (2022). Za století se dráhy téměř nezmění: sklon Ac
se posouvá s mediánem 9·10⁻⁴ °/rok (rozsah 16–84 % od −0,0011 do +0,0010 °/rok), asi o 0,04° za
40 let, což Brož et al. nemohli změřit. Ac proto zakrývá po staletí; Aa zakrývá ve 24 % izotropních
orientací svého uzlu. QZ Car je stabilní laboratoří pro studie časování.

## Vývoj do prvního přetečení Rocheova laloku {#qzcar-evolution}

**Otázka.** Hvězdy QZ Car jsou mladé (2–6 Myr), hmotné hvězdy se však vyvíjejí rychle. Která hvězda
vyplní svůj Rocheův lalok jako první a kdy?

**Vstup.** Tři ansámbly po 1000 členech s vývojovým backendem a vývojovými drahami MIST,
s hmotnostmi losovanými v rámci publikovaných chyb: (a) všechny hvězdy stejně staré, 2,5 Myr, na
nerotujících drahách; (b) totéž s rotujícími drahami ($v/v_\mathrm{crit}$ = 0,4), které mají
silnější vítr; (c) efektivní stáří zvolená tak, aby každá hvězda začínala se svým pozorovaným
poloměrem (hvězdy pak nejsou stejně staré).

![**Obrázek 8.** QZ Car: doba od dneška do prvního přetečení Rocheova laloku (vodorovná osa, Myr) ve třech ansámblech po 1000 členech (svislá osa: počet členů v intervalu 0,05 Myr). Barvy ukazují, která hvězda přeteče jako první: modře Aa1 (veleobr 69,8 M☉), zeleně Aa2, oranžově Ac1. Legenda uvádí počet členů pro každou hvězdu a v horní části členy, které skončí jinak. (a) Stejně staré hvězdy ve stáří 2,5 Myr, nerotující dráhy MIST. (b) Totéž s rotujícími drahami. (c) Efektivní stáří, která dávají každé hvězdě její pozorovaný poloměr.](../../figures/part1/qzcar_overflow.png)

**Výsledek.** Obrázek 8 ukazuje tři ansámbly. (a) Se stejně starými nerotujícími hvězdami končí
všech 1000 členů přetečením Rocheova laloku Aa1, po mediánu 0,44 Myr (16–84 %: 0,30–0,53 Myr);
vnější dráha se větrem mírně rozšíří (medián z 60,5 na 62,1 au). (b) S rotujícími drahami rozšiřuje
silnější vítr dráhy: přeteče jen 338 členů (Aa1 jako první ve 244, medián 0,92 Myr; Aa2 jako první
v 94, medián 1,10 Myr), zatímco 616 dosáhne konce svých vývojových drah bez přetečení, s Aa
rozšířenou na asi 2,3 au a vnější dráhou na asi 105 au; 46 členů je odmítnuto, protože Aa1 těžší
než asi 73,3 M☉ ve stáří 2,5 Myr na mřížce rotujících drah neexistuje. (c) S pozorovanými poloměry
končí přetečením 803 členů – Aa1 jako první v 626 (medián 0,38 Myr) a Ac1 jako první ve 177 (medián
0,12 Myr) – a 197 je odmítnuto na začátku, protože Ac1 těžší než 26,4 M☉ ve svém efektivním stáří
už přeplňuje svůj Rocheův lalok na pozorované dráze s periodou 6 d. Toto poslední číslo je samo
o sobě výsledkem: pozorovaný poloměr Ac1 a skutečnost, že Ac je oddělená, omezují hmotnost Ac1
shora, což je podmínka, kterou dynamické proložení neobsahuje. V panelu (c) obě dvojhvězdy
závodí: Ac1 vyhraje jen tehdy, přijde-li přetečení velmi brzy, u členů, v nichž je Ac1 vylosována
blízko této horní meze.

Publikovaná simulace budoucnosti QZ Car, se kterou by bylo možné srovnávat, neexistuje; vývoj je
předpovědí. Jeho složky jsou ověřeny samostatně: vývojové dráhy MIST hmotných hvězd proti ženevským
modelům a interpolaci s vynecháním jednoho bodu, vítr proti Vinkovi et al. (2001)
([kapitola 3](#massive-giants)), rozšiřování drah větrem proti Jeansovu zákonu (kapitoly 3
a [7](#catalogue)). Rozptyl mezi třemi ansámbly – dvojnásobek v době do přetečení – je nejistotou
modelů velmi hmotných hvězd, nikoli dynamiky.
