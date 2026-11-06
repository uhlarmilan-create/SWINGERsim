# Pozorovatelné veličiny: světelné křivky, okamžiky zákrytů a pulzace {#observables}

Astronomové nepozorují dráhy; pozorují světlo. Tato kapitola ověřuje, jak SWINGERsim převádí své
dráhy a hvězdy na to, co zaznamená dalekohled: světelné křivky zákrytových hvězd (kulových,
deformovaných a kontaktních), zeslabení světla prachovými disky, radiální rychlosti, okamžiky
minim zákrytů a jejich odchylky od přísné periody a časování pulzací. Všechny moduly této kapitoly
jsou vlastním kódem SWINGERsimu s výjimkou hybridního modelu PHOEBE, který pro světlo hvězd
deformovaných podle Rocheova modelu volá externí kód PHOEBE 2.5 (Prša et al. 2016).

## Světelné křivky kulových hvězd {#light-curves}

**Definice.** Základní fotometrický model SWINGERsimu (v programu nazývaný **tier A**) pracuje
s hvězdami jako s koulemi s radiálně symetrickým okrajovým ztemněním (čtyřkoeficientový zákon
Clareta 2000) a tok zakrytý během zákrytu počítá Gaussovou–Legendreovou integrací přes viditelný
disk. Okrajové ztemnění pochází z tabulek Clareta & Bloemenové (2011, atmosféry ATLAS) a pro hvězdy
teplejší než 15 000 K z tabulek Reeva & Howartha (2016, atmosféry TLUSTY bez LTE). Elipsoidální
proměnnost – změna jasnosti slapově protažené hvězdy – se přidává jako hlavní člen podle Morrise
(1985).

**Reference.** Analytický zákryt hvězdy s okrajovým ztemněním tmavým diskem (Mandel & Agol 2002)
v implementaci batman (Kreidberg 2015); hloubky zákrytů oddělené dvojhvězdy z koulí v ellc (Maxted
2016); elipsoidální proměnnost Rocheova modelu v ellc; PHOEBE 2.5 se stejnou fyzikou (koule, černá
tělesa, Claretovo okrajové ztemnění); publikované koeficienty TLUSTY a překryv TLUSTY a ATLAS
v rozmezí 15 000–30 000 K.

**Tabulka 23.** Světelné křivky kulových hvězd a okrajové ztemnění. Sloupce jsou definovány v oddílu [Jak číst tuto zprávu](#how-to-read).

| Test | Třída | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|---|
| `lc-occultation-small` | own | 0,987891 | 0,987891 | −1,9·10⁻¹² | 10⁻⁷ | PASS |
| `lc-occultation-limb` | own | 0,945236 | 0,945236 | +3,7·10⁻⁹ | 10⁻⁷ | PASS |
| `lc-occultation-large` | own | 0,128979 | 0,128979 | −5,5·10⁻⁹ | 10⁻⁷ | PASS |
| `lc-detached-primary` | own | 0,380821 | 0,380868 | +4,7·10⁻⁵ | 10⁻⁴ | PASS |
| `lc-detached-secondary` | own | 0,144658 | 0,144654 | −3,9·10⁻⁶ | 10⁻⁴ | PASS |
| `lc-ellipsoidal` | own | 0,0027949 | 0,00260619 | −6,75 % | 10 % | PASS |
| `phoebe-sphere-vs-tier-a` | own | 0 | 0,00037 | +0,00037 | 0,001 | PASS |
| `ld-tlusty-row` | own/tab | 1,18262 | 1,18262 | +0 | 10⁻⁴ | PASS |
| `ld-tlusty-vs-atlas` | own/tab | 1 | 0,997587 | −0,24 % | 2 % | PASS |

**Výsledek (tabulka 23).** Pro malé (poměr poloměrů 0,1), tečné a velké (1,2) zakrývající těleso
souhlasí tok s analytickým modelem s přesností lepší než 10⁻⁸: integrace je pro koule prakticky
přesná. Hloubky zákrytů oddělené dvojhvězdy souhlasí s ellc na 5·10⁻⁵, což je přesnost integrační
mřížky ellc. Proti PHOEBE 2.5 se stejnou fyzikou je největší rozdíl normovaných světelných křivek
přes celý oběh ($e$ = 0,1) 3,7·10⁻⁴ toku. Elipsoidální člen pro malé deformace leží při
$R/a$ = 0,15 o 6,8 % pod Rocheovým modelem ellc; při $R/a$ = 0,10, 0,20 a 0,25 je deficit 3, 11
a 17 %. Tento test vyznačuje mez kulového modelu: těsnější páry potřebují povrchové sítě
z následujícího oddílu. Koeficienty TLUSTY jsou reprodukovány v bodech mřížky a intenzity
středované přes disk z TLUSTY a ATLAS souhlasí na 0,2–0,7 % tam, kde platí obě, takže přechod mezi
nimi při 15 000 K je hladký.

## Deformované a kontaktní hvězdy na povrchových sítích {#meshes}

**Definice.** Pro hvězdy, které mají daleko ke kulovému tvaru, sestavuje SWINGERsim **povrchové
sítě**: povrch tuze rotující zploštělé hvězdy nebo hvězdy, která vyplňuje (nebo přeplňuje) svůj
Rocheův lalok, se rozdělí na elementy 64 × 128 směrů, každý s vlastní teplotou (gravitační
ztemnění), okrajovým ztemněním a promítnutou plochou; zákryty se počítají element po elementu.
Tytéž sítě dávají **střed světla** hvězdy: v kontaktní dvojhvězdě je radiální rychlost měřená ze
spektrálních čar tokem vážená rychlost viditelného povrchu, která se liší od rychlosti těžiště
(a během zákrytu vykazuje Rossiterův–McLaughlinův efekt; Rossiter 1924; McLaughlin 1924).

**Tabulka 24.** Povrchové sítě, kontaktní dvojhvězdy a rychlosti středu světla. Sloupce jsou definovány v oddílu [Jak číst tuto zprávu](#how-to-read).

| Test | Třída | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|---|
| `mesh-sphere-eclipse` | own | 0 | 0,00026 | +0,00026 | 0,001 | PASS |
| `contact-vs-phoebe` | own | 0 mag | 0,0030 mag | +0,0030 mag | 0,03 mag | PASS |
| `rv-light-vs-phoebe` | own | 0 | 0,022 | +0,022 | 0,05 | PASS |

**Výsledek (tabulka 24).** Sítě koulí reprodukují přesné zákryty koulí z předchozího oddílu na
2,6·10⁻⁴ toku; okraje hvězd v popředí jsou vyhlazeny (anti-aliasing) podle promítnuté velikosti
každého elementu. Kontaktní dvojhvězda ($q$ = 0,42, faktor přeplnění 0,3, $i$ = 80°) má stejnou
amplitudu ve filtru V jako v PHOEBE 2.5 s černými tělesy a stejným okrajovým ztemněním s přesností
0,003 mag. Rychlost středu světla její sekundární složky se od tokem vážené rychlosti z PHOEBE liší
nejvýše o 2,2 % amplitudy 52 km/s; největší rozdíly nastávají uvnitř zákrytů, kde sítě obou kódů
rozlišují zakrývaný okraj odlišně. Automatické testy k tomu ověřují, že mimo zákryt je posun kulové
hvězdy přesně nulový, že signál Rossiterova–McLaughlinova efektu je antisymetrický kolem středu
zákrytu s očekávanou velikostí $v_\mathrm{rot}(R_\mathrm{occ}/R)^2$ a že rychlosti kontaktní
dvojhvězdy spočtené se středem světla, ale proložené se středem hmoty, dávají $\chi^2$ více než
dvakrát větší a hmotnosti o několik procent příliš velké, kdežto proložení s
`fit.rv.centre: light` vrací hmotnosti v rámci 3σ.

## Prachové disky jako zakrývající tělesa {#discs}

**Definice.** Některé zákryty nezpůsobuje hvězda, ale tmavý prachový disk kolem průvodce nebo
prstenec kolem celé dvojhvězdy. SWINGERsim reprezentuje disk jako **desku** (slab) – válec
s vnějším poloměrem $r_\mathrm{out}$, volitelným vnitřním otvorem $r_\mathrm{in}$, tloušťkou $H$
a optickou hloubkou $\tau$. Hvězda za ním je zeslabena faktorem $e^{-\tau L/H}$, kde $L$ je délka
zorného paprsku uvnitř desky. Okraj může být ostrý nebo měkký, disk může vyzařovat a rozptylovat
světlo a může být pevný, stáčet se předepsanou rychlostí nebo se stáčet jako tuhé těleso působením
momentu síly své dvojhvězdy nebo průvodce. Disk kolem excentrické dvojhvězdy se může ustálit
v **polární** poloze, kolmo k dráze dvojhvězdy; nad kritickým sklonem libruje kolem polární
orientace (Farago & Laskar 2010; Martin & Lubow 2017).

**Reference.** Přesná geometrie pro propustnost; rychlost precese testovací částice pro úzký tuhý
prstenec (Zanazzi & Lai 2018, rov. 2); $\cos I_\mathrm{crit} = \sqrt{5e^2/(1+4e^2)}$ pro polární
libraci; a tři pozorované soustavy: ε Aur (veleobr třídy F zakrývaný každých 27,1 roku asi na dva
roky diskem kolem průvodce, hvězdy třídy B; dráha podle Stefanika et al. 2010, model s malou
hmotností podle Hoarda et al. 2010, disk podle Kloppenborga et al. 2015), HD 98800 B (polární
cirkumbinární disk skloněný o 26° k rovině oblohy; Kennedy et al. 2019; dráha podle
Zúñigy-Fernándeze et al. 2026) a EE Cep (hvězda třídy B zakrývaná každých 5,6 roku s hloubkami,
které se mění s periodou 11–12 oběhů, 62–67 let; Gałan et al. 2012; Pieńkowski et al. 2020).

**Tabulka 25.** Prachové disky. Sloupce jsou definovány v oddílu [Jak číst tuto zprávu](#how-to-read).

| Test | Třída | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|---|
| `disc-strip-uniform` | own | 0,623838 | 0,62475 | +0,000912 | 0,005 | PASS |
| `disc-face-on-tau` | own | 0,496585 | 0,496585 | +0 | 10⁻⁹ | PASS |
| `disc-inclined-tau` | own | 0,246597 | 0,246597 | −4·10⁻¹⁶ | 10⁻⁹ | PASS |
| `disc-ring-precession` | own | 1 | 1,00018 | +0,02 % | 0,1 % | PASS |
| `disc-polar-critical` | own | 37,7612° | 37,7612° | +6·10⁻⁹° | 0,2° | PASS |
| `epsaur-mid-eclipse` | own | JD 2455429,5 | JD 2455414,5 | −15 d | 30 d | PASS |
| `epsaur-duration` | own | 740 d | 638 d | −13,78 % | 15 % | PASS |
| `epsaur-depth` | own | 0,75 mag | 0,83134 mag | +0,0813 mag | 0,15 mag | PASS |
| `hd98800-polar-disc` | own | 26° | 28,4735° | +2,47° | 3° | PASS |
| `eecep-disc-precession` | own | 65 let | 183,163 roku | +181,79 % | 15 % | **EXPLAINED** |

**Výsledek: geometrie a dynamika (tabulka 25).** Neprůhledná deska viděná z boku s poloviční
tloušťkou 0,3 $R$ přes hvězdu s rovnoměrnou jasností propustí 0,6248 světla proti přesné hodnotě
0,6238: disk hvězdy je vzorkován slunečnicovou mřížkou 4000 bodů, která rozliší ostrý okraj asi
na 10⁻³ toku. Disky viděné shora a skloněné s optickou hloubkou 0,7 propouštějí přesně $e^{-0{,}7}$
a $e^{-1{,}4}$. Úzký tuhý prstenec ve vzdálenosti čtyř vzdáleností složek dvojhvězdy se stáčí
rychlostí testovací částice s přesností 0,02 % a rovnice tuhého disku přesně reprodukuje kritický
sklon polární librace pro $e$ = 0,5 (37,76°).

![**Obrázek 5.** Prachové disky jako zakrývající tělesa, s neprůhledným diskem s ostrým okrajem. Ve všech panelech je na svislé ose zeslabení v magnitudách, rostoucí směrem dolů, jak je astronomové vynášejí. *Vlevo* (ε Aur): spočtená světelná křivka zákrytu v letech 2009–2011 v závislosti na juliánském datu; čárkované svislé čáry: pozorovaný první a poslední kontakt; tečkovaná vodorovná čára: pozorovaná hloubka. *Uprostřed* (EE Cep, ilustrativní disk): hloubka každého zákrytu (body) během 70 let, kdy se disk stáčí s předepsanou periodou 65 let; šedý pás je pozorovaný rozsah hloubek. *Vpravo* (KH 15D, ilustrativní prstenec): 25 let světelné křivky mladé dvojhvězdy, která mizí za okrajem stáčejícího se cirkumbinárního prstence a znovu se vynořuje.](../../figures/discs.png)

**Výsledek: pozorované soustavy.** Obrázek 5 ukazuje tři disky s ostrým okrajem. U ε Aur (vlevo)
plyne okamžik středu zákrytu ze samotné spektroskopické dráhy a připadá 15 d před pozorovaný střed
kontaktů (16. srpna 2009 až 26. srpna 2011). Zákryt trvá 638 d proti pozorovaným 740 d (−14 %)
a je hluboký 0,83 mag proti 0,75 mag. Kratší trvání způsobuje ostrý okraj: pozorovaný vstup do
zákrytu trvá asi 145 d (Stencel 2013), takže skutečný okraj je měkký. Hloubka je o 0,08 mag příliš
velká, protože model vynechává světlo hvězdy třídy B uvnitř disku a světlo rozptýlené diskem.
U HD 98800 B je disk umístěný kolmo k dráze B a podél jejího vektoru excentricity – polární
orientace předpovídaná pro disky kolem excentrických dvojhvězd – skloněn o 28,5° k rovině oblohy,
spočteno z publikované dráhy bez jakéhokoli proložení; pozorovaný disk je skloněn o 26°.
V prostředním panelu dává ilustrativní disk EE Cep (průměr 1,6 au, $\tau$ = 1, náklon 15°), jehož
precese je *předepsána* s periodou 65 let, hloubky zákrytů 0,19–1,73 mag, blízké pozorovaným
0,5–2 mag. V pravém panelu reprodukuje ilustrativní prstenec KH 15D (neprůhledný, 3–5 au, náklon
20°, precese za 600 let) pozorovaný sled: nejprve jedna hvězda dvojhvězdy vykukuje poblíž apoastra
(špičky), pak jsou obě po léta skryty a pak se dvojhvězda znovu vynoří.

![**Obrázek 6.** Disky s měkkými okraji a tuhou dynamikou, vývoj do kontaktu a podíly Pejchy et al. *Vlevo nahoře* (ε Aur): zeslabení v závislosti na juliánském datu pro neprůhledný disk s ostrým okrajem (oranžově) a pro disk s měkkým okrajem a vlastním světlem (modře); šedá plocha: pozorovaný interval od prvního do posledního kontaktu. *Vpravo nahoře* (EE Cep): hloubka každého zákrytu (magnitudy, body spojené čarami) během 65 let pro disk s vnějším poloměrem 0,8 au, jehož precese s periodou 65 let je předepsána (oranžově), a pro tuhý disk s vnějším poloměrem 1,6 au, jehož precese buzená momentem síly má periodu 65 let (modře); šedý pás: pozorované hloubky. *Vlevo dole:* hmotnosti obou hvězd (plné čáry, levá osa) a jejich vzdálenost (čárkovaně, pravá osa) dvojhvězdy 1,0 + 0,5 M☉, která se vyvíjí do kontaktu a splyne (viz [kapitola 8](#life-cycles)). *Vpravo dole:* podíl sady Pejchy et al. (2013) (uspořádání Vokrouhlického, 400 členů), jejíž dvojhvězda A (modře) nebo B (oranžově) dosáhla $1-e<10^{-3}$ po sekulárních bězích dlouhých 11 kyr a 50 kyr; tečkovaně: publikované podíly.](../../figures/part1/v027_jobs.png)

**Výsledek: měkké okraje a momenty sil.** Obrázek 6 ukazuje, co se stane, když disky dostanou
měkký okraj, vlastní světlo a precesi spočtenou z momentu síly místo předepsané. U ε Aur (vlevo
nahoře) způsobí měkký okraj, že vstup do zákrytu trvá 141 d, jak se pozoruje (asi 145 d; ostrý
okraj dává 89 d), zákryt však trvá 829 d, o 12 % déle, než se pozoruje; šířka měkkého okraje je
předpokládaná, nikoli proložená, a disk přidává jen 0,1 % světla. U EE Cep se ilustrativní disk
s vnějším poloměrem 0,8 au, učiněný tuhým, stáčí působením momentu síly hvězdy třídy B s periodou
183 let – ne 65 let. Test je vysvětlen (EXPLAINED) jako omezení disku, nikoli jako chyba: rychlost
škáluje jako $r_\mathrm{out}^{3/2}$ a 65 let vyžaduje $r_\mathrm{out}\approx$ 1,6 au (nebo
hmotnější hvězdu B). S touto velikostí a ilustrativní tloušťkou a optickou hloubkou jsou však
zákryty 3,7–14 mag hluboké a 36–158 d dlouhé (vpravo nahoře, modře), daleko mimo pozorovaných
0,5–2 mag. Disk, který je tuhý a stáčí se s pozorovanou periodou, proto musí být také tenčí nebo
méně neprůhledný než ilustrativní. Prstenec KH 15D o poloměru 3–5 au, učiněný tuhým, se za 60 let
pootočí jen o 2,6° (perioda asi 6500 let) a skrývá dvojhvězdu 98 % času; změna pozorovaná během
desetiletí potřebuje 1000–3000 let, takže prstenec musí být blíže k dvojhvězdě. To jsou omezení
disků, která SWINGERsim odvozuje proto, že jeho dynamika disků prochází výše uvedenými
analytickými testy. Dolní panely obrázku 6 patří ke [kapitole 8](#life-cycles) a k ansámblům
Pejchy et al. v [kapitole 2](#pejcha2013): publikované podíly (tečkovaně) leží mezi hodnotami po
11 kyr a po 50 kyr.

## Okamžiky zákrytů a diagramy O−C {#timings}

**Definice.** Pozorovatel předpovídá okamžiky minim zákrytové dvojhvězdy s konstantní periodou,
$T_n = T_0 + nP$, a vynáší **diagram O−C** (pozorovaný minus vypočtený okamžik v závislosti na
čísle cyklu $n$). Jakákoli struktura v něm je signálem. **Efekt světelného času** (LTTE) třetího
tělesa posouvá minima až o $a_{12}\sin i/c$, kde $a_{12}$ je velká poloosa dráhy dvojhvězdy kolem
společného těžiště (Irwin 1952). **Dynamická** změna okamžiků zákrytů (ETV) je změna periody
dvojhvězdy působením přitažlivosti třetího tělesa. **Apsidální stáčení** posouvá primární
a sekundární minima opačnými směry. Uvnitř samotné dvojhvězdy je sekundární zákryt zpožděn
**Rømerovým zpožděním** $(P/\pi c)(K_2-K_1)$, protože obě hvězdy jsou při obou konjunkcích
v různých vzdálenostech od nás (Kaplan 2010). Ustálená změna periody $\dot P$ ohne O−C do paraboly,
$\tfrac12\dot P P n^2$ (Kwee 1958; Hilditch 2001).

SWINGERsim nabízí tři způsoby výpočtu minim: přesné konjunkce z N-částicové integrace včetně doby
šíření světla (`swingersim timings`); minima změřená ze syntetické světelné křivky tak, jak by je
změřil pozorovatel; a hybrid PHOEBE, který počítá světlo hvězd deformovaných podle Rocheova modelu
pomocí PHOEBE 2.5 po úsecích N-částicové dráhy. Pro pozorovaná minima poskytuje dvě proložení:
**krok 1**, analytický model O−C, který skupina SQUADRA používá v programu OCFit (Gajdoš &
Parimucha 2019) – LTTE, apsidální stáčení v řadě Giméneze & Bastera (1995) a kvadratický člen pro
každý pár – a **krok 2**, společný N-částicový model všech zákrytových párů (ověřovaný v částech
II–IV).

**Tabulka 26.** Okamžiky zákrytů a proložení v kroku 1. Sloupce jsou definovány v oddílu [Jak číst tuto zprávu](#how-to-read).

| Test | Třída | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|---|
| `timing-ltte-amplitude` | own | 41,5837 min | 41,5863 min | +0,01 % | 0,1 % | PASS |
| `phenom-romer-delay` | own | 27,4517 s | 27,3578 s | −0,34 % | 1,5 % | PASS |
| `phoebe-romer-delay` | ext | 27,452 s | 27,287 s | −0,60 % | 2 % | PASS |
| `pdot-parabola` | R | 1 | 1,00005 | +0,00 % | 0,1 % | PASS |
| `ocfit-gimenez-bastero` | own | 0 s | 0,00063 s | +0,00063 s | 0,01 s | PASS |
| `ocfit-lite-amplitude` | own | 1 | 1 | −0,00 % | 10⁻⁴ % | PASS |
| `ocfit-nbody-asini-A` | own | 1 | 0,99864 | −0,14 % | 1 % | PASS |
| `ocfit-nbody-asini-B` | own | 1 | 1,0024 | +0,24 % | 1 % | PASS |
| `ocfit-nbody-P3` | own | 1 | 0,99996 | −0,004 % | 0,5 % | PASS |

**Výsledek (tabulka 26).** Syntetický diagram O−C dvojhvězdy na kruhové vnější dráze viděné z boku
o velikosti 10 au (stejné hmotnosti) má poloviční amplitudu LTTE podle Irwina (1952) s přesností
10⁻⁴; dynamická zpoždění této široké soustavy jsou pod 0,01 min. Ve dvojhvězdě ZAMS 2 + 1 M☉
s $P$ = 5 d a $i$ = 89° přicházejí sekundární minima změřená ze světelné křivky 27,36 s po polovině
periody, proti Rømerovu zpoždění 27,45 s (−0,34 %); hybrid PHOEBE dává 27,29 s (−0,60 %).
Dvojhvězda, jejíž perioda se mění o 4·10⁻⁶ d/rok během 60 let, integrovaná N-částicově s impulzy
periody, dává kvadratický koeficient svých minim s přesností 5·10⁻⁵; proložení 160 takových minim
v kroku 2 nalezne $\dot P$ v rámci jeho chyby (pod 5 % hodnoty). V proložení kroku 1 souhlasí
apsidální člen SWINGERsimu spočtený z přesných konjunkcí s devítičlennou řadou Giméneze & Bastera
(1995), kterou používá OCFit, na 0,6 ms při $e$ = 0,3 (při $e$ = 0,6 se useknutá řada odchyluje
o 1,6 s), a jeho člen LTTE se rovná členu z OCFit. Konečně okamžiky zákrytů N-částicové čtyřhvězdy
2+2 s vnější dráhou s periodou 1500 d ($e$ = 0,3, $i$ = 60°), proložené v kroku 1 s jednou společnou
vnější dráhou, vracejí keplerovské $a\sin i$ obou párů na 0,14 % a 0,24 % a vnější periodu
na 4·10⁻⁵. Zbývajícím rozdílem je dynamické zpoždění, řádově $(3/8\pi)P_\mathrm{in}^2/P_\mathrm{out}$,
které analytický model kroku 1 neobsahuje.

## Časování pulzací {#pulsations}

**Definice.** Pulzující hvězda jsou hodiny. Ve dvojhvězdě přicházejí její pulzace dříve a později
o dobu šíření světla napříč její dráhou, přesně jako zákryty. SWINGERsim zahrnuje pulzace do svých
světelných křivek s periodami z publikovaných vztahů (vztah perioda–svítivost hvězd typu δ Scuti,
pulzační konstanta $Q$, modely cefeid) a zpoždění měří tak jako pozorovatel: prokládá fázi pulzací
v krátkých úsecích světelné křivky (Murphy et al. 2014).

**Reference.** Pro kruhovou dráhu je amplituda časových zpoždění $a_1\sin i/c$ (Murphy et al. 2014,
oddíl 2.7, rov. 10); pro excentrickou dráhu závisí poměr druhé a první harmonické zpoždění jen na
excentricitě; KIC 9651065, hvězda typu δ Scuti na dráze s periodou 272,7 d (Murphy et al. 2014,
tabulka A1); vztah perioda–svítivost 2375 galaktických hvězd typu δ Scuti (Poro et al. 2024,
rov. 10) se základní pulzační konstantou $Q$ = 0,0329 d (Breger 1979); a vztah perioda–poloměr
rotujících ženevských modelů cefeid (Anderson et al. 2016, tabulka 5) s periodami cefeid podle
Bona et al. (2000).

**Tabulka 27.** Časování pulzací a periody pulzací. Sloupce jsou definovány v oddílu [Jak číst tuto zprávu](#how-to-read).

| Test | Třída | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|---|
| `pulsation-time-delay-circular` | own | 53,859 s | 53,859 s | +0,00 % | 0,01 % | PASS |
| `murphy2014-kic9651065-a1` | own | 165,57 s | 167,767 s | +1,33 % | 6 % | PASS |
| `murphy2014-kic9651065-a2-a1` | own | 0,231986 | 0,219651 | −0,0123 | 0,025 | PASS |
| `delta-scuti-pl-1p8` | own/tab | 0 mag | −0,0252 mag | −0,0252 mag | 0,25 mag | PASS |
| `delta-scuti-pl-2p0` | own/tab | 0 mag | −0,115 mag | −0,115 mag | 0,25 mag | PASS |
| `cepheid-period-radius-6` | own/tab | 0 dex | −0,0077 dex | −0,0077 dex | 0,05 dex | PASS |
| `cepheid-period-radius-7` | own/tab | 0 dex | +0,0039 dex | +0,0039 dex | 0,05 dex | PASS |

**Výsledek (tabulka 27).** Pulzující hvězda podobná δ Scuti (1,7 M☉, perioda pulzací 0,1 d)
s průvodcem 0,5 M☉ ($P$ = 100 d, $i$ = 60°), vzorkovaná jako data Kepleru s dlouhou kadencí
a analyzovaná ve dvoudenních úsecích, dává amplitudu zpoždění $a_1\sin i/c$ = 53,859 s s přesností
3·10⁻⁷ poté, co se středování přes úsek opraví příslušným faktorem sinc. Pro KIC 9651065
s publikovanou dráhou ($a_1\sin i$ = 0,37 au, $e$ = 0,468, $\varpi$ = 2,01 rad) je první harmonická
zpoždění 167,8 s proti naměřeným 165,57 ± 2,03 s (+1,3 %; tolerance 6 % je nejistota
$a_1\sin i$) a poměr harmonických 0,220 proti 0,232 ± 0,013, což ověřuje strukturu zpoždění na
excentrické dráze. Hvězdy MIST o hmotnostech 1,8 a 2,0 M☉ leží v místě, kde na hlavní posloupnosti
procházejí $T_\mathrm{eff}$ = 7500 K, 0,03 a 0,12 mag od vztahu perioda–svítivost hvězd δ Scuti
(tolerance 0,25 mag odráží nejistotu $Q$, až 18 %). Cefeidy MIST o hmotnostech 6 a 7 M☉ na svých
modrých smyčkách mají poloměry v rozmezí 0,008 a 0,004 dex od ženevského vztahu perioda–poloměr;
tolerance 0,05 dex je rozdíl, který Anderson et al. (2016) zjišťují mezi předpovězenými
a empirickými poloměry, a oba řetězce modelů nesdílejí žádný kód. Světelné křivky a zpoždění
skutečných pulzujících hvězd ve vícenásobných soustavách ukazují [kapitola 7](#catalogue-pulsations)
a [kapitola 8](#life-pulsators).
