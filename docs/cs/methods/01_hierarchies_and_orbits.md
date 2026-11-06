---
title: "SWINGERsim"
subtitle: "Fyzika a metody"
---

# Hierarchie a dráhy {#hierarchies-and-orbits}

Tento svazek popisuje fyziku a numerické metody SWINGERsimu, simulátoru hierarchických
vícenásobných hvězd. U každé fyzikální složky uvádí, co se modeluje, jaké jsou rovnice
a předpoklady, jakou numerickou metodou se počítá, které konfigurační klíče ji řídí, jak byla
ověřena a kde přestává platit. Samotné ověření – srovnání každého modulu s publikovanými výsledky –
je dokumentováno ve svazku *Testování SWINGERsimu, část I: Ověření proti literatuře*; zde je
u každé metody jen jednovětné shrnutí jejího testu. Jak napsat konfiguraci, spouštět příkazy
a číst výstupy, popisuje *Uživatelská příručka* (kapitoly Konfigurace, Příkazový řádek a Grafické
rozhraní).

Tato první kapitola definuje geometrii: jak se vícenásobná hvězda popisuje jako strom drah, jaké
souřadnice a dráhové elementy SWINGERsim používá, jak do popisu vstupuje pozorovatel a které pomalé
(sekulární) jevy mění dráhy v čase.

## Hierarchické vícenásobné hvězdy {#hierarchical-multiple-stars}

*Vícenásobná hvězda* je gravitačně vázaná soustava tří nebo více hvězd. Zhruba polovina hvězd
podobných Slunci a téměř všechny hmotné hvězdy mají alespoň jednoho průvodce a mnoho soustav
obsahuje tři, čtyři i více hvězd (Tokovinin 2018). Taková soustava přežije dlouho jen tehdy, je-li
*hierarchická*: lze ji rozdělit na dvojice, jejichž dráhy mají velmi rozdílné velikosti. Dvě blízké
hvězdy tvoří *dvojhvězdu*; z vzdálenosti mnohem větší, než je její rozměr, působí dvojhvězda na
ostatní tělesa téměř jako jediná hvězda se součtovou hmotností. Může proto obíhat kolem jiné hvězdy
nebo jiné dvojhvězdy po mnohem širší dráze. Každá dráha je pak téměř pevná elipsa (*keplerovská*
dráha), pomalu rušená ostatními. Bez tohoto oddělení škál si hvězdy během několika oběhů vymění
partnery a navzájem se vymrští (kapitola [Stabilita a osudy](#stability-and-fates)).

Hierarchie je **binární strom**: každý uzel je dráha, která spojuje dvě složky, a každá složka je
buď hvězda (list), nebo další dráha. SWINGERsim přijme libovolný takový strom zhruba do deseti
hvězd. Běžné stromy jsou předvolbami klíče `topology`; jakýkoli jiný strom se zadá explicitně
pomocí `topology: custom` a `hierarchy`, tedy zobrazení z názvu každé dráhy na její dvě složky.
Předvolby uvádí tabulka 1.

**Tabulka 1.** Hierarchie dostupné jako předvolby. *Předvolba* je hodnota klíče `topology`;
*Hvězdy* počet hvězd; *Dráhy* názvy, pod nimiž dráhy vystupují v konfiguraci a ve všech výstupech;
*Popis* uspořádání. Zápis „2+2“ znamená dvě dvojhvězdy obíhající kolem sebe, „3+1“ trojhvězdu
(dvojhvězdu s blízkou třetí hvězdou), kolem níž obíhá vzdálená čtvrtá hvězda; závorky označují
podsoustavu uvnitř větší soustavy.

| Předvolba | Hvězdy | Dráhy | Popis |
|---|---|---|---|
| `custom` | 2–3 | volné názvy, např. `A`, `out` | dvojhvězda nebo trojhvězda („2+1“) |
| `2+2` | 4 | `A`, `B`, `out` | dvě dvojhvězdy A = (m1, m2) a B = (m3, m4) na společné dráze |
| `3+1` | 4 | `inner`, `mid`, `out` | dvojhvězda (m1, m2), hvězda m3 kolem ní, hvězda m4 kolem trojhvězdy |
| `2+2+1`, `3+1+1` | 5 | jako výše a navíc `out5` | čtyřhvězda se vzdálenou pátou hvězdou |
| `2+(2+1)` | 5 | `A`, `B`, `Bt`, `out` | dvojhvězda A a trojhvězda tvořená dvojhvězdou B a blízkou hvězdou m5 |
| `2+2+2` | 6 | `A`, `B`, `C`, `AB`, `out` | čtyřhvězda 2+2 AB a třetí dvojhvězda C |
| `3+3` | 6 | `A`, `B`, `T1`, `T2`, `out` | dvě trojhvězdy na společné dráze |
| `custom` | až ~10 | volné názvy | libovolný strom vnořených dvojic |

![**Obrázek 1.** Hierarchie zakreslené jako stromy. Každý panel je jedna hierarchie; její název
(dvojhvězda, trojhvězda 2+1, 2+2, 3+1, 2+2+1, 2+(2+1), 2+2+2, 3+3 a vlastní strom deseti hvězd)
je v záhlaví. Oranžové kroužky dole jsou hvězdy m1, m2, …; zelené čtverce jsou dráhy označené
svými názvy (A, B, out, inner, mid, …); šedé čáry spojují každou dráhu s jejími dvěma složkami.
Výška dráhy ve stromu udává její úroveň: výše položená dráha je širší a pomalejší. Obrázek ukazuje,
že každá hierarchie, jakkoli velká, je složena ze stejného prvku – dráhy dvou složek –, a proto
všem slouží jedna sada rovnic.](../../figures/concepts/hierarchies.png)

Obrázek 1 ukazuje stromy předvoleb. Každý panel se čte zdola: dvě hvězdy se spojí do dráhy
(například A), dráha se spojí s další složkou do širší dráhy (out) a tak dále až nahoru. V panelu
2+2 jsou obě složky dráhy `out` dvojhvězdy; v panelu 3+1 je strom řetězcem, v němž každá dráha
obsahuje předchozí. Poslední panel ukazuje, že vlastní strom může oba vzory kombinovat. Protože
model pracuje se stromem, a nikdy s pevným počtem hvězd, jsou dynamika, diagnostiky, pozorovatelné
veličiny i proložení pro každou hierarchii stejné; předvolby poskytují jen názvy.

Zvlášť cennou třídou jsou dvojitě zákrytové soustavy 2+2, v nichž zakrývají obě dvojhvězdy,
protože zákryty omezují orientaci obou vnitřních drah (Zasche et al. 2019, 2023b, 2024; Kolář et
al. 2026). Mnoho příkladů v tomto svazku se jich týká, metody však na ně omezeny nejsou.

## Jacobiho (hierarchické) souřadnice {#jacobi-coordinates}

Každá dráha stromu popisuje **relativní pohyb své druhé složky kolem první**, přičemž každá složka
je zastoupena svým těžištěm. To jsou *Jacobiho* (neboli *hierarchické*) souřadnice. V soustavě 2+2
je dráha `out` pohybem těžiště B kolem těžiště A; v soustavě 3+1 je dráha `mid` pohybem m3 kolem
těžiště (m1, m2). Gravitační parametr dráhy je $\mu = G(M_1 + M_2)$, kde $M_1$ a $M_2$ jsou
celkové hmotnosti jejích dvou složek a $G$ gravitační konstanta.

**Počáteční podmínky.** Konfigurace zadává keplerovské elementy každé dráhy. SWINGERsim je převádí
na polohy a rychlosti rekurzivně od nejvnější dráhy dovnitř: nejvnější dráha umístí těžiště svých
dvou složek, potom každá vnitřní dráha umístí své dvě složky kolem vlastního těžiště, až k hvězdám.
Nakonec se celá soustava posune tak, aby její barycentrum bylo v klidu v počátku souřadnic. Zpětný
převod výsledného stavu na elementy vrátí vstup s přesností $10^{-12}$ (viz *Testování, část I*).

**Oskulační a střední elementy.** Elementy, které hlásí N-částicový backend, jsou *oskulační*:
keplerovská elipsa, po níž by se obě složky pohybovaly, kdyby v daném okamžiku všechny ostatní
hvězdy zmizely. Protože ostatní hvězdy nezmizí, oskulační elementy kmitají na časové škále oběhů,
i když se v průměru nic nemění. Sekulární backendy ([Dynamika](#dynamics)) pracují se *středními*
(přes oběh středovanými) elementy, které takové kmity neobsahují. Rozdíl je podstatný, když se jako
výchozí bod použije publikované fotodynamické řešení: jeho oskulační elementy v dané epoše dají
v jiném kódu mírně odlišnou střední periodu, která se nasčítá do driftu okamžiků zákrytů
([Prokládání](#anchoring)).

## Dráhové elementy a vztažná soustava {#orbital-elements}

Každou dráhu popisuje šest elementů uvedených v tabulce 2.

**Tabulka 2.** Dráhové elementy jedné dráhy. *Symbol* je označení v tomto svazku; *Klíč* název
v oddílu `orbits.<name>` konfigurace; *Význam* definice. Velikost dráhy se zadává buď periodou,
nebo velkou poloosou, nikdy oběma.

| Symbol | Klíč | Význam |
|---|---|---|
| $P$ nebo $a$ | `P` [d] nebo `a` [au] | perioda, nebo velká poloosa; svázány třetím Keplerovým zákonem |
| $e$ | `e` | excentricita, $0 \le e < 1$ |
| $i$ | `i` [deg] | sklon roviny dráhy k vztažné rovině, 0°–180° |
| $\Omega$ | `Omega` [deg] | délka výstupného uzlu ve vztažné rovině |
| $\omega$ | `omega` [deg] | argument periastra, měřený v rovině dráhy od uzlu |
| $M_0$ | `M0` [deg] | střední anomálie v čase $t = 0$ |

Třetí Keplerův zákon spojuje periodu a velkou poloosu,

$$ a^3 = \frac{G\,(M_1 + M_2)}{4\pi^2}\,P^2 , $$

kde $M_1 + M_2$ je celková hmotnost obou složek dráhy. Dvě Slunce s periodou jednoho roku jsou od
sebe vzdálena $2^{1/3} \approx 1{,}26$ au; s periodou 10 d asi 0,11 au, tj. 24 R☉.

![**Obrázek 2.** Dráhové elementy jedné dráhy. Šedá rovina je vztažná rovina (rovina oblohy, jsou-li
elementy zadány vůči obloze). Černá svislá osa je její normála $z$, která míří **od pozorovatele**;
černá přímka v rovině je vztažný směr (osa $x$), od něhož se měří uzel. Modrá elipsa je dráha, plnou
čarou nad vztažnou rovinou a čárkovaně pod ní. Tečkovaná přímka je uzlová přímka; oranžový oblouk
$\Omega$ je délka výstupného uzlu, zelená přímka normála dráhy skloněná o sklon $i$ od $z$
a fialový oblouk $\omega$ argument periastra od výstupného uzlu k periastru (fialový bod).
Dohromady $i$ a $\Omega$ určují orientaci roviny dráhy, $\omega$ orientaci elipsy v této rovině
a $a$, $e$ její velikost a tvar.](../../figures/concepts/orbit_elements.png)

Geometrii ukazuje obrázek 2. SWINGERsim používá pro všechny dráhy jednu **vztažnou soustavu**:
rovina $x$–$y$ je vztažná rovina, $z$ její normála a rotace z dráhové soustavy do vztažné soustavy
je $R = R_z(\Omega)\,R_x(i)\,R_z(\omega)$. Jsou-li elementy zadány vůči obloze – což je běžný případ
pozorované soustavy –, je vztažnou rovinou rovina oblohy a **zorný paprsek je osa $z$, mířící od
pozorovatele**. Důsledky této konvence, používané v celém svazku:

- radiální rychlost hvězdy je $v_z$, kladná, když se hvězda vzdaluje;
- světlo vyslané v hloubce $z$ dorazí později o $z/c$, takže pozorovaný okamžik děje
  v dynamickém čase $t$ je $t_\mathrm{obs} = t + z/c$, s barycentrem soustavy v $z = 0$
  ([Pozorovatelné veličiny](#light-travel-time));
- $i = 90°$ znamená, že dráhu vidíme z boku a zákryty jsou možné; $i = 0°$ nebo $180°$ znamená
  pohled shora.

**Rozsah sklonu 0°–180° a sklápění do 0°–90°.** V konfiguraci leží $i$ mezi 0° a 180°: dráha
s $i < 90°$ obíhá při pohledu ze strany pozorovatele proti směru hodinových ručiček, dráha
s $i > 90°$ po směru. Oba smysly oběhu jsou dynamicky podstatné, protože na nich závisí vzájemný
sklon dvou drah (další oddíl). Pozorovatel je však z fotometrie nerozliší: zákrytová dvojhvězda se
sklonem $i$ a $180° - i$ má stejnou světelnou křivku, neboť zákryty závisejí jen na $|\cos i|$.
SWINGERsim proto

- v dynamice zachovává celý rozsah 0°–180°;
- sklápí sklony do 0°–90° vždy, když je porovnává s fotometrickými měřeními
  ($i \to 180° - i$ pro $i > 90°$; [Prokládání](#inclination-data));
- v ansámblech losuje oba smysly oběhu se stejnou pravděpodobností, je-li u sklonu nastaveno
  `mirror: true` ([Ansámbly](#what-is-varied));
- uvádí *odklon* $|i - 90°|$ zákrytové dvojhvězdy, jehož rychlost změny má stejné znaménko pro oba
  smysly oběhu ([Ansámbly](#rates-today)).

Fyzikální smysl má jen *relativní* smysl oběhu dvou drah; jednu dráhu (podle konvence první
dvojhvězdu) lze vždy ponechat s $i < 90°$.

**Invariantní rovina.** Kromě elementů vůči obloze ukládá SWINGERsim sklon a uzel každé dráhy také
vůči *invariantní rovině*, rovině kolmé k zachovávajícímu se celkovému momentu hybnosti soustavy
(výstupní sloupce `inc_inv`, `Omega_inv`). V této rovině jsou formulovány sekulární teorie a v ní
se měří i periody precese uzlů uvedené v souhrnu běhu.

## Vzájemný sklon {#mutual-inclination}

**Vzájemný sklon** $i_\mathrm{mut}$ dvou drah je úhel mezi jejich rovinami, spočtený z vektorů
dráhového momentu hybnosti $\mathbf{h}_1$ a $\mathbf{h}_2$:

$$ \cos i_\mathrm{mut} = \frac{\mathbf{h}_1 \cdot \mathbf{h}_2}{|\mathbf{h}_1|\,|\mathbf{h}_2|} . $$

Na rozdíl od $i$ nezávisí na vztažné rovině a právě on řídí dynamiku. Hodnoty pod 90° znamenají, že
dráhy obíhají ve stejném smyslu (prográdně), nad 90° v opačném (retrográdně). Každá dvojice drah má
svůj vlastní vzájemný sklon; SWINGERsim uvádí všechny (sloupce `imut_<pair>`).

Pro dvě dráhy se sklony k obloze $i_1$, $i_2$ a uzly $\Omega_1$, $\Omega_2$ platí

$$ \cos i_\mathrm{mut} = \cos i_1 \cos i_2 + \sin i_1 \sin i_2 \cos(\Omega_1 - \Omega_2) . $$

Zákryty měří $i_1$ a $i_2$ (až na výše popsané sklápění), nikoli však uzly na obloze, takže
vzájemný sklon dvou zákrytových dvojhvězd je znám jen v rozmezí, které připouští neznámý rozdíl
$\Omega_1 - \Omega_2$. To je hlavní důvod, proč se pozorované soustavy studují pomocí ansámblů
([Ansámbly a inference](#ensembles-and-inference)).

## Precese uzlů {#nodal-precession}

Nejsou-li dvě dráhy koplanární, působí vzdálený průvodce na vnitřní dvojhvězdu momentem síly. Ten
způsobuje **precesi** roviny vnitřní dráhy: její vektor momentu hybnosti se pohybuje po kuželu kolem
celkového momentu hybnosti a $i_\mathrm{mut}$ zůstává téměř stálý. V nejnižším (kvadrupólovém) řádu
a pro kruhovou vnitřní dráhu je perioda precese

$$ P_\mathrm{prec} \approx \frac{4}{3}\,\frac{P_\mathrm{out}^2}{P_\mathrm{in}}\,
\frac{m_\mathrm{tot}}{m_\mathrm{pert}}\,\left(1-e_\mathrm{out}^2\right)^{3/2}
\frac{1}{\cos i_\mathrm{mut}} , $$

kde $P_\mathrm{in}$ a $P_\mathrm{out}$ jsou vnitřní a vnější perioda, $e_\mathrm{out}$ vnější
excentricita, $m_\mathrm{tot}$ celková hmotnost vnitřní dráhy a jejího rušitele a $m_\mathrm{pert}$
hmotnost rušitele (v soustavě 2+2 druhé dvojhvězdy). Odvození pro trojhvězdy pochází od
Söderhjelma (1975); případ 2+2 zpracoval Vokrouhlický (2016). SWINGERsim vypisuje tento odhad pro
každou dvojici drah (`t_prec`) vedle periody naměřené v běhu z uzlu v invariantní rovině. Odhad
slouží k orientaci; výsledkem je simulace.

![**Obrázek 3.** Precese uzlů a její pozorovatelný důsledek (ilustrativní soustava 2+2). Vlevo:
geometrie. Zelená šipka je celkový moment hybnosti; modré šipky jsou tři polohy normály dvojhvězdy
B, která se kolem něj otáčí po kuželu (čárkovaná elipsa). Vpravo: sklon B k obloze (svislá osa,
stupně) v závislosti na čase (vodorovná osa, roky). Modrá křivka kmitá s periodou precese; oranžový
pás, zhruba 83°–97°, vyznačuje sklony, při nichž tato dvojhvězda zakrývá. Zákryty jsou vidět, jen
dokud křivka prochází pásem – zde po několik desetiletí dvakrát za periodu precese –, a mezi tím
chybějí.](../../figures/concepts/precession_eclipses.png)

Obrázek 3 ukazuje, proč je precese pro pozorovatele důležitá. Kužel vlevo si zachovává vrcholový úhel
(vzájemný sklon), orientace B vůči zornému paprsku se však mění, takže sklon k obloze vpravo kmitá
mezi $90° - \theta$ a $90° + \theta$ (pro osu kužele ležící v rovině oblohy; $\theta$ je úhel
otevření kužele). Dvojhvězda zakrývá jen uvnitř úzkého oranžového pásu, jehož šířku určuje podmínka
zákrytu z kapitoly [Pozorovatelné veličiny](#eclipse-geometry). Pomalu precedující dvojhvězda proto
vykazuje zákryty, které se během let až staletí prohlubují, mělčí a mizí; takové změny jsou
pozorovány u několika soustav 2+2 (Zasche et al. 2025). V soustavě 2+2 precedují obě dvojhvězdy,
obecně různou rychlostí, takže kmitá i jejich vzájemný sklon. Ve větších hierarchiích se totéž děje
na každé úrovni, u širších drah pomaleji.

## Cykly ZLK {#zlk-cycles}

Leží-li vzájemný sklon vnitřní a vnější dráhy zhruba mezi 39,2° a 140,8°, moment síly rovinu dráhy
jen neotáčí: vyměňuje sklon za excentricitu. Vnitřní excentricita roste k maximu, zatímco
$i_\mathrm{mut}$ klesá k 39,2°, a pak se obě veličiny vracejí. To jsou **cykly von Zeipela–Lidova–
Kozaie (ZLK)** (Kozai 1962; Lidov 1962; přehled Naoz 2016). Okno plyne z podmínky
$\cos^2 i_\mathrm{mut} < 3/5$, tj. $39{,}23° < i_\mathrm{mut} < 140{,}77°$.

Pro zpočátku kruhovou vnitřní dráhu, vzdáleného rušitele a kvadrupólový řád závisí maximální
excentricita jen na počátečním vzájemném sklonu $i_0$,

$$ e_\mathrm{max} = \sqrt{1 - \tfrac{5}{3}\cos^2 i_0} , $$

a cykly se opakují na časové škále řádu

$$ t_\mathrm{ZLK} \sim \frac{P_\mathrm{out}^2}{P_\mathrm{in}}\,
\frac{m_\mathrm{tot}}{m_\mathrm{pert}}\,\left(1-e_\mathrm{out}^2\right)^{3/2} , $$

se symboly z předchozího oddílu. Přesná perioda závisí na počátečních podmínkách (Antognini 2015).
SWINGERsim uvádí pro každou dvojici drah, zda začíná v okně ZLK, kvadrupólové $e_\mathrm{max}$,
$t_\mathrm{ZLK}$ a z dvojnásobně středovaných rovnic Naoze et al. (2013), platných pro libovolný
poměr hmotností, také předpovězený rozsah $e$ a $i_\mathrm{mut}$ a periodu cyklu
([Dynamika](#secular-prediction)).

![**Obrázek 4.** Cykly ZLK vnitřní dvojhvězdy spočtené sekulárním backendem. Nahoře: excentricita
vnitřní dvojhvězdy v závislosti na čase v letech. Dole: vzájemný sklon ve stupních na téže časové
ose. Modře: počáteční vzájemný sklon 70°; šedě: 20°. Červená tečkovaná čára vyznačuje 39,2°, dolní
okraj okna ZLK. Při 70° vystoupá excentricita z 0,05 na zhruba 0,86, zatímco sklon klesne k čáře
39,2°, a obě veličiny se pak vrátí – sklon se vyměňuje za excentricitu. Při 20°, mimo okno, zůstávají
obě téměř konstantní.](../../figures/concepts/zlk.png)

Obrázek 4 ukazuje tuto výměnu. Oba panely sdílejí časovou osu a maxima excentricity modré křivky
se shodují s minimy sklonu; minimum sklonu se dotýká čáry 39,2°, jak předpovídá kvadrupólová teorie.
Maximální excentricita souhlasí s $\sqrt{1 - (5/3)\cos^2 70°} = 0{,}81$ až na vliv malé počáteční
excentricity a konečného momentu hybnosti vnější dráhy. Cykly ZLK jsou důležité, protože při
maximální excentricitě se hvězdy v periastru k sobě přibližují: slapy pak mohou dráhu zmenšit
(Fabrycky & Tremaine 2007), hvězdy mohou vyplnit své Rocheovy laloky nebo se srazit.

**Oktupólový řád.** Liší-li se hmotnosti vnitřní dvojhvězdy a je-li vnější dráha excentrická,
naruší další (oktupólový) člen rozvoje pravidelnost cyklů a může vyhnat excentricitu k hodnotám
blízkým 1 a dráhu překlopit (Naoz 2016). Jeho síla se měří veličinou

$$ \varepsilon_\mathrm{oct} = \frac{|m_a - m_b|}{m_a + m_b}\,\frac{a_\mathrm{in}}{a_\mathrm{out}}\,
\frac{e_\mathrm{out}}{1 - e_\mathrm{out}^2} , $$

kde $m_a$, $m_b$ jsou hmotnosti obou složek vnitřní dráhy a $a_\mathrm{in}$, $a_\mathrm{out}$
vnitřní a vnější velká poloosa. SWINGERsim uvádí $\varepsilon_\mathrm{oct}$ pro každou dvojici;
N-částicový backend obsahuje všechny řády automaticky.

**Potlačení silami krátkého dosahu.** Jakékoli další apsidální stáčení vnitřní dráhy, které je
rychlejší než cyklus ZLK – relativita, slapová a rotační deformace hvězd, blízká třetí hvězda –,
cykly potlačí, protože periastrum se pootočí dřív, než moment síly stihne excentricitu vybudovat
(Liu, Muñoz & Lai 2015; Naoz 2016). Souhrn běhu dvojici označí, je-li relativistická apsidální
perioda kratší než $t_\mathrm{ZLK}$. Slapová a rotační deformace bývá silnější: v modelu V994 Her
podle Vokrouhlického (2016) dosáhne dvojhvězda B s hmotnými body $e = 0{,}63$, se započtenou
deformací hvězd však zůstává blízko počáteční hodnoty $e = 0{,}13$ (viz *Testování, část I*). Model
čistě hmotných bodů proto u těsných dvojhvězd nadhodnocuje růst excentricity; pro ně má být zapnut
klíč `physics.tides` ([Slapy a rotace](#tides-and-rotation)).

## Apsidální stáčení {#apsidal-motion}

Periastrum excentrické dráhy se ve své rovině pomalu otáčí: **apsidální stáčení**, měřené jako
rychlost $\dot\omega$ nebo jako apsidální perioda $U = 2\pi/\dot\omega$. U zákrytové dvojhvězdy je
přímo pozorovatelné, protože primární a sekundární minima se v diagramu O−C posouvají opačnými směry
([Pozorovatelné veličiny](#apsidal-motion-in-o-c)). Přispívají tři příčiny:

- **Obecná relativita.** V prvním postnewtonovském řádu (1PN)
  $$ \dot\omega_\mathrm{GR} = \frac{3\,G\,(m_1+m_2)\,n}{c^2\,a\,(1-e^2)} , \qquad n = \frac{2\pi}{P}, $$
  kde $m_1$, $m_2$ jsou hmotnosti, $a$ velká poloosa, $e$ excentricita, $n$ střední pohyb a $c$
  rychlost světla.
- **Slapová a rotační deformace hvězd** (Sterne 1939), která závisí na apsidálních konstantách
  $k_2$ a na poloměrech hvězd ([Slapy a rotace](#classical-apsidal-motion)).
- **Ostatní hvězdy** soustavy, jejichž sekulární moment síly periastrum rovněž stáčí.

Pro každou dvojhvězdu ze dvou hvězd udává SWINGERsim klasické očekávání (Sterne plus relativita)
a rychlost naměřenou v běhu z délky periastra $\varpi = \Omega + \omega$; rozdíl je příspěvek
ostatních hvězd. V kompaktní čtyřhvězdě může být tento příspěvek stejně velký jako klasický.
Příkaz `swingersim apsidal` rozloží rychlost na slapovou, rotační, relativistickou a dynamickou
část a udá faktor, o který by se $k_2$ musela změnit, aby souhlasila s naměřenou rychlostí.

## Rezonance {#resonances}

Pro soustavy 2+2 jsou typické dva druhy téměř souměřitelnosti:

- **Rezonance středních pohybů (period).** Jsou-li periody A a B blízké poměru malých celých čísel,
  mohou být dvojhvězdy zachyceny v rezonanci. Případ 1:1 zpracovali Breiter & Vokrouhlický (2018);
  rezonance 3:2 je pozorována v CzeV343 (Pejcha et al. 2022) a její statistický nadbytek mezi
  soustavami 2+2 diskutují Kolář et al. (2025). SWINGERsim uvádí nejbližší poměr malých celých
  čísel obou period. Rezonance středních pohybů v sekulárních backendech ze své podstaty chybějí
  a sleduje je jen N-částicový backend.
- **Precesní (sekulární) rezonance.** Precedují-li A a B podobnou rychlostí, sčítají se jejich
  vzájemné momenty sil koherentně a mohou vybudit mnohem větší kmity excentricity než v kterékoli
  dvojhvězdě samotné (Pejcha et al. 2013; Hamers & Lai 2017).

Parametr $\beta$ podle Hamerse & Laie (2017, jejich rov. 32–33) porovnává rychlost precese dráhy
rušitele s rychlostí ZLK vnitřní dráhy. Pro $\beta \ll 1$ platí klasický obraz ZLK; pro
$\beta \gg 1$ cítí vnitřní dráha rušitele středovaného přes jeho precesi; pro $\beta \sim 0{,}3$–3
může být vnitřní excentricita vybuzena silně a chaoticky, a to i při počátečních sklonech hluboko
pod oknem ZLK. SWINGERsim vypisuje $\beta$ pro A a B (2+2) nebo pro vnitřní dráhu (3+1).

## Vazba rovin drah A a B {#plane-coupling}

Kdyby veškerý moment hybnosti nesla vnější dráha, precedovaly by roviny A a B kolem ní nezávisle
a úhel $J_{AB}$ mezi nimi by se měnil od $|J_A - J_B|$ do $J_A + J_B$, kde $J_A$ a $J_B$ jsou jejich
sklony k vnější dráze. Ve skutečnosti precese i vnější dráha a dvojhvězdy jsou skrze ni svázány.
Pro kruhové dráhy A a B v kvadrupólovém řádu našel Vokrouhlický (2016, odd. 2.3) dva integrály
pohybu, které omezují $J_{AB}$ maximem $J_{AB}^\mathrm{max}$, často mnohem menším než $J_A + J_B$:
v jeho modelu V994 Her je $J_A = 11{,}6°$ a $J_B = 15{,}8°$, ale $J_{AB}^\mathrm{max} = 5{,}2°$.
Roviny dvou zákrytových dvojhvězd tak mohou zůstat dlouho téměř rovnoběžné, což zvyšuje
pravděpodobnost, že zakrývají obě. SWINGERsim vypisuje $J_{AB}^\mathrm{max}$ pro každé jádro 2+2,
u něhož je tato mez informativní (`J_AB_max_analytic_deg`). Mez předpokládá kruhové dvojhvězdy;
excentrické ji mohou mírně překročit: pro pozorované excentricity V994 Her (obr. 3 článku) nachází
Vokrouhlický $J_{AB}$ mezi 1° a 6° a N-částicová integrace téhož modelu ve SWINGERsimu dosahuje
6,3° (viz *Testování, část I*).

## Větší hierarchie a pátá hvězda {#larger-hierarchies}

Skutečné čtyřhvězdy bývají součástí pětihvězd a neviděná pátá hvězda může dynamiku ovlivnit oběma
směry:

- **Destabilizace.** Pátá hvězda na široké skloněné dráze pohání cykly ZLK dráhy obou dvojhvězd;
  excentricita této dráhy roste a její periastrum se k nim přibližuje, což může čtyřhvězdu, která
  je sama o sobě stabilní, posunout přes mez stability. Zpráva o stabilitě pro tento účel uvádí
  vnořenou trojici dráhy `out` uvnitř `out5`.
- **Stabilizace.** Jakákoli precese rychlejší než cyklus ZLK jej potlačí (Liu, Muñoz & Lai 2015).
  Blízká třetí hvězda u jedné dvojhvězdy způsobí rychlou precesi této dvojhvězdy a odstíní ji;
  obdobně kvadrupól obou dvojhvězd způsobuje precesi jejich společné dráhy a omezuje cyklus ZLK
  poháněný pátou hvězdou. V testovacím modelu se skutečnými dvojhvězdami (1 + 1 M☉, $a$ = 1 au)
  místo hmotných bodů dosáhne excentricita dráhy `out` hodnoty 0,73 místo 0,80 pro hmotné body.
- **Pozorovatelné důsledky.** Průvodce jen jedné dvojhvězdy mění jen sklon a okamžiky zákrytů této
  dvojhvězdy, jak bylo navrženo pro některé dvojitě zákrytové soustavy, jejichž dvojhvězdy se
  chovají odlišně (Zasche et al. 2025).

Každá metoda tohoto svazku pracuje s libovolným stromem. SWINGERsim přitom zůstává kódem pro
*hierarchické* soustavy: nehierarchickou (demokratickou) konfiguraci sleduje N-částicový backend,
dokud se nerozpadne, jako v rozptylových experimentech ([Stabilita a osudy](#scattering)), není však
cílem diagnostik.

## Implementace {#implementation-orbits}

Vše v této kapitole je vlastní kód SWINGERsimu; kromě NumPy a SciPy pro lineární algebru
a integraci obyčejných diferenciálních rovnic se neuplatňuje žádný externí balík. Vlastní
implementace: strom hierarchie a rekurzivní převod elementů všech drah na polohy a rychlosti
a zpět (Keplerova rovnice řešená Newtonovou iterací na strojovou přesnost); vzájemné sklony
a invariantní rovina; analytické diagnostiky každé dvojice drah (kvadrupólové $e_\mathrm{max}$,
$t_\mathrm{ZLK}$, $t_\mathrm{prec}$, $\varepsilon_\mathrm{oct}$, $\varepsilon_\mathrm{SA}$, parametr
$\beta$ podle Hamerse & Laie 2017, poměry period) naprogramované podle citovaných vzorců; integrace
dvojnásobně středovaných kvadrupólových rovnic Naoze et al. (2013, rov. A27, A29, A36) integrátorem
DOP853 ze SciPy; a mez $J_{AB}^\mathrm{max}$ naprogramovaná podle Vokrouhlického (2016, rov. 26–36).
