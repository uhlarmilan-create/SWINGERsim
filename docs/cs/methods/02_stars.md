# Hvězdy {#stars}

Dynamika vícenásobné hvězdy potřebuje od každé hvězdy víc než jen hmotnost. Její poloměr rozhoduje
o srážkách, přetečení Rocheova laloku a zákrytech; její vnitřní stavba prostřednictvím apsidální
konstanty $k_2$ a momentu setrvačnosti rozhoduje o tom, jak silně ji deformují slapy a rotace; její
zářivý výkon a efektivní teplota určují její světlo. To vše se s vývojem hvězdy mění. Tato kapitola
popisuje, odkud SWINGERsim tyto veličiny bere, jak se vyvíjejí a jak se velikost hvězdy srovnává
s jejím Rocheovým lalokem. Ztrátou a přenosem hmoty, které mění samotné hmotnosti, se zabývá kapitola
[Přenos hmoty a vývoj](#mass-transfer-and-evolution).

## Hertzsprungův–Russellův diagram {#hr-diagram}

*Hertzsprungův–Russellův (HR) diagram* zobrazuje zářivý výkon $L$ hvězdy (energii, kterou vyzáří za
sekundu, v jednotkách zářivého výkonu Slunce L☉) v závislosti na její efektivní teplotě
$T_\mathrm{eff}$ (teplotě černého tělesa se stejným povrchovým tokem; horké hvězdy vlevo). Obě
veličiny souvisejí s poloměrem $R$ vztahem $L = 4\pi R^2 \sigma T_\mathrm{eff}^4$, kde $\sigma$ je
Stefanova–Boltzmannova konstanta. Většina hvězd leží na **hlavní posloupnosti**, kde v jádře spalují
vodík: čím hmotnější, tím teplejší, jasnější a kratčeji žijící. Před ní se mladá hvězda smršťuje ve
**fázi před hlavní posloupností**; po ní se jádro smršťuje a obálka nabobtná v **obra**, později
v hvězdu **asymptotické větve obrů (AGB)**. Hmotné hvězdy ztrácejí vodíkové obálky silným větrem
a stávají se horkými **Wolfovými–Rayetovými** hvězdami.

![**Obrázek 5.** Kde v HR diagramu mohou ležet členové modelované soustavy. Vodorovná osa: efektivní
teplota v kelvinech, klesající doprava (logaritmicky); svislá osa: zářivý výkon v L☉
(logaritmicky). Černé čáry: vývojové dráhy MIST hvězd o hmotnostech 0,5, 1, 2, 5, 15 a 60 M☉
(popsané u svých bodů na hlavní posloupnosti nulového stáří), od hlavní posloupnosti do konce
dráhy; šedé čárkované čáry: tytéž hvězdy ve fázi před hlavní posloupností; červeně: fáze
Wolfovy–Rayetovy hvězdy u hvězdy 60 M☉. Oranžový pás je klasický pás nestability, v němž hvězdy
pulzují. Rámečky pojmenovávají jednotlivé oblasti. Každá hvězda modelované soustavy může ležet
kdekoli na těchto drahách, takže jedna soustava může kombinovat například mladou dvojhvězdu
s obrem.](../../figures/concepts/hr_zoo.png)

Obrázek 5 shrnuje hvězdnou stránku modelu. Každá černá dráha ukazuje, jak se jedna hvězda v diagramu
s věkem pohybuje: na hlavní posloupnosti téměř stojí, pak se jako obr posouvá doprava a nahoru.
Poloměr, který je pro dráhy nejdůležitější, vzroste podél větve obrů a AGB stokrát i více. Pás
nestability vyznačuje oblast, kde z hvězdy dělají hodiny pulzace
([Pozorovatelné veličiny](#pulsating-stars)).

## Modely hvězd {#stellar-models}

*Model hvězdy* přiřazuje hmotnosti hvězdy (a u vyvíjejících se modelů i jejímu věku a metalicitě)
poloměr, zářivý výkon, efektivní teplotu, apsidální konstantu $k_2$ a konstantu setrvačnosti
$r_g^2 = I/(MR^2)$ ($I$ je moment setrvačnosti). Volí se klíčem `stellar_model.kind`; možnosti
uvádí tabulka 3.

**Tabulka 3.** Modely hvězd. *Druh* je hodnota klíče `stellar_model.kind`; *Zdroj* model nebo data,
na nichž je založen; *Rozsah* pokryté hmotnosti, metalicity a fáze; *Stavba* odkud pocházejí $k_2$
a $r_g^2$. Hodnoty zadané pro jednotlivou hvězdu v oddílu `stars.<name>` (`R`, `Teff`,
`apsidal_k2`, `rg2`, `lag_s`) mají vždy přednost před modelem.

| Druh | Zdroj | Rozsah | Stavba ($k_2$, $r_g^2$) |
|---|---|---|---|
| `zams` (výchozí) | poloměr na hlavní posloupnosti nulového stáří podle Touta et al. (1996) | 0,1–100 M☉, sluneční metalicita, bez stárnutí | polytropa, $n = 1{,}5$ pod 0,35 M☉, $n = 3$ nad |
| `sse` | vzorce vývoje osamocených hvězd podle Hurleyho, Polse & Touta (2000) | 0,1–8 M☉, $Z$ = 0,004–0,03; hlavní posloupnost, Hertzsprungova mezera, větev obrů | Claret (2023) pro 0,2–8 M☉, mimo ně polytropa |
| `mist` | vývojové dráhy MIST v1.2 (Choi et al. 2016; Dotter 2016) | 0,1–300 M☉, [Fe/H] −1,5 až +0,5, $v/v_\mathrm{crit}$ = 0 nebo 0,4; od fáze před hlavní posloupností do konce AGB | Claret (2023) pro 0,2–8 M☉, Claret (2004) nad 8 M☉ |
| `table` | tabulka CSV uživatele (hmotnost, poloměr, $k_2$, $r_g^2$, volitelně zpoždění) | podle tabulky, lineárně v hmotnosti | z tabulky |
| `python` | funkce uživatele vracející tytéž veličiny | libovolný | z funkce |

S `structure: mesa:<grid>` pocházejí strukturní konstanty místo toho z mřížky vlastních modelů MESA
uživatele (Jermyn et al. 2023) (oddíl [Vlastní mřížky MESA](#own-mesa-grids)).

**Poloměr, teplota a zářivý výkon jednotlivé hvězdy.** Je-li zadán jen poloměr hvězdy, zachová se
zářivý výkon modelu a efektivní teplota plyne ze vztahu $L = 4\pi R^2 \sigma T_\mathrm{eff}^4$; jsou-li
zadány $R$ i $T_\mathrm{eff}$, určují zářivý výkon. Klíč `radius_inflation` zvětší hvězdu vůči modelu
při stálém zářivém výkonu. Magneticky aktivní hvězdy malých hmotností jsou větší, než předpovídají
standardní modely – CM Dra o 3–4 % při sluneční metalicitě a o 5–7 % při nižší metalicitě (Feiden &
Chaboyer 2014) –, a MIST vykazuje stejný deficit.

### Výchozí volba: hvězdy hlavní posloupnosti nulového stáří {#zams-model}

Model `zams` přiřadí každé hvězdě poloměr na hlavní posloupnosti nulového stáří (ZAMS) podle
analytické aproximace Touta et al. (1996) pro sluneční metalicitu (1 M☉ → 0,89 R☉, 3 M☉ → 1,98 R☉).
Hvězdy nestárnou. Strukturní konstanty jsou konstanty polytropy – hvězdy, v níž jsou tlak a hustota
svázány vztahem $P \propto \rho^{1 + 1/n}$ – s indexem $n = 1{,}5$ pro plně konvektivní hvězdy pod
0,35 M☉ a $n = 3$ nad touto hmotností. SWINGERsim integruje Laneovu–Emdenovu rovnici polytropy spolu
s Radauovou rovnicí pro její slapovou deformaci a získá $k_2$ a $r_g^2$; reprodukuje analytické
hodnoty pro $n = 1$ a tabelované hodnoty $k_2 = 0{,}1433$ ($n = 1{,}5$) a 0,0144 ($n = 3$) podle Brookera
& Olleho (1955) (viz *Testování, část I*).

Polytropy jsou hrubou aproximací. Skutečné hvězdy hlavní posloupnosti nad zhruba 1,5 M☉ jsou
více koncentrované ke středu, s $k_2$ řádu několika $10^{-3}$, a $k_2$ s vývojem hvězdy dále klesá;
polytropa $n = 3$ nadhodnocuje $k_2$ hvězd 2–3 M☉ asi trojnásobně. Pro jakékoli kvantitativní
apsidální stáčení má stavba pocházet z modelů hvězd (`sse`, `mist`, tabulka nebo hodnoty pro
jednotlivé hvězdy).

### Vyvíjející se hvězdy: SSE {#sse-model}

S `kind: sse` mají hvězdy věk `age_myr` (počítaný od ZAMS; všechny hvězdy jsou stejně staré, pokud
některou neposune `stars.<name>.age_offset_myr`) a metalicitu `z` (v této kalibraci je sluneční
$Z = 0{,}02$). Poloměr, zářivý výkon, efektivní teplota a konvektivní obálka plynou ze vzorců Hurleyho,
Polse & Touta (2000), aproximujících podrobné modely Polse et al. (1998), pro hlavní posloupnost,
Hertzsprungovu mezeru a první větev obrů. Hmotnost a hloubka konvektivní obálky pocházejí
z rozšíření pro dvojhvězdný kód podle Hurleyho, Touta & Polse (2002), které je kalibrováno pro
$Z = 0{,}02$.

Vzorce se za běhu nevyhodnocují. Byly jednou vyhodnoceny publikovanými rutinami SSE na mřížce
320 hmotností a 5 metalicit; každá fáze je vzorkována ve stejných zlomcích fáze (*ekvivalentní
vývojové body*, jako v MIST) a hvězda libovolné hmotnosti a metalicity se interpoluje při stálém
zlomku fáze, lineárně v $\log M$ a $\log Z$. Vůči přímým výpočtům SSE je interpolace přesná na 0,8 %
v poloměru na hlavní posloupnosti, 1,7 % v Hertzsprungově mezeře a 4 % na větvi obrů (95. percentily
pro 0,8–3 M☉, $Z$ = 0,01–0,03) a na 0,2 % v hranicích fází; u hmotných obrů s nízkou metalicitou
dosahují chyby 30 %. Za zapálením helia (hoření helia v jádře, AGB, pozůstatky) se hvězdy SSE
nesledují; běh tam skončí. Pro Slunce (1 M☉, $Z$ = 0,02, 4,57 Gyr) dává model $R$ = 0,985 R☉,
$L$ = 0,958 L☉ a bázi konvektivní obálky v 0,734 $R$ oproti helioseismické hodnotě 0,713 $R$
(Christensen-Dalsgaard et al. 1991).

### Vyvíjející se hvězdy: MIST {#mist-model}

S `kind: mist` používá SWINGERsim jako nezávislou rodinu modelů vývojové dráhy MESA z MIST v1.2
(Choi et al. 2016; Dotter 2016). Pokrývají 0,1–300 M☉, [Fe/H] od −1,5 do +0,5 v krocích po 0,25
a obě počáteční rotace MIST, $v/v_\mathrm{crit}$ = 0 a 0,4 (výchozí, standard MIST). Každá dráha je
uložena v *ekvivalentních vývojových bodech* (EEP) 1–1409: fáze před hlavní posloupností (1–201),
hlavní posloupnost (202–453), větev podobrů a rudých obrů (454–604), heliový záblesk a hoření helia
v jádře (605–707), raná AGB (708–807) a tepelně pulzující AGB až po opuštění AGB (808–1409); hmotné
hvězdy končí hořením uhlíku v jádře. Hvězdy po AGB a bílí trpaslíci se nesledují. Hvězda se
interpoluje při stálém EEP mezi sousedními hmotnostmi (lineárně v $\log M$) a metalicitami, stejně
jako v izochronách MIST, a ve věku mezi body EEP. Při dvojnásobném kroku mřížky reprodukuje
interpolace poloměr na konci hlavní posloupnosti na 1 % a věky na 1–2 %. Věky se počítají od
začátku modelů MIST před hlavní posloupností, nikoli od ZAMS jako u SSE.

Dráhy MIST neobsahují stavbu konvektivní obálky, kterou potřebuje slapová disipace
([Slapy a rotace](#dissipation-from-structure)). SWINGERsim ji odvozuje: obři na větvi obrů a na AGB
a hvězdy hořící helium v jádře chladnější než 5500 K mají konvektivní obálku obsahující veškerou
hmotu nad heliovým jádrem, sahající do hloubky celé hvězdy (Hurley et al. 2002), která se na modré
smyčce mezi 5500 a 7000 K lineárně zmenšuje k nule; hvězdy hlavní posloupnosti a Hertzsprungovy
mezery přebírají podíly hmotnosti a hloubky obálky od hvězdy SSE téže hmotnosti ve stejném zlomku
její fáze; hvězdy před hlavní posloupností o hmotnostech 0,1–2,5 M☉ jsou na Hayashiho dráze plně
konvektivní a od minima své efektivní teploty (kde vzniká zářivé jádro) se obálka lineárně v EEP
zmenšuje na svou velikost na ZAMS. Pro Slunce dává dráha s [Fe/H] = 0 ve 4,57 Gyr $R$ = 1,025 R☉
a $L$ = 1,11 L☉, dráha SSE 0,985 a 0,958: obě rodiny Slunce uzavírají mezi sebe.

**Systematické rozdíly mezi rodinami.** Ve srovnání s MIST jsou doby setrvání na hlavní posloupnosti
podle SSE asi o 10 % delší a poloměry na konci hlavní posloupnosti až o 10 % menší. To je
systematická nejistota každého věku odvozeného pomocí jedné rodiny.

**Hmotné hvězdy.** Nad 8 M☉ jsou pro dynamiku důležité tři vlastnosti MIST. (1) Větry jsou silné
a závisejí na rotaci: hvězda 100 M☉ ([Fe/H] = 0) opouští hlavní posloupnost s 68 M☉ bez rotace
a s 25 M☉ při $v/v_\mathrm{crit}$ = 0,4, takže `vvcrit` je volbou historie hmotnosti
([Přenos hmoty a vývoj](#winds)). (2) Doby hoření vodíku souhlasí se ženevskými modely Ekströmové et al.
(2012) na 4–12 %, hmotnosti na konci hlavní posloupnosti se však liší až o 30 %; Choi et al. (2016)
označují 1D modely obálek velmi hmotných hvězd, v nichž dominuje záření, za nejisté. Hmotnosti
vyvinutých hmotných hvězd jsou proto volbou modelu, nikoli předpovědí. (3) Od pozdní hlavní
posloupnosti se obálky hvězd 40–70 M☉ nafukují a smršťují v rozmezí několika M☉ počáteční
hmotnosti; zde mřížka poloměr nerozliší (chyby 20–40 %) a tam, kde na poloměru záleží, je třeba
použít hmotnosti z mřížky.

**Současné hmotnosti.** Hmotnosti z dráhového řešení jsou současné hmotnosti, kdežto vývojové
dráhy jsou označeny počátečními hmotnostmi. S `physics.mass_loss: {source: track, masses: current}`
najde SWINGERsim počáteční hmotnost každé hvězdy, jejíž dráha má ve věku `age_myr` danou hmotnost.
Hmotné hvězdy mohou být nejednoznačné (ve 2,5 Myr je hvězda 26 M☉ buď hvězdou hlavní posloupnosti
s 26,4 M☉, nebo hvězdou 257 M☉ obnaženou na Wolfovu–Rayetovu hvězdu); standardně se s varováním
použije nejméně vyvinuté řešení a hmotnost, které v daném věku nedosáhne žádná dráha, se odmítne
jako rozpor.

## Strukturní konstanty: $k_2$ a $r_g^2$ {#structure-constants}

**Apsidální konstanta** $k_2$ měří, jak silně hvězda reaguje na slapovou nebo rotační poruchu; je
rovna polovině potenciálového *Loveova čísla* $k_{2,\mathrm{Love}} = 2k_2$ používaného ve slapové
teorii. Hvězda koncentrovaná ke středu (hmota soustředěná v jádře) má malé $k_2$: polytropa $n = 3$
má 0,0144, homogenní těleso 0,75. **Konstanta setrvačnosti** $r_g^2 = I/(MR^2)$ určuje, kolik momentu
hybnosti rotující hvězda uchovává. Obě veličiny vstupují do slapových sil, vývoje rotace
i apsidálního stáčení ([Slapy a rotace](#tides-and-rotation)).

U `sse` a `mist` pocházejí obě z modelů MESA podle Clareta (2023): 69 vývojových drah, 0,2–8 M☉,
[Fe/H] = −0,5, 0, +0,5, od fáze před hlavní posloupností po bílého trpaslíka. Podle Claretova
doporučení se tabulky interpolují v povrchovém tíhovém zrychlení $\log g$ hvězdy – které nese
vývojový stav nezávisle na kódu vývojové dráhy – podél monotónního úseku každé dráhy (fáze před
hlavní posloupností, hlavní posloupnost a větev obrů, hoření helia v jádře), pak lineárně v $\log M$
a [Fe/H]; [Fe/H] hvězdy SSE je $\log(Z/0{,}02)$. Hodnoty reprodukují řádky tabulek s publikovanou
přesností. Na hlavní posloupnosti $k_2$ klesá, jak se hvězda stává koncentrovanější ke středu (pro
2 M☉ z 0,005 na 0,003); na větvi obrů je hluboká konvektivní obálka zvýší o více než řád. Nad 8 M☉
(MIST) pocházejí konstanty od Clareta (2004) (0,8–125 M☉, $Z$ = 0,02, od ZAMS po vyčerpání uhlíku),
interpolované stejným způsobem a souhlasící s Claretem (2023) při 8 M☉ na 1 %; nad 125 M☉ se použije
dráha 125 M☉. Hvězdy AGB přebírají hodnoty na vrcholu větve obrů, protože Claretovy mřížky modely
AGB neobsahují.

**Rotace.** Rotující hvězda je koncentrovanější ke středu než statická. S `physics.k2_rotation`
opraví SWINGERsim modelové konstanty,
$\delta\log k_2 = -0{,}87\,\lambda_s$ s $\lambda_s = 2\Omega^2R^3/(3GM) = 2V^2/(3gR)$, kde $\Omega$ je
úhlová rychlost rotace, $V$ rovníková rotační rychlost, $g$ povrchové tíhové zrychlení, $R$ poloměr
a $M$ hmotnost (Claret 1999; Claret & Giménez 2010 tuto opravu používají u pozorovaných dvojhvězd),
s rotací každé hvězdy v daném okamžiku.

### Vlastní mřížky MESA {#own-mesa-grids}

Pro strukturní veličiny, které publikované mřížky neposkytují, umí SWINGERsim řídit místní instalaci
MESA (Paxton et al. 2011; Jermyn et al. 2023). Soubor mřížky uvádí hmotnosti, metalicity a fyziku
(směšovací délku, překmit konvekce, koncový bod); každá dráha začíná ve fázi před hlavní posloupností
a zapisuje profily stavby. Redukce spočte z každého profilu bázi konvektivní obálky, konvektivní
obrátkový čas, $k_2$ z Radauovy rovnice (metodou Clareta 2023), $I/(MR^2)$ a koeficient disipace
inerciálních vln ([Slapy a rotace](#inertial-waves)). Na polytropách zapsaných ve formátu MESA
reprodukuje redukce $k_2$ na 1 %, $I/(MR^2)$ na 0,3 % a koeficient inerciálních vln plně konvektivní
hvězdy na 1 %.

## Rocheovy laloky {#roche-lobes}

Ve dvojhvězdě, která rotuje spolu se svou dráhou, má efektivní potenciál (gravitace obou hvězd plus
odstředivý člen) mezi hvězdami sedlový bod, vnitřní Lagrangeův bod L1. Ekvipotenciální plocha
procházející L1 uzavírá dvě oblasti tvaru kapky, **Rocheovy laloky**: hmota uvnitř laloku patří jeho
hvězdě. Hvězda, která nabobtná za svůj lalok, přelévá plyn přes L1 na svého průvodce.

![**Obrázek 6.** Rocheovy laloky kruhové synchronizované dvojhvězdy s poměrem hmotností $q = 0{,}5$,
řez v rovině dráhy. Černé body jsou středy hmotnější (vlevo) a méně hmotné (vpravo) hvězdy. Šedé
čárkované čáry jsou ekvipotenciály efektivního potenciálu; silná červená čárkovaná osmička je
kritická ekvipotenciála procházející vnitřním Lagrangeovým bodem L1 (červený křížek), která ohraničuje
oba Rocheovy laloky. Velké křivky u okrajů jsou ekvipotenciály daleko od dvojhvězdy, kde převládá
odstředivý člen. Hvězda, jejíž povrch dosáhne laloku, předává hmotu průvodci přes
L1.](../../figures/concepts/roche.png)

Geometrii ukazuje obrázek 6. Blízko každé hvězdy jsou ekvipotenciály téměř kulové, takže hvězda
mnohem menší než její lalok je téměř kulová. Blízko laloku se plochy směrem k L1 zašpičatují: hvězda,
která svůj lalok téměř vyplňuje, je viditelně deformovaná, což je důležité pro její světelnou křivku
([Pozorovatelné veličiny](#distorted-stars)). Lalok hmotnější hvězdy je větší.

SWINGERsim používá *objemově ekvivalentní* poloměr laloku $R_L$ (poloměr koule se stejným objemem
jako lalok) podle Eggletona (1983),

$$\frac{R_L}{d} = \frac{0{,}49\,q^{2/3}}{0{,}6\,q^{2/3} + \ln(1 + q^{1/3})}, \qquad
q = \frac{M}{M_\mathrm{companion}},$$

kde $d$ je vzdálenost hvězd, $M$ hmotnost hvězdy a $M_\mathrm{companion}$ hmotnost jejího průvodce;
pro stejné hmotnosti je $R_L/d \approx 0{,}38$. Pro každou dvojhvězdu ze dvou hvězd se *faktor
vyplnění* $R/R_L$ vyhodnocuje v periastru, $d = a(1-e)$ oskulační dráhy, jak je u excentrických drah
obvyklé (Sepinsky et al. 2007). Hvězda, která obíhá kolem podsoustavy – obří třetí hvězda kolem
vnitřní dvojhvězdy –, má lalok vůči celkové hmotnosti této podsoustavy (Eggleton 1983). Ukládá se
největší faktor vyplnění každé hvězdy; s `run.stop_on_roche_overflow` se běh zastaví při
$R/R_L = 1$ a s `physics.mass_transfer` se místo toho sleduje přenos
([Přenos hmoty a vývoj](#roche-lobe-overflow)).

**Meze.** Vzorec předpokládá kruhovou synchronizovanou dvojhvězdu; na excentrické nebo asynchronní
dráze se lalok podél dráhy mění a jeho hodnota v periastru je přibližná (Sepinsky et al. 2007).
Pro tvary deformovaných hvězd ve světelných křivkách řeší SWINGERsim místo toho úplný Rocheův
potenciál ([Pozorovatelné veličiny](#distorted-stars)).

## Druhy soustav {#kinds-of-systems}

Tentýž aparát popisuje soustavy ve velmi odlišných životních stadiích. Tabulka 4 je vodítkem k fyzice,
kterou každý druh potřebuje; kapitoly uvedené v posledním sloupci ji popisují.

**Tabulka 4.** Druhy vícenásobných hvězd a fyzika, kterou potřebují. *Druh* popisuje soustavu;
*Model* doporučený `stellar_model.kind`; *Fyzika* přepínače v oddílu `physics` (a klíče jednotlivých
hvězd), na nichž záleží; *Backend* integrační metodu z kapitoly [Dynamika](#dynamics); *Kapitola*
místo, kde je fyzika popsána.

| Druh | Model | Fyzika | Backend | Kapitola |
|---|---|---|---|---|
| mladé, dosud v plynu | `mist`, věk < 1 Myr | `gas`, `tides`, `tidal_lag: hurley2002` | sekulární | Dynamika, Slapy |
| podobné Slunci, hlavní posloupnost | `zams` nebo `mist` | `gr`, `tides`, `rotation`; `magnetic_braking` | kterýkoli | Slapy, Přenos hmoty |
| těsné zákrytové dvojhvězdy ve vícenásobné hvězdě | `mist` nebo `sse` | `gr`, `tides`, `rotation` | N-částicový pro okamžiky minim | Slapy, Pozorovatelné veličiny |
| obři, AGB | `mist` | `evolution`, `mass_loss`, `tides`, `tidal_lag: hurley2002` | sekulární | Přenos hmoty |
| algolidy, obří třetí hvězdy | `sse` nebo `mist` | `evolution`, `mass_transfer` | sekulární | Přenos hmoty |
| kontaktní dvojhvězdy, rotující hvězdy | `mist` nebo `sse` | `mass_transfer.contact`; `shape`, `contact_fillout` | sekulární pro vývoj, N-částicový pro světlo | Přenos hmoty, Pozorovatelné veličiny |
| pulzující složky | `mist` | `stars.<m>.pulsation`, `evolution` | N-částicový pro světlo | Pozorovatelné veličiny |
| hmotné hvězdy, Wolfovy–Rayetovy | `mist`, `vvcrit: 0.4` | `evolution`, `mass_loss`, `tides` | sekulární | Přenos hmoty |
| velmi široké (> 1000 au) | libovolný | `galactic_tide`, `encounters` | kterýkoli | Dynamika |
| kompaktní pozůstatky | hmotné body | `gr`, `gw` | kterýkoli | Dynamika |

Užitečný zvyk převzatý z dokumentace programu PHOEBE: zapínat efekty po jednom a porovnávat běhy.
Každý přepínač ukáže, co daná fyzikální složka se soustavou dělá, a výsledek, který se výrazně mění
s efektem pokládaným za zanedbatelný, si zaslouží bližší pohled.

## Implementace {#implementation-stars}

**Externí data.** Vzorce SSE byly jednou, mimo běh, vyhodnoceny publikovanými fortranskými rutinami
SSE Hurleyho, Polse & Touta (2000) do uložené mřížky; vývojové dráhy MIST v1.2 (Choi et al. 2016;
Dotter 2016) jsou převedeny z publikovaných archivů EEP; strukturní konstanty jsou tabulky Clareta
(2023) a Clareta (2004); MESA (Jermyn et al. 2023) je externí program, který uživatel spouští pro
vlastní mřížky. **Vlastní kód.** Aproximace poloměru na ZAMS podle Touta et al. (1996) a polytropy
(Laneova–Emdenova a Radauova rovnice integrované pomocí SciPy); interpolace v EEP a ve zlomku fáze
pro obě rodiny vývojových drah; interpolace Claretových tabulek v $\log g$ podél monotónních úseků;
konvektivní obálky hvězd MIST (obři, modrá smyčka, fáze před hlavní posloupností) odvozené výše
popsaným způsobem; hledání počátečních hmotností ze současných hmotností; rotační oprava $k_2$;
redukce profilů MESA ($k_2$ z Radauovy rovnice, $I/(MR^2)$, obrátkový čas, koeficient inerciálních
vln); a Eggletonovy (1983) laloky vyhodnocované v periastru pro dvojhvězdy i pro hvězdy obíhající
kolem podsoustav.
