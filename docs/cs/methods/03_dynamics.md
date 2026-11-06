# Dynamika {#dynamics}

SWINGERsim počítá pohyb hvězd dvěma vzájemně se doplňujícími způsoby, nazývanými *backendy*.
**N-částicový backend** integruje pohybové rovnice všech hvězd přímo; je přesný v rámci zahrnuté
fyziky, vidí každou dráhu a každé těsné přiblížení a je pomalý, když jsou vnitřní dráhy krátké a běh
dlouhý. **Sekulární backendy** středují každou dráhu přes její oběh a sledují jen pomalé změny
středovaných elips; jsou tisíckrát rychlejší a dosáhnou miliard let, předpokládají však, že
hierarchie zůstane neporušená. Backend se volí klíčem `run.backend` (`nbody` nebo `secular`). Tato
kapitola popisuje oba, diagnostiky počítané během běhu a slabé vnější poruchy (plyn, slapy Galaxie,
procházející hvězdy, gravitační vlny).

## Přímá N-částicová integrace {#nbody}

Všechny hvězdy (2 až zhruba 10) se integrují jako hmotné body – při zapnutých slapech jako rozlehlá
tělesa ([Slapy a rotace](#tides-and-rotation)) – integrátorem IAS15 z REBOUND (Rein & Liu 2012;
Rein & Spiegel 2015). IAS15 je Gaussův–Radauův integrátor 15. řádu s adaptivním krokem. Jeho chyba
zůstává na úrovni zaokrouhlovací chyby aritmetiky s pohyblivou řádovou čárkou, asi $10^{-15}$ na
krok, a bez zvláštního nastavení zvládá těsná přiblížení i vysoké excentricity. Cenou je rychlost:
krok se řídí nejkratší dráhou, takže výpočetní čas roste s $t_\mathrm{max}/P_\mathrm{in}$, délkou
běhu dělenou nejkratší periodou.

Protože se řeší úplná N-částicová úloha, jsou automaticky zahrnuty všechny řády interakce
(kvadrupól, oktupól i vyšší), krátkoperiodické členy, které středování přes oběh odstraní,
rezonance středních pohybů, vazba mezi oběma dvojhvězdami soustavy 2+2 i jakýkoli rozpad hierarchie.

## Relativita {#relativity}

S `physics.gr: true` (výchozí) přidá SWINGERsim pro všechna tělesa pohybové rovnice prvního
postnewtonovského řádu (1PN) pomocí efektu `gr_full` z REBOUNDx (Tamayo et al. 2020). Jednodušší
efekty REBOUNDx `gr` a `gr_potential` předpokládají jednu dominantní centrální hmotu, což ve
vícenásobné hvězdě neplatí, a nepoužívají se. Hlavním důsledkem je stáčení periastra
$\dot\omega_\mathrm{GR}$ každé dráhy ([Hierarchie a dráhy](#apsidal-motion)); simulované stáčení
izolované dvojhvězdy (1,5 M☉, $a$ = 0,05 au, $e$ = 0,3) souhlasí s analytickým vzorcem na 1 % (viz
*Testování, část I*).

Pro většinu hvězdných dvojhvězd je relativita malou opravou, s apsidálními periodami tisíců až
desetitisíců let při periodách řádu dní. Rozhodující se stává u dvojhvězd s periodami kolem jednoho
dne nebo kratšími a s pomalu působícím vzdáleným průvodcem, kde může cykly ZLK úplně potlačit
([Hierarchie a dráhy](#zlk-cycles)).

## Kontroly zákonů zachování {#energy-error}

Relativní chyba energie $|E/E_0 - 1|$ se zaznamenává v každém výstupu. Bez relativity zůstává blízko
$10^{-14}$. S relativitou je kontrolovanou veličinou zachovávající se 1PN hamiltonián všech těles,
který již obsahuje newtonovskou kinetickou a potenciální energii; zůstává na $10^{-13}$–$10^{-11}$,
i při vysokých excentricitách a se slapy. Se slapy jsou zahrnuty rotační a slapové energie. Celkový
moment hybnosti se kontroluje v newtonovském tvaru, který rovnice 1PN přesně nezachovávají:
s relativitou kolísá na úrovni asi $10^{-7}$, aniž by rostl. Běh se zastaví s událostí
`energy_error`, když chyba překročí `run.max_energy_error` (výchozí $10^{-5}$), což je vysoko nad
oběma úrovněmi. Když se energie ztrácí fyzikálně – slapovou disipací, ztrátou hmoty, gravitačními
vlnami –, energie se vypisuje, ale jako kritérium zastavení se již nepoužívá; se slapovou disipací
se místo ní kontroluje zachování celkového momentu hybnosti drah a rotací (`angmom_error`).

## Srážky {#collisions}

Každá hvězda má poloměr, z modelu hvězdy nebo z konfigurace. S `run.collisions` (standardně
zapnuto) testuje REBOUND v každém kroku, zda se povrchy dvou hvězd dotýkají; běh se pak zastaví
s událostí `collision`, časem a zúčastněnými hvězdami – nebo s `physics.mergers` pokračuje se
splynulou hvězdou ([Přenos hmoty a vývoj](#mergers)). Výchozí poloměry modelu `zams` jsou poloměry
hvězd hlavní posloupnosti nulového stáří; pro každou soustavu, v níž na srážkách záleží, je třeba
použít poloměry vyvinutých hvězd.

## Rozpoznání hierarchie {#hierarchy-detection}

V každém výstupu SWINGERsim určí, které hvězdy jsou vázány na které. Vychází z toho, že každá hvězda
je samostatný objekt, a opakovaně spojuje dvojici, jejíž relativní dvoučásticová dráha je vázaná
a má nejmenší velkou poloosu, nahradí tuto dvojici jejím těžištěm a skončí, když nezbývá žádná
vázaná dvojice. Stabilní soustava 2+2 dá jeden strom, `((m1, m2), (m3, m4))`. Nová struktura musí
přetrvat `run.hierarchy_persistence` po sobě jdoucích výstupů (výchozí 3), než je potvrzena, takže
jediný neobvyklý snímek poblíž těsného průchodu se nepočítá. Potvrzené změny se klasifikují jako
výměna, přestavba, vymrštění, rozpad nebo zachycení ([Stabilita a osudy](#fates)). Rozpoznání je
čistě energetické (vázáno, nebo ne); neposuzuje dlouhodobou stabilitu, kterou ukáže teprve další
dynamika.

## Diagnostiky dvojic drah {#pair-diagnostics}

Pro každou dráhu, která má ve stromu nadřazenou dráhu (A a B v `out` u soustavy 2+2; `inner` v `mid`
a `mid` v `out` u soustavy 3+1), porovnává SWINGERsim analytická očekávání s během. Uvádí je
tabulka 5.

**Tabulka 5.** Diagnostiky dvojice drah (vnitřní dráhy a dráhy kolem ní). *Veličina* je to, co se
vypisuje; *Základ* vzorec nebo teorie; *Význam* jak ji číst. Všechny analytické hodnoty se počítají
z počátečních elementů a jsou odhady pro orientaci; výsledkem jsou simulované hodnoty.

| Veličina | Základ | Význam |
|---|---|---|
| okno ZLK, $e_\mathrm{max}$, $t_\mathrm{ZLK}$ | kvadrupólová teorie testovací částice ([Hierarchie a dráhy](#zlk-cycles)) | zda a jak silné cykly ZLK lze očekávat |
| rozsah $e$, $i_\mathrm{mut}$, perioda cyklu | dvojnásobně středované kvadrupólové rovnice, libovolný poměr hmotností (Naoz et al. 2013) | předpověď sekulárního cyklu |
| $t_\mathrm{prec}$ | kvadrupólová precese uzlů ([Hierarchie a dráhy](#nodal-precession)) | očekávaná perioda precese |
| relativistická apsidální perioda, příznak potlačení ZLK | stáčení periastra v řádu 1PN | relativita potlačí ZLK, je-li rychlejší |
| $\varepsilon_\mathrm{oct}$ | Naoz (2016) | síla oktupólového členu |
| $\varepsilon_\mathrm{SA}$ | Luo, Katz & Dong (2016) | platnost středování přes oběh; varování nad 0,01 |
| poměry podle Mardlinga–Aarsetha a Vynatheyi et al. | kritéria stability ([Stabilita a osudy](#criteria)) | > 1 stabilní |
| $\beta$ | Hamers & Lai (2017) | vazba precese a ZLK |
| $J_{AB}^\mathrm{max}$, poměr period | Vokrouhlický (2016); nejbližší poměr malých celých čísel | vazba rovin v 2+2, rezonance |
| naměřený rozsah $e$, rozsah $i_\mathrm{mut}$, počet a perioda cyklů ZLK, naměřená perioda precese | běh | simulovaný výsledek |

### Sekulární předpověď v kvadrupólovém řádu {#secular-prediction}

Pro každou dvojici drah integruje SWINGERsim přes oběh středované kvadrupólové rovnice Naoze et al.
(2013, jejich rov. A27, A29 a A36), s rušitelem pojatým jako hmotný bod. Na rozdíl od vzorce pro
testovací částici platí tyto rovnice pro libovolný poměr hmotností, přesně zachovávají celkový moment
hybnosti a dovolují vnější dráze měnit sklon. Předpovídají rozsah $e$ a $i_\mathrm{mut}$ a periodu
cyklu excentricity. Implementace reprodukuje publikovanou křivku Naoze et al. (2013, obr. 4) na 0,04°
a Kozaiovu mez na $5\times10^{-5}$ (viz *Testování, část I*).

### Extrémy mezi výstupy {#extremes}

Maxima excentricity v cyklech ZLK mohou být mnohem kratší než rozestup výstupů; v excentrickém testu
ZLK klesne $1-e$ na několik staletí až na zhruba $10^{-5}$. N-částicový backend proto interně vzorkuje
stav každých `run.extremes_every_orbits` period nejkratší vnitřní dráhy (výchozí 10) a uchovává pro
každou dráhu největší excentricitu a nejmenší vzdálenost periastra $q = a(1-e)$ s jejich časy a pro
každou dvojhvězdu ze dvou hvězd největší faktor vyplnění Rocheova laloku v periastru. Rozestup
několika oběhů stačí, protože tato maxima jsou sekulární a trvají mnoho vnitřních oběhů; srážky se
zachycují zvlášť v každém kroku integrátoru. Navýšení výpočetní náročnosti činí několik procent.

### Indikátor chaosu MEGNO {#megno}

S `run.megno: true` počítá SWINGERsim ze variačních rovnic REBOUND veličinu MEGNO (Mean Exponential
Growth factor of Nearby Orbits; Cincotta & Simó 2000). MEGNO směřuje k 2 pro pravidelný
(kvaziperiodický) pohyb a pro chaotický pohyb roste lineárně s časem; její směrnice udává
Ljapunovův čas, dobu, za niž vzroste vzdálenost sousedních trajektorií $e$-krát. Chaos znamená, že
jednotlivý běh nelze předpovídat za hranici několika Ljapunovových časů a smysl mají jen statistiky
přes mnoho běhů ([Ansámbly a inference](#ensembles-and-inference)). Variační rovnice REBOUND
pokrývají jen newtonovskou gravitaci, takže MEGNO se počítá v paralelní newtonovské kopii běhu
s hmotnými body, bez relativity, slapů a ztráty hmoty (MEGNO spolu se ztrátou hmoty se odmítne);
výpočetní náročnost se zhruba zdvojnásobí.

## Integrace zpět v čase {#backward}

Záporné `run.t_max_yr` integruje do minulosti, s oběma backendy. Gravitace i relativistické členy
jsou vratné v čase: integrace běhu dopředu a jeho koncového stavu zpět vrátí počáteční elementy
s přesností $10^{-12}$. Pro minulost platí tři omezení:

- **Chaos.** Za hranicí několika Ljapunovových časů narostou nejistoty pozorovaných elementů tak, že
  jediná minulost nemá smysl; rozdělení možných minulostí pak dá jen ansámbl běhů v rámci nejistot.
- **Slapová disipace** běžící zpět dělá dráhy širšími a excentričtějšími; výsledek přímo závisí na
  špatně známé rychlosti disipace.
- **Hvězdný vývoj** běží také zpět: hvězdy se zmenšují směrem k ZAMS (SSE) nebo k začátku fáze před
  hlavní posloupností (MIST), kde se běh zastaví (`stellar_birth`).

Přežití v obou směrech času je podmínkou používanou v ansámblech ([Stabilita a osudy](#survival)).

## Sekulární backendy {#secular}

Středování přes oběh nahradí každou keplerovskou dráhu prstencem hmoty rozložené podél ní úměrně
době, kterou tam hvězda stráví. Interakční energie středovaná přes vnitřní i vnější dráhu
(*dvojnásobné středování*) závisí jen na tvarech a orientacích drah, nikoli na polohách hvězd, takže
velké poloosy jsou (v úloze hmotných bodů) konstantní a rovnice obsahují jen pomalé, sekulární
časové škály. Krok integrátoru se pak řídí těmito časovými škálami místo oběžných period.

### Vektorové elementy a pohybové rovnice {#vector-elements}

Každá dráha X je popsána dvěma vektory, normovaným momentem hybnosti
$\mathbf{K}_X = \sqrt{1-e_X^2}\,\mathbf{k}_X$ a Laplaceovým vektorem
$\boldsymbol{\pi}_X = e_X\,\mathbf{p}_X$, kde $\mathbf{k}_X$ je jednotková normála dráhy
a $\mathbf{p}_X$ jednotkový vektor k periastru (Vokrouhlický 2016). Tyto elementy nemají singularity
při nulové excentricitě ani sklonu. Jejich pohybové rovnice (jeho rov. 12–17) jsou

$$\frac{d\mathbf{K}_X}{dt} = -\frac{1}{\Lambda_X}\left(\boldsymbol{\pi}_X\times\frac{\partial H}{\partial\boldsymbol{\pi}_X} + \mathbf{K}_X\times\frac{\partial H}{\partial\mathbf{K}_X}\right), \qquad
\frac{d\boldsymbol{\pi}_X}{dt} = -\frac{1}{\Lambda_X}\left(\mathbf{K}_X\times\frac{\partial H}{\partial\boldsymbol{\pi}_X} + \boldsymbol{\pi}_X\times\frac{\partial H}{\partial\mathbf{K}_X}\right),$$

kde $\Lambda_X = \mu_X n_X a_X^2$ je moment hybnosti kruhové dráhy ($\mu_X$ její redukovaná hmotnost,
$n_X$ její střední pohyb, $a_X$ její velká poloosa) a $H$ středovaná interakční energie celé
soustavy. Rovnice se integrují Rungeovou–Kuttovou metodou osmého řádu DOP853 ze SciPy s relativní
tolerancí `run.secular_rtol` (výchozí $10^{-10}$).

### Čtyřhvězdy 2+2 {#secular-2p2}

Pro soustavu 2+2 je $H$ součtem kvadrupólových a oktupólových členů dvou trojhvězd: dvojhvězdy A
rušené těžištěm B a dvojhvězdy B rušené těžištěm A (Vokrouhlický 2016, rov. 18–25). Oktupólové členy
lze vypnout pomocí `run.secular_octupole: false`. K $H$ jsou přidány členy, které závisejí jen na
excentricitě každé dráhy a reprodukují relativistické stáčení apsid (s `physics.gr`) a stáčení apsid
ze slapové a rotační deformace hvězd podle Sterna (1939) (s `physics.tides`). Protože všechny tyto
členy jsou konzervativní, zůstává $H$ přesně zachováno a slouží jako kontrola integrace. Se slapovou
disipací se velké poloosy, excentricity a rychlosti rotace A a B navíc vyvíjejí podle teorie
konstantního časového zpoždění ([Slapy a rotace](#tidal-evolution)). Běh z obr. 2 Vokrouhlického
(2016) přes 50 tisíc let trvá asi 5 s místo zhruba hodiny s N-částicovým backendem.

### Soustavy 3+1 a další hierarchie {#secular-nested}

Pro každou jinou hierarchii (3+1 a stromy s pěti a více hvězdami) používá SWINGERsim obecnou
formulaci vnořených dvojic podle Hamerse et al. (2015) a Hamerse & Portegiese Zwarta (2016),
zapsanou ve stejných vektorech a pohybových rovnicích. $H$ je součtem členů *dvojice drah*
(binary-pair), jednoho pro každou dráhu $k$ a každou dráhu $l$, která ji obsahuje:

$$H_{kl} = -\frac{G\,\mu_k M_o a_k^2}{8\,a_l^3}\,f_2(\mathbf{K}_k,\boldsymbol{\pi}_k,\mathbf{K}_l)
- \frac{15}{64}\,\frac{G\,\mu_k M_o a_k^3}{a_l^4}\,\frac{m_{k1}-m_{k2}}{M_k}\,s\,f_3(\mathbf{K}_k,\boldsymbol{\pi}_k,\mathbf{K}_l,\boldsymbol{\pi}_l),$$

kde $\mu_k$ je redukovaná hmotnost dráhy $k$, $a_k$ a $a_l$ velké poloosy, $M_o$ hmotnost složky
dráhy $l$, která neobsahuje $k$, $m_{k1}$, $m_{k2}$ hmotnosti obou složek $k$ a $M_k$ jejich součet,
$s = \pm 1$ udává, ve které složce $l$ dráha $k$ leží, a $f_2$, $f_3$ jsou kvadrupólová a oktupólová
funkce vektorových elementů (Hamers et al. 2015). Pro čtyřhvězdu 2+2 jsou to přesně členy backendu
2+2 a oba backendy dávají stejný výsledek. Nejsou zahrnuty členy *trojice drah* (binary-triplet;
oktupólový řád, $(a_k/a_m)^2(a_m/a_l)^3$ pro $k$ uvnitř $m$ uvnitř $l$) ani členy hexadekapólového
a vyššího řádu; srovnání s kódem MSE (Hamers et al. 2021), který je obsahuje, ukazuje pro testované
soustavy rozdíly v $e_\mathrm{max}$ pod 0,001. Relativita, deformace hvězd a slapová disipace
v každé dráze dvou osamocených hvězd se přidávají stejně jako u backendu 2+2; pro soustavu 2+2 se
oba backendy shodují v rámci integrační tolerance a rychlost cirkularizace je rovna rychlosti podle
Huta (1981).

### Výstup sekulárních backendů {#secular-output}

Výstup má stejnou strukturu jako výstup N-částicového backendu. Elementy jsou střední elementy;
střední anomálie se nesleduje a ukládá se jako 0. Hierarchie se nemůže změnit a srážky nemohou
nastat. Rotace hvězd zůstávají souosé s jejich drahami; u rotujících hvězd to způsobuje malou chybu
v kontrole momentu hybnosti, v typickém příkladu řádu $10^{-5}$.

### Kdy středování přes oběh selhává {#secular-validity}

Středování přes oběh předpokládá, že vnější dráha je mnohem pomalejší než vnitřní. Velikost
zanedbaných členů měří parametr jednoduchého středování (Luo, Katz & Dong 2016)

$$ \varepsilon_\mathrm{SA} = \left(\frac{a_\mathrm{in}}{a_\mathrm{out}(1-e_\mathrm{out}^2)}\right)^{3/2}
\frac{m_\mathrm{pert}}{\sqrt{m_\mathrm{in}\,m_\mathrm{tot}}} , $$

kde $a_\mathrm{in}$, $a_\mathrm{out}$ jsou velké poloosy, $e_\mathrm{out}$ vnější excentricita,
$m_\mathrm{in}$ hmotnost vnitřní dráhy, $m_\mathrm{pert}$ hmotnost rušitele a $m_\mathrm{tot}$ jejich
součet. Tabulka 6 ukazuje tento efekt pro jednu dvojhvězdu a rušitele při rozšiřování vnější dráhy.

**Tabulka 6.** Perioda cyklu ZLK téže dvojhvězdy a rušitele spočtená N-částicovým backendem
a dvojnásobně středovanou teorií pro tři vnější velké poloosy. *$a_\mathrm{out}$* je uvedena vůči
nejtěsnějšímu případu, který má $\varepsilon_\mathrm{SA}$ = 0,03; *rozdíl* je (N-částicová −
středovaná)/středovaná. Rozdíl rychle klesá, jak hierarchie sílí.

| $a_\mathrm{out}$ | N-částicová perioda | středovaná perioda | rozdíl |
|---|---|---|---|
| 1 × | 12,8 let | 18,9 let | −32 % |
| 2 × | 153,8 let | 151,4 let | +1,5 % |
| 4 × | 1211 let | 1213 let | −0,2 % |

Při $\varepsilon_\mathrm{SA}$ = 0,03 rozsah excentricity ještě souhlasí, perioda cyklu je však
o třetinu kratší; při dvojnásobné vzdálenosti ($\varepsilon_\mathrm{SA}$ osmkrát menší) je rozdíl na
úrovni procent. Souhrn běhu proto varuje, je-li $\varepsilon_\mathrm{SA} > 0{,}01$. Blízko velmi
vysokých excentricit mohou krátkoperiodické členy zanedbané středováním vyhnat oskulační
excentricitu výše než střední (Luo, Katz & Dong 2016); při reprodukci populační studie Pejchy et al.
(2013) našel N-částicový backend třikrát více trojhvězd dosahujících $1-e < 10^{-3}$ než sekulární,
kdežto pro čtyřhvězdy se oba shodly.

**Shrnutí platnosti.** Sekulární backendy jsou sekulární teorie prvního řádu. Jsou spolehlivé pro
silné hierarchie ($\varepsilon_\mathrm{SA}$ hluboko pod 0,01) daleko od rezonancí středních pohybů.
Zanedbávají krátkoperiodické členy a v soustavě 2+2 přímou vazbu A a B řádu $\varepsilon^4$;
nedokážou popsat těsná přiblížení, výměny, vymrštění ani srážky. V chaotických oblastech sledují
N-částicové řešení jen po určitou dobu: pro model z obr. 4 Vokrouhlického (2016) se oba shodují asi
4 tisíce let a pak se rozejdou, stejně jako v článku.

![**Obrázek 7.** Tatáž trojhvězda spočtená oběma backendy. Vodorovná osa: čas v letech; svislá osa:
excentricita vnitřní dvojhvězdy. Tenká šedá křivka je N-částicové řešení, které ukazuje každý
krátkoperiodický kmit oskulační excentricity; silná modrá křivka je sekulární (přes oběh středované)
řešení. Obě ukazují tytéž cykly ZLK s periodou asi 70 let a maximem kolem 0,77; sekulární křivka
sleduje průměr N-částicové křivky.](../../figures/concepts/nbody_vs_secular.png)

Obrázek 7 ukazuje vztah obou backendů pro dobře oddělenou trojhvězdu. Šedá křivka kmitá na časové
škále vnější dráhy kolem modré, ale žádná z nich neujíždí: v tomto režimu reprodukuje sekulární
backend pomalý vývoj za zlomek výpočetní ceny. Tabulka 7 shrnuje, kdy použít který backend.

**Tabulka 7.** Volba backendu. *Situace* je otázka nebo soustava; *Backend* doporučený backend;
*Důvod* fyzika, na níž volba stojí.

| Situace | Backend | Důvod |
|---|---|---|
| stabilita, rozpad, výměny, srážky | N-částicový | jen on sleduje nehierarchický pohyb |
| dráhy nedostatečně oddělené ($\varepsilon_\mathrm{SA} > 0{,}01$, poměr stability blízko 1) | N-částicový | středování u slabých hierarchií selhává |
| okamžiky zákrytů, světelné křivky, radiální rychlosti | N-částicový | potřebuje polohy na každé dráze |
| vnitřní periody řádu dní, běhy na miliony let a déle | sekulární | N-částicový běh by trval roky výpočtů |
| vývoj po dobu života hvězd (slapy, větry, přenos hmoty) | sekulární | miliardy let |
| velké ansámbly stabilních soustav | sekulární | sekundy na člena |
| nejistota | oba | spočítat jednu soustavu oběma a porovnat |

## Vnější poruchy {#external}

### Plynná fáze vzniku {#gas}

Vícenásobné hvězdy vznikají v plynu a jejich dráhy se mění, dokud je plyn přítomen. SWINGERsim
popisuje tuto fázi parametrickým tlumením podle Papaloizoua & Larwooda (2000): pro každou dráhu
uvedenou v `physics.gas.orbits` platí

$$\frac{da}{dt} = -\frac{f(t)\,a}{\tau_a}, \qquad \frac{de}{dt} = -\frac{f(t)\,e}{\tau_e}, \qquad
\frac{d\mathbf{k}}{dt} = \frac{f(t)}{\tau_i}\left[\mathbf{n} - (\mathbf{n}\cdot\mathbf{k})\,\mathbf{k}\right],$$

kde $\tau_a$, $\tau_e$, $\tau_i$ jsou časové škály migrace, cirkularizace a srovnání rovin
(`tau_a_yr`, `tau_e_yr`, `tau_i_yr`), $\mathbf{k}$ normála dráhy, $\mathbf{n}$ normála plynného disku
(počáteční invariantní rovina nebo zadaná rovina) a $f = 1$ až do `t_end_yr`,
$f = \exp[-(t - t_\mathrm{end})/t_\mathrm{decay}]$ potom. Poslední rovnice natáčí rovinu dráhy
k rovině plynu, aniž by měnila $e$, takže $\tan(i/2)$ vůči rovině plynu klesá jako $e^{-t/\tau_i}$.
Velké poloosy se stávají stavovými proměnnými. Tato fáze běží v sekulárním backendu. Časové škály jsou
vstupy: závisejí na hmotnosti a struktuře plynu, které SWINGERsim nemodeluje.

### Slapy Galaxie {#galactic-tide}

Pro dráhy širší než asi 1000 au jsou důležité slapy galaktického disku. V okolí Slunce převládá jejich
svislá složka s potenciálem $2\pi G\rho_0 z_\mathrm{G}^2$ ($z_\mathrm{G}$ je výška podél galaktického
pólu, $\rho_0$ místní hustota disku, výchozí 0,185 M☉ pc⁻³; Heisler & Tremaine 1986; Grishin & Perets
2022, rov. 3–4). Po středování přes dráhu s redukovanou hmotností $\mu$

$$H = \pi G\rho_0\,\mu\,a^2\left[(1-e^2) - (\mathbf{j}\cdot\hat{\mathbf{z}})^2 + 5(\mathbf{e}\cdot\hat{\mathbf{z}})^2\right],$$

kde $\mathbf{j} = \sqrt{1-e^2}\,\mathbf{k}$, $\mathbf{e}$ je vektor excentricity a $\hat{\mathbf{z}}$
galaktický pól. Slapy pohánějí cykly podobné ZLK nad kritickým sklonem $\cos^2 i < 4/5$ až
k $e_\mathrm{max} = \sqrt{1 - (5/4)\cos^2 i}$, na časové škále
$\tau_\mathrm{GT} = 0{,}89\ \mathrm{Gyr}\,(M/\mathrm{M}_\odot)^{1/2}(a/10^4\,\mathrm{au})^{-3/2}$
(Grishin & Perets 2022, rov. 6–7). V sekulárních backendech se tyto členy přidávají k $H$ každé
dráhy; v N-částicovém backendu se slapy uplatňují jako rychlostní impulzy
$\Delta\mathbf{v}_i = -4\pi G\rho_0\,[(\mathbf{r}_i - \mathbf{r}_\mathrm{cm})\cdot\hat{\mathbf{z}}]\,\hat{\mathbf{z}}\,\Delta t$
každou 1/400 vnější periody (operátorové štěpení) a jejich potenciální energie vstupuje do kontroly
energie. Široká dvojhvězda se sklonem 60° dosáhne $e_\mathrm{max}$ s odchylkou do $10^{-4}$ od
analytické hodnoty.

### Procházející hvězdy {#encounters}

Hvězda procházející ve vzdálenosti pericentra $Q$ mnohem větší než dráha mění excentricitu
a orientaci dráhy, její energii však téměř ne. V prvním řádu poruchy plynou změny $\mathbf{e}$
a $\mathbf{j}$ z práce Hamerse & Samsinga (2019, rov. 9) s malým parametrem
$\epsilon = [M^2/(m(m+M))\,(a/Q)^3]^{1/2}(1+E)^{-3/2}$ ($M$ je hmotnost rušitele, $m$ hmotnost dráhy,
$E$ excentricita hyperboly rušitele). S `physics.encounters` losuje SWINGERsim poissonovskou
posloupnost setkání (četnost $n\pi b_\mathrm{max}^2\langle v\rangle$, relativní rychlosti
z Maxwellova rozdělení váženého $v$, srážkové parametry rovnoměrně v ploše, izotropní orientace).
Sekulární backendy uplatňují impulzy prvního řádu na každou dráhu, pro niž je průchod vzdálený
($Q > 3a(1+e)$), a bližší průchody jen započítají, aniž by je uplatnily. N-částicový backend
integruje každou procházející hvězdu jako částici od 20 pericentrových vzdáleností před pericentrem
do 20 za ním, takže sleduje i blízké průchody. Při $Q = 20a$ souhlasí impulz prvního řádu s přímým
N-částicovým průchodem na 0,6 %.

### Gravitační vlny {#gravitational-waves}

S `physics.gw` ztrácí každá dráha dvou osamocených hvězd energii a moment hybnosti gravitačním
zářením, středovaným přes oběh podle Peterse (1964):

$$\dot a = -\frac{64}{5}\frac{G^3 m_1 m_2 M}{c^5 a^3 (1-e^2)^{7/2}}\left(1+\tfrac{73}{24}e^2+\tfrac{37}{96}e^4\right),\qquad
\dot e = -\frac{304}{15}\frac{e\,G^3 m_1 m_2 M}{c^5 a^4 (1-e^2)^{5/2}}\left(1+\tfrac{121}{304}e^2\right),$$

kde $m_1$, $m_2$ jsou hmotnosti a $M = m_1 + m_2$. Sekulární backendy integrují tyto rychlosti spolu
s drahami; N-částicový backend je uplatní jednou za výstupní krok při pevné orientaci a střední
anomálii. Zkracování dráhy Hulseova–Taylorova pulzaru je reprodukováno na 0,02 % (Weisberg & Huang
2016; viz *Testování, část I*). Efekt je důležitý jen pro kompaktní pozůstatky a velmi těsné
dvojhvězdy.

## Implementace {#implementation-dynamics}

**Externí knihovny.** REBOUND (Rein & Liu 2012) s integrátorem IAS15 (Rein & Spiegel 2015),
hledáním srážek a variačními rovnicemi pro MEGNO; REBOUNDx (Tamayo et al. 2020) s efekty `gr_full`
(1PN pro všechna tělesa), `tides_spin` (Lu et al. 2023; další kapitola) a operátorem `modify_mass`
(ztráta hmoty); integrátor DOP853 ze SciPy pro sekulární rovnice; Numba, je-li nainstalována,
překládá sekulární jádro 2+2 SWINGERsimu. **Vlastní kód.** Sekulární pohybové rovnice a středované
energie pro 2+2 (Vokrouhlický 2016, rov. 12–25) a pro libovolný vnořený strom (členy dvojic drah
podle Hamerse et al. 2015 a Hamerse & Portegiese Zwarta 2016) byly naprogramovány pro SWINGERsim,
spolu s členy hamiltoniánu relativity a deformace hvězd závislými jen na excentricitě, které
zachovávají $H$, disipací s konstantním časovým zpožděním v sekulárních rovnicích, tlumením plynem,
slapy Galaxie (sekulární členy a N-částicové impulzy), posloupností setkání a impulzy prvního řádu,
přes oběh středovanými rychlostmi ztrát gravitačními vlnami, rozpoznáním hierarchie, sledováním
extrémů mezi výstupy, diagnostikami dvojic, bilancí energie (1PN hamiltonián, rotační a slapová
energie, moment hybnosti) a během zpět v čase.
