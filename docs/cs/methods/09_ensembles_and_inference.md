# Ansámbly a inference {#ensembles-and-inference}

Parametry skutečné vícenásobné hvězdy nejsou nikdy známy všechny. Fotometrie určí sklon zákrytové
dvojhvězdy, nikoli však orientaci jejího uzlu na obloze; vnější dráha bývá známa jen ze signálu
v okamžicích minim; hmotnosti mají nejistoty. Jediná simulace odpovídá na otázku „co kdyby“; aby
odpověděl na otázku „jak pravděpodobně“, spouští SWINGERsim **ansámbly Monte Carlo** – mnoho členů,
z nichž každý je úplným během s parametry vylosovanými z rozdělení – a podmiňuje je měřeními. Tato
kapitola popisuje ansámbly, pozorovací váhy a jejich interpretaci jako bayesovských veličin
a pomocné nástroje strojového učení postavené na ansámblech.

SWINGERsim je zaměřen na **konkrétní soustavy**: konfigurace popisuje jednu pozorovanou vícenásobnou
hvězdu a ansámbl zkoumá nejistotu právě této soustavy. Emulátory a simulační inference (konec kapitoly)
jsou pomocné nástroje, které práci s takovými ansámbly zlevňují; nejsou aparátem pro populační
syntézu.

![**Obrázek 20.** Ansámbl Monte Carlo ve třech krocích (ilustrativní). Vlevo: histogram neznámého
parametru vylosovaného pro členy (zde vnější perioda ve dnech; svislá osa: počet členů). Uprostřed:
excentricita vnitřní dvojhvězdy v závislosti na čase v letech pro vzorek členů, každý jednou barevnou
čarou – každý člen je úplný běh. Vpravo: histogram odpovědi na položenou otázku (zde největší dosažená
excentricita; svislá osa: počet členů). Výsledkem je rozdělení vpravo, nikoli kterýkoli jednotlivý
běh.](../../figures/concepts/ensemble_idea.png)

Obrázek 20 ukazuje logiku postupu: nejisté vstupy (vlevo) se stanou svazkem možných historií
(uprostřed) a potom rozdělením zkoumané veličiny (vpravo). Každé tvrzení odvozené z ansámblu je
tvrzením o podílu nebo kvantilu takového rozdělení, se statistickou chybou.

## Co se obměňuje {#what-is-varied}

Konfigurace popisuje nominální soustavu; `ensemble.vary` přiřazuje parametrům (`masses.m1`,
`orbits.out.P`, `orbits.B.Omega`, `stars.m2.R`, …) rozdělení, která tvoří **apriorní rozdělení**
ansámblu. Rozdělení uvádí tabulka 14.

**Tabulka 14.** Rozdělení v `ensemble.vary`. *Zápis* je způsob, jak se zadává; *Rozdělení* co se losuje;
*Typické použití* parametry, pro něž se hodí. Každé rozdělení lze oříznout pomocí `min`, `max` (hodnota
mimo meze se losuje znovu).

| Zápis | Rozdělení | Typické použití |
|---|---|---|
| číslo | pevná hodnota | – |
| `uniform: [a, b]` | rovnoměrné mezi $a$ a $b$ | excentricity, fáze |
| `loguniform: [a, b]` | rovnoměrné v $\log x$ | periody, vzdálenosti přes několik řádů |
| `normal: [μ, σ]` | gaussovské | naměřené hmotnosti, periody, poloměry |
| `isotropic` nebo `isotropic: [i1, i2]` | $\cos i$ rovnoměrně (mezi mezemi) | sklon dráhy neznámé orientace |
| `angle` | rovnoměrné v 0°–360° | uzly, argumenty periastra, střední anomálie |
| `choice: [x1, x2, …]` | jedna z hodnot se stejnou pravděpodobností | diskrétní alternativy |
| `mirror: true` (u sklonu) | $i \to 180° - i$ s pravděpodobností 1/2 | smysl oběhu zákrytové dvojhvězdy |

**Úhly v soustavě oblohy.** Pro dráhu neznámé orientace dává izotropní sklon spolu s rovnoměrným uzlem
izotropní rozdělení její normály a vzájemné sklony plynou z geometrie
([Hierarchie a dráhy](#mutual-inclination)). U změřené zákrytové dvojhvězdy se $i$ losuje kolem
naměřené hodnoty a `mirror: true` přidává neznámý smysl oběhu; záleží jen na relativním smyslu dvou
drah, takže první dvojhvězdu lze ponechat s $i < 90°$. V populační studii (například Pejcha et al. 2013)
je vztažnou rovinou místo toho vnější dráha.

**Požadavky.** `require` znovu losuje členy, které nesplňují kritérium: `stable_ma01`, `stable_v22`
([Stabilita a osudy](#criteria)) nebo `hierarchical` (každé vnitřní apocentrum leží uvnitř pericentra
dráhy kolem něj). Požadavek ořezává apriorní rozdělení; výsledky jsou pak na něm podmíněny.
`ma01_unstable_keep` a `ma01_window` utvářejí trénovací množiny pro klasifikátory stability (jako
Vynatheya et al. 2023 ponechali 0,3 losování nestabilních podle MA01 pro 2+2 a 0,2 pro 3+1).

**Směry.** S `directions: [forward, backward]` se každý člen integruje oběma směry ze stejného
současného stavu: běh zpět říká, jak dlouho soustava vypadá tak, jak vypadá nyní (například od kdy
zakrývají obě dvojhvězdy), běh dopředu, co přijde dál, a přežití v obou směrech je silnou podmínkou
([Stabilita a osudy](#survival)).

**Souhrn a reprodukovatelnost.** Každý člen se zredukuje na jeden řádek: osud, extrémy excentricity,
vzájemné sklony a jejich rozsahy, podíly času se zákryty a změny zákrytů, apsidální periody, vyplnění
Rocheova laloku, verdikty stability (MA01, V22, volitelně klasifikátor MLP a test duchem) a současné
rychlosti změn. Členy závisejí jen na konfiguraci a na `seed`, takže ansámbl lze rozšířit nebo znovu
spočítat člen po členu jiným backendem.

### Současné rychlosti změn {#rates-today}

Nejinformativnějšími měřeními precedující soustavy jsou její současné rychlosti změn, například
rychlost změny sklonu zákrytové dvojhvězdy změřená z měnících se hloubek zákrytů (Zasche et al. 2025).
SWINGERsim počítá `<orbit>_didt0_deg_yr`, rychlost změny sklonu k obloze, a `<orbit>_dtiltdt0_deg_yr`,
rychlost změny *odklonu* $|i - 90°|$, lineárním proložením od $t = 0$ přes dvě periody rušící dráhy
(nejméně šest výstupů) nebo přes `run.rate_window_yr`, aby se kmity oskulačních elementů s vnější
periodou vystředovaly. U zákrytové dvojhvězdy je užitečná rychlost změny odklonu: prohlubující se
zákryty znamenají zápornou rychlost změny odklonu bez ohledu na neznámý smysl oběhu, kdežto znaménko
$di/dt$ na něm závisí. Běh dopředu dává jednostrannou rychlost dopředu, běh zpět dozadu; pozorovatel
měří rychlost na základně se středem v současnosti, takže podmínění standardně používá *centrální*
rychlost, průměr obou. Když sklon rychle kmitá, mohou se jednostranné rychlosti lišit o několik stupňů
za rok; `rate_window_yr` má být zhruba polovinou pozorovací základny.

## Statistická chyba podílu {#wilson}

Podíl $p$ odhadnutý z $n$ členů má binomickou nejistotu přibližně $\sqrt{p(1-p)/n}$: ±2,5 % pro
$p$ = 0,3 a $n$ = 400. SWINGERsim uvádí Wilsonův interval (Wilson 1927) na úrovni 68 %, který zůstává
uvnitř 0–1 i pro podíly blízké 0 nebo 1. Pro spojité sloupce se uvádějí kvantily (medián, 16 % a 84 %).
Členy se slabou hierarchií bývají chaotické: dva členy s téměř stejnými parametry skončí odlišně, takže
smysl mají jen rozdělení přes mnoho členů.

## Podmínění pozorováními {#weights}

Bez měření toho, jak se soustava mění, jsou ansámbly pozorované soustavy dopředu a zpět statisticky
stejné. Měření tuto symetrii narušují. Každý člen $i$ dostane váhu gaussovské věrohodnosti

$$ w_i = \exp\left[-\frac{1}{2}\sum_k \left(\frac{x_{ik} - x_k^\mathrm{obs}}{\sigma_k}\right)^2\right] , $$

kde $x_{ik}$ je hodnota souhrnného sloupce $k$ pro člen $i$ (například rychlost změny odklonu B nebo
apsidální perioda) a $x_k^\mathrm{obs} \pm \sigma_k$ měření (`ensemble.observed`, položky
`{value, sigma}`). Pevné meze (`{min, max}` nebo `{equals: …}` pro textové sloupce, jako je osud)
nastaví $w_i = 0$ pro členy, které odporují tomu, co je s jistotou známo – například
`backward:m1_roche_fill_max: {max: 1}` (hvězdy pozorované oddělené dvojhvězdy v minulosti nevyplnily
své Rocheovy laloky) nebo `backward:fate: {equals: intact}`. Člen bez konečné hodnoty pozorovaného
sloupce dostane váhu 0. Názvy sloupců mohou nést předponu `forward:`, `backward:` nebo `central:`,
která říká, ze kterého běhu hodnota pochází.

Všechny statistiky zprávy jsou pak vážené, oba směry používají váhu téhož člena a vylosované parametry
se uvádějí před vážením (apriorní) a po něm (podmíněné). Vážené kvantily se počítají z kumulativních
vah.

**Interpretace.** Jde o *vzorkování podle důležitosti* (importance sampling) aposteriorního rozdělení
s vylosovanými rozděleními jako apriorním rozdělením: podle Bayesovy věty
$p(\theta \mid x^\mathrm{obs}) \propto L(x^\mathrm{obs} \mid \theta)\,\pi(\theta)$ a členy vylosované
z apriorního rozdělení $\pi$ a vážené věrohodností $L$ reprezentují aposteriorní rozdělení. Obsaženy
jsou dva předpoklady: chyby měření jsou gaussovské a nezávislé a simulované souhrnné hodnoty se berou
bez chyby modelu (rozptyl mezi členy reprezentuje neznámé parametry, nikoli šum modelu).

**Efektivní počet členů.** Vážení soustředí informaci do menšího počtu členů. Efektivní počet
(Kish 1965)

$$ n_\mathrm{eff} = \frac{\left(\sum_i w_i\right)^2}{\sum_i w_i^2} $$

nahrazuje $n$ ve statistických chybách (Wilsonovy intervaly se počítají s ním). Je-li malý (pod
zhruba 20), vybírá měření jen malou část ansámblu a výsledky jsou zašuměné: je třeba vylosovat více
členů nebo zúžit apriorní rozdělení kolem oblasti, která data vystihuje.

![**Obrázek 21.** Podmínění váhami důležitosti v modelovém příkladu (nový obrázek pro tento svazek;
syntetická čísla, nikoli skutečná soustava). Parametr θ – zde v roli vnější periody ve dnech – se
losuje logaritmicky rovnoměrně pro 4000 členů; pozorovatelná veličina každého člena je mocninnou
funkcí θ s 30% náhodným rozptylem (zastupujícím neznámé úhly). (a) Simulovaná pozorovatelná veličina
v závislosti na θ, každý člen bod obarvený podle své váhy $w$ (bílá 0, tmavě modrá 1); oranžová čára
a pás jsou měření −1,4 ± 0,2. (b) Hustota pravděpodobnosti θ pro apriorní rozdělení (šedý histogram,
všechny členy se stejnou váhou) a po podmínění (modrý obrys, vážený). Text uvádí $n$, $n_\mathrm{eff}$
a střední váhu $Z$.](../../figures/methods/weights_toy.png)

Obrázek 21 ukazuje, co váhy dělají. V panelu (a) dostanou vysokou váhu jen členy, jejichž simulovaná
hodnota padne do oranžového pásu; kvůli rozptylu pocházejí z celého rozsahu θ, nikoli z jedné hodnoty.
Panel (b) ukazuje výsledek: široké apriorní rozdělení (šedě) se stane úzkým podmíněným rozdělením
(modře) kolem θ ≈ 100–150. Z 4000 členů nese váhu asi 245 ($n_\mathrm{eff}$), takže podmíněné
statistiky mají přesnost ansámblu o zhruba 245 členech, nikoli 4000.

### Střední váha jako bayesovská evidence {#evidence}

Střední váha ansámblu,

$$ Z = \frac{1}{n}\sum_{i=1}^{n} w_i , $$

je odhadem Monte Carlo integrálu věrohodnosti přes apriorní rozdělení,
$\int L(x^\mathrm{obs}\mid\theta)\,\pi(\theta)\,d\theta$, tj. *bayesovské evidence* (marginální
věrohodnosti) hypotézy, kterou ansámbl reprezentuje, až na konstantní normování
$\prod_k (2\pi\sigma_k^2)^{-1/2}$ gaussovské věrohodnosti. Členy odstraněné pevnými mezemi se počítají
s $w = 0$; členy znovu losované kvůli `require` se nepočítají vůbec, takže požadavek je součástí
hypotézy. Její statistická chyba je asi $\mathrm{std}(w)/\sqrt{n}$ a je-li $n_\mathrm{eff}$ malé, je
zašuměná.

**Srovnávání hypotéz.** Pro dvě hypotézy $H_1$, $H_2$ – například čtyřhvězdu 2+2 a pětihvězdu
2+(2+1) nebo dva rozsahy neznámé vnější periody –, podmíněné **týmiž pozorováními se stejnými chybami
a stejnými pevnými mezemi**, se normování vykrátí a poměr

$$ B_{12} = \frac{Z_1}{Z_2} $$

je *Bayesův faktor*: faktor, o který data mění poměr šancí $H_1$ proti $H_2$. Při stejných apriorních
šancích je roven aposteriornímu poměru šancí. V aplikaci na CzeV4315 ve svazku *Testování SWINGERsimu,
část IV* (jedna dvojhvězda s naměřenou rychlostí změny sklonu, druhá s trvalými zákryty) dosáhl model
2+(2+1) hodnoty $Z = 0{,}022$ s $n_\mathrm{eff} = 46$ (1500 N-částicových členů). Model 2+2 s periodou
119 d jako jejich společnou dráhou nedal mezi 1000 členy žádný s váhou nad $10^{-12}$, tj.
$Z \lesssim 10^{-3}$. Bayesův faktor ve prospěch pětihvězdy je proto alespoň asi 20, formálně $10^{13}$.
Když jedna hypotéza nemá (téměř) žádný vyhovující člen, je třeba uvést horní mez $Z < 1/n$, nikoli
formální hodnotu. Platí tři výhrady: (1) evidence závisí na apriorním rozdělení – široké apriorní
rozdělení, které plýtvá členy na oblasti vyloučené daty, snižuje $Z$, což je bayesovská penalizace za
neurčitost, srovnání je však jen tak smysluplné, jak fyzikálně podložená jsou apriorní rozdělení;
(2) malé $n_\mathrm{eff}$ činí samotné $Z$ nejistým až o velký faktor; (3) $Z$ ze sekulárního ansámblu
může být nadhodnoceno, když jsou vyhovující členy na delších časových škálách chaotické nebo
nestabilní, což sekulární backend nevidí – doporučuje se N-částicová kontrola nejlepších členů. Zpráva
ansámblu vypisuje $n_\mathrm{eff}$ a $Z$ s chybou $\mathrm{std}(w)/\sqrt{n}$ a upozorní na případ
méně než dvou vyhovujících členů, pro který platí horní mez.

**Uvádění výsledků.** Podmíněný výsledek má uvádět apriorní rozdělení (rozdělení v `vary` a požadavky),
pozorování a jejich chyby, $n$, $n_\mathrm{eff}$, vážený medián a rozsah 16–84 % každé veličiny,
podíly s jejich Wilsonovými intervaly a při srovnávání hypotéz obě hodnoty $Z$ s jejich
$n_\mathrm{eff}$.

## Emulátory {#emulators}

Ansámbl zobrazuje parametry svých členů na souhrnné sloupce za cenu sekund až hodin na člena. Neuronový
*emulátor* se toto zobrazení naučí z hotového ansámblu a vyhodnotí je za mikrosekundy. Sítě jsou
vícevrstvé perceptrony (scikit-learn, jako u Vynatheyi et al. 2022, 2023); vstupy a regresní cíle se
standardizují (volitelně po $\log_{10}$), úhly lze zobrazit na $(\cos\varphi, \sin\varphi)$ a známý
mocninný zákon lze z cíle vydělit. *Hluboký soubor sítí* (deep ensemble) s různými semínky (výchozí 5;
Lakshminarayanan et al. 2017) udává rozptyl předpovědi; 20 % členů se odloží a slouží k výpočtu skóre.
Dvě použití:

- **Inference.** Vylosovat $10^6$–$10^7$ sad parametrů z rozdělení ansámblu, předpovědět pozorované
  sloupce a každé losování zvážit věrohodností měření, v níž se k rozptylu měření přičte rozptyl
  emulátoru (rozptyl hlubokého souboru sítí plus jeho reziduum na odložených členech). Je to výše
  popsané podmínění s mnohem více členy; přesná měření již nezhroutí $n_\mathrm{eff}$.
- **Klasifikace.** Binární sloupec (verdikt testu duchem, `fate==intact`) se naučí klasifikátor, který
  vrací pravděpodobnost – například stabilitu čtyřhvězdy z jejích bezrozměrných parametrů, podobně jako
  klasifikátory Vynatheyi et al. (2022, 2023).

Emulátor je jen tak dobrý jako jeho ansámbl: mimo vzorkované rozsahy a tam, kde je členů málo,
extrapoluje, což rozptyl hlubokého souboru sítí ukazuje jen částečně.

## Simulační inference {#sbi}

Neuronový odhad aposteriorního rozdělení (NPE; Papamakarios & Murray 2016; Greenberg et al. 2019)
s balíkem `sbi` (Tejero-Cantero et al. 2020) natrénuje podmíněnou hustotu $q(\theta \mid x)$ na dvojicích
vylosovaných parametrů $\theta$ a simulovaných pozorovatelných veličin $x$ ansámblu. Chyby měření
vstupují jako gaussovský šum přičtený k simulovaným pozorovatelným veličinám a pevné meze vybírají
trénovací členy. Poté se aposteriorní rozdělení pro libovolné naměřené $x$ navzorkuje za méně než
sekundu, bez nových simulací. Odpovídá na tutéž otázku jako váhy důležitosti, s rozděleními ansámblu
jako apriorním rozdělením; vážení nepotřebuje trénink, jeho $n_\mathrm{eff}$ se však pro přesná měření
zhroutí, kdežto NPE mezi členy interpoluje. Natrénovaná síť se může mýlit tam, kde má ansámbl málo
členů, a proto se kontroluje proti vážení a testem pokrytí: podílem odložených členů, jejichž skutečné
hodnoty leží uvnitř jejich 68% intervalů (ideálně 0,68), a simulační kalibrací (Talts et al. 2018),
při níž je pořadí skutečné hodnoty mezi vzorky aposteriorního rozdělení pro kalibrované aposteriorní
rozdělení rovnoměrné.

## Implementace {#implementation-ensembles}

**Externí knihovny.** Vícevrstvé perceptrony ze scikit-learn pro emulátory a klasifikátory; balík `sbi`
(Tejero-Cantero et al. 2020) s PyTorch pro neuronový odhad aposteriorního rozdělení; natrénované
klasifikátory Vynatheyi et al. (2023), jsou-li nainstalovány. **Vlastní kód.** Losování členů (všechna
rozdělení, ořezání, zrcadlení, požadavky, utváření trénovacích množin) s deterministickými semínky;
paralelní spouštění s možností navázání; redukce každého běhu na jeho souhrnný řádek (včetně
centrálních rychlostí změn); Wilsonovy intervaly; pozorovací váhy s gaussovskými členy, pevnými mezemi
a předponami směru; $n_\mathrm{eff}$ a vážené kvantily; vrstva emulátorů (transformace vstupů a cílů,
hluboké soubory sítí, přetrénování sítí uvízlých na plató, rozptyl emulátoru ve věrohodnosti); a rozhraní
SBI s přidáváním šumu a testy pokrytí a kalibrace.
