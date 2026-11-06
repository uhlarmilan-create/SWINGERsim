# Stabilita: test duchem, kritéria a klasifikátory {#stability}

Ansámbl vícenásobných hvězd obsahuje mnoho členů, které by nepřežily: jejich dráhy jsou příliš
blízko u sebe a během několika vnějších oběhů se rozpadnou. Aby měl dlouhodobý vývoj člena vůbec
smysl, musí SWINGERsim nejprve rozhodnout, zda je člen stabilní. Nabízí tři nástroje s rostoucí
výpočetní cenou:

1. **Algebraická kritéria** uplatněná na vnořené trojice hierarchie: MA01 (Mardling & Aarseth 2001)
   a kritérium Vynatheyi et al. (2022, V22), které zahrnuje závislost na vzájemném sklonu. Stojí
   mikrosekundy a jsou vlastním kódem SWINGERsimu.
2. **Neuronové klasifikátory** (vícevrstvé perceptrony, MLP) natrénované na čtyřhvězdách 2+2 a 3+1.
   SWINGERsim si trénovací množiny sestavuje sám (N-částicové testy duchem 10⁴ čtyřhvězd pro každou
   hierarchii) a síť trénuje knihovnou scikit-learn; knihovna neuronových sítí je externí, trénovací
   data a použití jsou vlastní SWINGERsimu.
3. **Test duchem**, definice stability podle Vynatheyi et al. (2023, oddíl 2): soustava a její kopie,
   „duch“, jejíž vnitřní velká poloosa je větší o 10⁻⁶ au, se integrují N-částicově po dobu 100
   vnějších oběhů (newtonovské hmotné body); soustava je nestabilní, stane-li se nevázanou nebo
   liší-li se vnitřní velké poloosy obou kopií o více než 1 %. To je reference.

## Reprodukce Vynatheyi et al. (2023) {#vynatheya2023}

**Reference.** Vynatheya et al. (2023) natrénovali klasifikátory MLP na 4·10⁵ čtyřhvězdách
označených testem duchem a publikovali **skóre** – podíl testovacích členů klasifikovaných správně
– pro MA01 a pro svůj MLP pro čtyřhvězdy 2+2 a 3+1 (jejich tabulky 3–4). Níže se objevují ještě
dvě míry: **úplnost** (recall) stabilní třídy (podíl skutečně stabilních členů, které klasifikátor
označí za stabilní) a **přesnost** (precision; podíl členů označených za stabilní, které jsou
skutečně stabilní).

**Vstup.** Dvě trénovací množiny po 10⁴ čtyřhvězdách, losované z apriorního rozdělení článku (jejich
oddíl 3.1: hmotnosti, poměry period, excentricity, sklony), klasifikované testem duchem s limitem
30 min na člena; MLP s 11 vstupními parametry a architekturou článku (4 skryté vrstvy po 50
logistických jednotkách, optimalizátor Adam, regularizace L2 10⁻⁴; jejich oddíl 3.2), průměr pěti
sítí, trénováno na 80 % členů a hodnoceno na zbývajících 20 %. Třetí množina 10⁴ čtyřhvězd 3+1 byla
losována, jako v článku, jen blízko hranice stability MA01.

**Tabulka 43.** Vynatheya et al. (2023): skóre proti testu duchem. Tyto čtyři testy jsou ruční;
MLP je třídy ext (síť scikit-learn na trénovacích množinách SWINGERsimu), MA01 třídy own. Sloupce jsou definovány v oddílu [Jak číst tuto zprávu](#how-to-read).

| Test | Třída | Očekáváno | Spočteno | Rozdíl | Tolerance | Verdikt |
|---|---|---|---|---|---|---|
| `vyn23-2p2-ma01-score` | own | 0,83 | 0,869 | +0,039 | 0,05 | PASS |
| `vyn23-2p2-mlp-score` | ext | 0,94 | 0,9375 | −0,0025 | 0,03 | PASS |
| `vyn23-3p1-ma01-score` | own | 0,56 | 0,7896 | +0,23 | 0,05 | **EXPLAINED** |
| `vyn23-3p1-mlp-score` | ext | 0,93 | 0,9195 | −0,0105 | 0,03 | PASS |

**Tabulka 44.** Tři trénovací množiny. Sloupce 2–4: podíl členů, které jsou stabilní podle testu
duchem, podle MA01 a podle V22; sloupce 5–7: skóre MA01, V22 a MLP proti testu duchem.

| Množina | Stabilní (duch) | Stabilní (MA01) | Stabilní (V22) | Skóre MA01 | Skóre V22 | Skóre MLP |
|---|---|---|---|---|---|---|
| 2+2, apriorní rozdělení článku | 28,2 % | 39,3 % | 34,9 % | 0,869 | 0,895 | 0,938 |
| 3+1, apriorní rozdělení článku | 9,5 % | 28,5 % | 26,3 % | 0,790 | 0,812 | 0,920 |
| 3+1, blízko hranice MA01 | 21,0 % | – | – | 0,509 | 0,568 | – |

![**Obrázek 22.** Skóre tří klasifikátorů – MA01 (oranžově), V22 (zeleně) a neuronové sítě MLP (modře) – proti N-částicovému testu duchem pro množinu 2+2 (vlevo) a první množinu 3+1 (vpravo). Svislá osa je skóre, podíl testovacích čtyřhvězd, jejichž stabilita je předpovězena správně (1 = dokonale); číslo nad každým sloupcem je opakuje. Černé čárkované čáry jsou skóre publikovaná Vynatheyou et al. (2023), kteří pro V22 žádné neuvádějí.](../../figures/part1/stability.png)

**Výsledek (tabulky 43 a 44).** Na obrázku 22 dosahují modré sloupce v obou množinách čárkovaných
publikovaných čar: MLP dosahuje skóre 0,938 (2+2) a 0,920 (3+1) proti publikovaným 0,94 a 0,93 –
v rozmezí 0,01, přestože byl natrénován na 40krát menším počtu soustav. MA01 dosahuje na množině
2+2 skóre 0,869 proti publikovaným 0,83 (+0,04, v rámci tolerance 0,05 dané rozdílem mezi
trénovacími množinami). Na první množině 3+1 však MA01 dosahuje 0,79 proti 0,56. Tento test je
vysvětlen (EXPLAINED) vzorkováním testovací množiny, a to ze tří důvodů. Za prvé, při pouhých 9,5 %
stabilních soustav dosáhne klasifikátor, který označí za nestabilní *všechno*, skóre 0,905; MA01
označí za stabilní 28,5 % – příliš mnoho –, rozpozná však většinu z mnoha nestabilních členů, takže
jeho skóre je vysoké. Za druhé, článek losoval své členy blízko hranice MA01 a jeho vlastní
tabulka 4 (přesnost a úplnost) implikuje skóre 0,69, nikoli 0,56. Za třetí, množina losovaná blízko
hranice jako v článku dává 0,51, blízko publikovaných 0,56 (a V22 0,57). Tato množina z okna má
méně stabilních členů (21 %) než asi 36 %, které implikují publikované přesnosti a úplnosti, což
snižuje přesnost stabilní třídy (0,29 proti publikovaným 0,54; úplnost 0,94 proti 0,95). Skóre
závisí na tom, jak je testovací množina losována, a má smysl jen spolu se vzorkováním. Žádný člen
z 3·10⁴ testů duchem nedosáhl časového limitu.

![**Obrázek 23.** Diagnostika neuronové sítě pro 3+1. *Vlevo nahoře:* trénovací ztráta každé z pěti sítí v závislosti na trénovací epoše. *Nahoře uprostřed:* skóre na validačních členech odložených během trénování v závislosti na epoše. *Vlevo dole:* křivka ROC – pro každý práh předpovězené pravděpodobnosti podíl skutečně stabilních členů označených za stabilní (svisle) v závislosti na podílu nestabilních členů označených za stabilní (vodorovně); tečkovaně: náhodný klasifikátor; plocha pod křivkou je 0,914. *Dole uprostřed:* spolehlivost – z členů s předpovězenou pravděpodobností $p$ podíl, který je skutečně stabilní (tečkovaně: ideál). *Vpravo dole:* histogram předpovězených pravděpodobností odložených členů podle jejich skutečné třídy (modře stabilní, oranžově nestabilní).](../../figures/studies/vyn23_3p1_mlp.png)

**Kalibrace.** Obrázek 23 ukazuje, že se síť pro 3+1 naučila (ztráta klesá a ustálí se; validační
skóre roste z 0,90, skóre klasifikace všeho jako nestabilního, na 0,92–0,93) a že její
pravděpodobnosti lze brát doslova: body spolehlivosti leží na diagonále. Histogram ukazuje
omezení: málo stabilních soustav 3+1 dostane pravděpodobnost nad 0,5, takže úplnost stabilní třídy
při prahu 0,5 je jen 0,20 (přesnost 0,65; pro 2+2: úplnost 0,90, přesnost 0,88). U soustav 3+1 je
třeba používat samotnou pravděpodobnost, nikoli označení při prahu 0,5.

**Zjištění o kritériích.** Na širokých množinách označují MA01 a V22 za stabilní více soustav než
test duchem (2+2: 39,3 % a 34,9 % proti 28,2 %; 3+1: 28,5 % a 26,3 % proti 9,5 %). Jsou optimistická:
soustava, kterou označí za stabilní, nemusí být stabilní po dobu 100 vnějších oběhů. V ansámblech
proto SWINGERsim uvádí algebraické verdikty jako diagnostiku a tam, kde na stabilitě záleží,
používá N-částicový backend (nebo test duchem), jako v [kapitole 6](#czev4315).

## Mapy chaosu a simulační inference {#chaos-sbi}

**Mapy chaosu.** SWINGERsim počítá MEGNO na mřížce dvou dráhových parametrů, aby zmapoval
pravidelné a chaotické oblasti. Srovnávacím testem je vnitřní mez stability cirkumbinární dráhy
kolem kruhové dvojhvězdy se stejnými hmotnostmi, publikovaná jako mřížka Quarlesem et al. (2018):
2,27 vzdálenosti složek dvojhvězdy (aproximace Holmana & Wiegerta 1999 dává 2,39). SWINGERsim, který
vyžaduje, aby bylo všech 8 oběžných fází pravidelných, nachází 2,30
(`chaosmap-circumbinary-limit`, třída R; tolerance 0,06, krok mřížky 0,05 plus hrubší vzorkování
8 fází a 10⁴ oběhů dvojhvězdy místo 10⁵).

**Simulační inference.** Pro aposteriorní rozdělení, kterých ansámbl nemůže dosáhnout vážením, umí
SWINGERsim natrénovat neuronový odhad aposteriorního rozdělení pomocí externího balíku `sbi`
(Tejero-Cantero et al. 2020) na svých vlastních simulacích. Srovnávacím testem je lineární
gaussovský model s analytickým aposteriorním rozdělením: s 2500 trénovacími simulacemi je šířka
aposteriorního rozdělení 0,291 proti přesné hodnotě 0,309 (−5,8 %, tolerance 10 %;
`sbi-linear-gaussian-std`, třída ext). Tento test ověřuje, že SWINGERsim předává své simulace
balíku správně; samotný algoritmus inference je externí.
