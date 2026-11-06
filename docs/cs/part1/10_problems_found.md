# Problémy zjištěné během testování a jak se s nimi zachází {#problems-found}

Běh srovnávacích testů, katalogu a ansámblů až s 5·10⁴ členy ověřil program stejně jako fyziku. Tato
kapitola shrnuje, co testy odhalily o numerice a praktickém používání SWINGERsimu, popsané jako
vlastnosti programu, které uživatel potřebuje znát: co se může pokazit, jak to SWINGERsim zjistí
nebo tomu zabrání a co by měl uživatel udělat. Druhá skupina zjištění se týká samotných
publikovaných referencí; je uvedena na konci, protože ovlivňuje, jak by měl být s literaturou
srovnáván jakýkoli program.

Tabulka 45 zjištění shrnuje; každý řádek je podrobně rozebrán v oddílu uvedeném v jeho posledním
sloupci.

**Tabulka 45.** Shrnutí zjištění. Sloupec 2 říká, co se může stát, sloupec 3 jak s tím SWINGERsim
zachází, sloupec 4 kde v této zprávě se zjištění objevuje.

| Oblast | Co se může stát | Jak s tím SWINGERsim zachází | Výskyt |
|---|---|---|---|
| nestabilní členy | těsné setkání způsobí zhroucení časového kroku a běh nikdy neskončí | hledání srážek platné v obou směrech času; limit skutečného času kontrolovaný po každém kroku; osud `wall_time` | [kap. 6](#czev4315) |
| krátké vývojové fáze | úseky vývojového backendu přeskočí tepelné pulzy | hranice úseků v uzlech vývojových drah, zjemněné na 1 % v poloměru | [kap. 8](#giant-tertiary) |
| oskulační elementy | oskulační periody nejsou periodami zákrytů | střední pohyb přes oběh pro úseky světelné křivky; ukotvení period zákrytů při prokládání | [kap. 5](#qzcar-minima) |
| přenos hmoty na excentrických drahách | třetí hvězda zvyšuje excentricitu páru, v němž probíhá přenos | varování; slapy dráhu zakulatí; kontrola počátku počítá s vynucenou excentricitou | [kap. 8](#algol-triple) |
| chaotické a špatně podmíněné běhy | dvě správné integrace se rozbíhají | srovnávat statistiky ansámblů, nikoli jednotlivé členy | [kap. 2](#vokrouhlicky2016) |
| sekulární backend blízko meze stability | středování přes oběh selhává | diagnostika stability v každém běhu; N-částicová kontrola nejlepších členů | [kap. 6](#j075203) |
| váhy ansámblu | celou váhu nese jeden člen | efektivní počet členů uváděný s každou evidencí | [kap. 6](#czev4315) |
| nedefinované veličiny | apsidální perioda téměř kruhové dráhy | uváděna jako neměřitelná | [kap. 7](#catalogue-precession) |
| vývojové dráhy hvězd | současné hmotnosti dosáhnou dvě vývojové dráhy; hmotnost mimo mřížku; hvězda přeplňující svůj lalok na začátku | varování a nejméně vyvinutá dráha; člen odmítnut s uvedením důvodu | [kap. 5](#qzcar-evolution) |
| tabulky atmosfér | teploty nebo tíhová zrychlení mimo mřížky okrajového ztemnění | nejbližší uzel s varováním | [kap. 8](#life-pulsators) |
| zachovávané veličiny | newtonovská energie se s relativitou nezachovává | chyby měřené vůči hamiltoniánu 1PN | [kap. 2](#nbody-basics) |

## Nestabilní členy a limity skutečného času {#wall-clock}

Ansámbl losovaný ze širokého apriorního rozdělení vždy obsahuje členy blízko meze stability nebo za
ní. V takovém členu mohou dvě hvězdy padat k sobě. Je-li srážka zjištěna, běh končí osudem
*srážka*. Testy ukázaly, že výchozí („direct“) hledání srážek v REBOUND přehlíží páry, které se
k sobě přibližují, když integrace běží v čase dozadu; dvě hvězdy pak projdou jedna druhou,
adaptivní časový krok IAS15 se neomezeně zmenšuje (až na 10⁻⁴⁸ roku) a běh neskončí; s přímým
hledáním se takto chová jedenáct z 10⁴ běhů ansámblu CzeV4315 se slapy. SWINGERsim proto pro běhy
dozadu používá hledání srážek, které funguje v obou směrech času (hledání „line“), a limit
skutečného (wall-clock) času (`run.max_wall_s`) kontroluje po *každém* integračním kroku, nejen mezi
výstupními okamžiky. Člen, který limitu dosáhne, končí osudem `wall_time` a je tak i uveden; není
potichu vypuštěn. Při integraci s hledáním „line“ končí všech jedenáct těchto členů srážkou, což je
jejich správný osud. Testy duchem v [kapitole 9](#stability) mají vlastní limit 30 minut na člena;
žádný z 3·10⁴ testů duchem jej nedosáhl. Co by měl uživatel udělat: limit ponechat (chrání ansámbl,
ne fyziku), číst osudy a členy `wall_time` považovat za nerozhodnuté, nikoli za stabilní.

## Délka kroku vývojového backendu {#segments}

Vývojový backend posouvá dráhy v úsecích, v nichž se poloměry hvězd mění o méně než 1 %. Na tepelně
pulzující AGB vyskočí poloměr hvězdy při každém pulzu 2–4krát jen na několik set let, a to každých
10⁴–10⁵ let. Rovnoměrná časová mřížka pulzy přeskočí: na takové mřížce dosáhne obr jako třetí
hvězda z [kapitoly 8](#giant-tertiary) nejvýše 313 R☉ a 73 % svého Rocheova laloku a skončí na konci
své vývojové dráhy, zatímco sama vývojová dráha MIST dosahuje asi 500 R☉. SWINGERsim proto klade
hranice úseků do uzlů vývojové dráhy MIST a zjemňuje je, dokud se mezi dvěma hranicemi žádný
poloměr nezmění o více než 1 %; tak tentýž obr přeteče svůj lalok při 452 R☉ na vrcholu pulzu.
Stejné pravidlo chrání každou krátkou vývojovou fázi (Hertzsprungovu mezeru, heliové záblesky).
Cenou je doba běhu na AGB; zjemnění nelze vypnout.

## Oskulační elementy a periody zákrytů {#osculating}

Dráhové elementy SWINGERsimu jsou **oskulační**: spočtené z okamžitých poloh a rychlostí.
Ve vícenásobné soustavě a ve dvojhvězdě slapově deformovaných hvězd není oskulační perioda
periodou, s níž se zákryty opakují. V testech se projevily dva důsledky. (i) Hybrid PHOEBE počítá
světlo deformovaných hvězd po úsecích N-částicové dráhy; kdyby každý úsek používal oskulační
periodu slapově deformované dvojhvězdy (u QZ Car Ac odchylka 137 s za oběh), jeho minima by se
posouvala o desítky sekund. SWINGERsim používá střední pohyb přes jeden oběh a minima z PHOEBE
souhlasí s kulovým modelem na několik sekund ([kapitola 5](#qzcar-minima)). (ii) Publikované elementy
zadané doslova do N-částicového modelu dávají zákryty, které se během let posunou o hodiny; část II
ukazuje, jak SWINGERsim *ukotvuje* periodu a fázi každého zákrytového páru k pozorovaným minimům.
Související bod: dobu šíření světla je třeba na polohy uplatňovat konzistentně. S přesnými
retardovanými polohami sledují minima QZ Car barycentrické konjunkce posunuté jen o Rømerovo
zpoždění uvnitř dvojhvězdy (±16 s); oprava poloh v prvním řádu by je posunula asi o 200 s. Hlídají
to testy Rømerova zpoždění v [kapitole 4](#timings).

## Přenos hmoty na excentrické dráze {#eccentric-transfer}

Přenos hmoty ve SWINGERsimu je formulován pro kruhové dráhy, stejně jako v BSE. Třetí hvězda však
může páru, v němž probíhá přenos, vnutit malou excentricitu nebo ji zvyšovat cykly ZLK. Z toho
plynuly dva jevy. Za prvé, hvězda na mírně excentrické dráze vyplní svůj lalok v periastru dříve
než na kruhové (v trojhvězdě typu Algol o 11 kyr); kontrola počátku s vynucenou excentricitou
počítá, takže běh se nezastaví před přenosem. Za druhé, bez slapů zvýší třetí hvězda v trojhvězdě
typu Algol excentricitu páru, v němž probíhá přenos, na 0,44 (0,33 na konci běhu), mimo platnost
kruhové formulace; SWINGERsim vypíše varování, které uvádí dvojhvězdu a excentricitu. Se zapnutými
slapy, jako v běhu z [kapitoly 8](#algol-triple), zůstává pár kruhový ($e$ < 4·10⁻⁴) a hmotnosti na
konci jsou stejné (0,454 a 3,172 M☉). Co by měl uživatel udělat: zapnout slapy pro každou
dvojhvězdu, která ve vícenásobné soustavě přenáší hmotu.

## Chaotické běhy a špatně podmíněné začátky {#chaos}

V chaotických konfiguracích se dvě správné integrace po konečné době rozejdou: N-částicový
a sekulární běh uspořádání z obr. 4 Vokrouhlického (2016) asi po 4 kyr (článek uvádí totéž pro své
dva integrátory), SWINGERsim a MSE po třech cyklech ZLK čtyřhvězdy 2+2 a jednotlivé členy ansámblů
Pejchy et al. Podobná citlivost se projevuje u pravidelných drah začínajících těsně vedle
separatrisy: soustava ZLK začínající s $e_0$ = 0,01 se po dvou cyklech liší v časování svých cyklů
mezi N-částicovým a sekulárním backendem o 10 % (0,5 % pro konstantní hmotnosti daleko od
separatrisy). Ani jedno není chybou. Co by měl uživatel udělat: srovnávat chaotické soustavy
statisticky (podíly ansámblů, jako v [kapitole 2](#pejcha2013), kde oba backendy souhlasí v rámci
chyb) a vyhýbat se srovnávání jednotlivých členů nebo okamžiků za Ljapunovovou dobou.

## Sekulární backend blízko meze stability {#secular-limit}

Sekulární backendy středují přes dráhy, a proto předpokládají stabilní hierarchii s dobře
oddělenými periodami. Blízko meze stability to selhává: v ASASSN-V J075203 (poměr period 6–9) je
N-částicová apsidální perioda A asi poloviční oproti sekulární a v CzeV4315 se polovina
sekulárních řešení 2+2 N-částicově ukáže jako chaotická nebo se srazí
([kapitola 6](#squadra-systems)). Sekulární běh nestabilitu vidět nemůže; ve svém souhrnu uvádí
poměry MA01 a V22 a oktupólový parametr a parametr „jednoduchého středování“ každé vnořené trojice
a označuje mírné hierarchie („secular estimates approximate“). Co by měl uživatel udělat: použít
sekulární backend k předvýběru a pak členy, na nichž záleží, ověřit N-částicovým backendem
a MEGNO.

## Co znamená váha ansámblu soustředěná v jednom členu {#weights}

Evidence modelu je střední váhou jeho členů a její spolehlivost se měří efektivním počtem členů
$n_\mathrm{eff}$ ([kapitola 6](#squadra-systems)). Jsou-li pozorování přesná a apriorní rozdělení
široké, může téměř celá váha připadnout jedinému členu. ASASSN-V J075203 s izotropním apriorním
rozdělením dává $n_\mathrm{eff}$ = 1,0: jeden člen z 20 000 nese vše, takže „evidence“ je váhou
jediného šťastného losu a o modelu nic neříká. Čtyřhvězda 2+2 CzeV4315 podmíněná zákryty A má
$n_\mathrm{eff}$ = 4,8 a srovnání modelů se stejným apriorním rozdělením 119 d mělo
$n_\mathrm{eff}$ = 1,6, a je proto jen zašuměným odhadem. SWINGERsim vypisuje $n_\mathrm{eff}$
vedle každé evidence a každého váženého podílu. Co by měl uživatel udělat: věřit evidenci jen při
$n_\mathrm{eff}$ v řádu desítek či více; jinak spustit více členů, zúžit apriorní rozdělení na
oblast, která data prokládá (jako to udělal cílený ansámbl J075203, $n_\mathrm{eff}$ = 53), nebo
použít emulátor či simulační inferenci.

## Veličiny, které nejsou definovány {#undefined}

Některé veličiny, na které se uživatel může ptát, pro danou dráhu neexistují. Nejzřetelnějším
případem je apsidální perioda téměř kruhové dráhy: v TIC 120362137 nutí třetí hvězda excentricitu
páru s periodou 3,3 d kmitat mezi 0,0001 a 0,016, oskulační $\omega$ skáče a „změřená“ apsidální
perioda nemá smysl (článek uvádí jen teoretickou hodnotu). Totéž platí pro Algol ($e\approx0{,}001$)
a pro dráhy s $e$ < 10⁻⁶. SWINGERsim takovou periodu uvádí jako „v běhu neměřitelnou“, místo aby
vypsal číslo. Podobně špatně definované jsou dlouhodobé extrémy silně interagujících členů, v nichž
rostou malé numerické rozdíly; jako výsledky se nepoužívají.

## Vývojové dráhy hvězd na hranicích svého rozsahu {#track-limits}

Když konfigurace požaduje od vývojových drah něco, co nemohou dát, nastávají tři situace. (i) U velmi
hmotných hvězd se silným větrem může současné hmotnosti hvězdy v daném stáří dosáhnout dvě
vývojové dráhy s různou počáteční hmotností (například 79 M☉ ve stáří 2,5 Myr dráhy asi 108
a 139 M☉ v ansámblech QZ Car); SWINGERsim vezme nejméně vyvinutou hvězdu a vypíše varování.
(ii) Hmotnost, která v požadovaném stáří na zvolené mřížce neexistuje (Aa1 soustavy QZ Car nad
73,3 M☉ na rotujících drahách ve stáří 2,5 Myr), nelze spustit; člen končí osudem *chyba*
s uvedením důvodu. (iii) Hvězda, která už na začátku přeplňuje svůj Rocheův lalok, je odmítnuta
s uvedením důvodu. Odmítnuté členy jsou informací, nikoli selháním: v QZ Car omezují hmotnost Ac1
shora ([kapitola 5](#qzcar-evolution)). Chyby se v souhrnu každého ansámblu počítají zvlášť
a uživatel je může po změně spustit znovu pomocí `--resume --retry-errors`.

## Tabulky modelových atmosfér {#atmospheres}

Tabulky okrajového ztemnění Clareta & Bloemenové (2011) končí na 3500 K a pokrývají $\log g$ od 0
do 5. Dlouhoperiodická proměnná s teplotou 2867 K nebo veleobr s $\log g$ < 0 použije nejbližší
uzel, s varováním, které uvádí hvězdu a hodnoty. Pro relativní světlo takové hvězdy ve vícenásobné
soustavě to stačí; pro profil jejího vlastního zákrytu ne a tabulky pro chladnější hvězdy jsou mimo
rozsah SWINGERsimu.

## Zachovávané veličiny jako diagnostika {#conservation}

Každý běh uvádí své relativní chyby energie a momentu hybnosti. S relativitou je zachovávanou
veličinou hamiltonián 1PN, nikoli newtonovská energie; SWINGERsim měří chybu vůči hamiltoniánu
efektu `gr_full` z REBOUNDx, jinak by správný relativistický běh vykazoval falešný drift.
S disipací nebo ztrátou hmoty se energie fyzikálně nezachovává a diagnostikou je pak adiabatický
invariant ($aM$ pro vítr, ověřeno na 4·10⁻⁹ v [kapitole 3](#mass-loss)). Další diagnostikou
konzervativní dynamiky je symetrie v čase: selhání dopředu a dozadu v čase souhlasí v ansámblech
CzeV4315 na 0,02 % ([kapitola 6](#czev4315)).

## Problémy v publikovaných referencích {#published-problems}

Sledování každého rozdílu zpět k jeho článku odhalilo problémy v samotných referencích. Uvádíme je
zde, protože na ně narazí jakýkoli program srovnávaný s těmito články.

Tabulka 46 je uvádí; každý řádek jmenuje referenci, co je v ní chybné nebo nejednoznačné a jak se
s tím srovnávací test vypořádává.

**Tabulka 46.** Problémy v publikovaných referencích zjištěné při jejich reprodukci. Sloupec 2
popisuje problém, sloupec 3 co používá srovnávací test SWINGERsimu.

| Reference | Problém | Co používá srovnávací test |
|---|---|---|
| Naoz et al. (2013), obr. 4 | popisek (úhly 6,75° a 63,25°) a uvedené $a_2$ = 50 au implikují různé poměry momentů hybnosti (7,60 a 6,79); publikovaná křivka patří k 7,60 | $a_2$ = 39,9 au, což dává 7,60 |
| Veras et al. (2011), rov. 17 | jmenovatel $1-e_0\cos f$ neplyne z jejich rov. 4 | tvar odvozený z rov. 4, $1+e_0\cos f$ |
| Antognini & Thompson (2016), tabulka 3 | uvedené normování $\pi(a_1^2+a_2^2)$; čísla odpovídají $\pi a^2$ | $\pi a^2$ |
| Zasche et al. (2023), V994 Her C | uvedený relativistický podíl 6 %; jejich tabulka 3 dává 8,1 % | EXPLAINED s 8,1 % |
| Hurley et al. (2002), příklad Algolu | publikovaná historie by vyžadovala, aby dráha v Hertzsprungově mezeře získala moment hybnosti | EXPLAINED; konzervativní fáze obra je samostatným testem |
| Portegies Zwart & Meinen (1993) | poloměry záchytu citované z Raye et al. (1987) patří k jiným konvencím | EXPLAINED; referencí je McMillan et al. (1987) |
| Vynatheya et al. (2023), MA01 pro 3+1 | uvedené skóre 0,56 neplyne z jejich vlastní tabulky 4 (0,69) | EXPLAINED; množina losovaná blízko hranice |
| Pejcha et al. (2013), Vokrouhlický (2016) | podíly rostou s délkou běhu, která není uvedena dost přesně | N-částicově PASS na úrovni 3σ; sekulárně EXPLAINED |
| Kemp et al. (2024), KIC 4150611 | dynamický poměr hmotností odporuje poměrům velikostí a period ze stejného článku | hmotnosti konzistentní s poměry |
| Pribulla et al. (2020), VW LMi | uvedená apsidální perioda 78,6 roku je proložením tří let dat; vlastní modely článku dávají 95,5–96,2 roku | srovnání s N-částicovým modelem článku |
| MSE (Hamers et al. 2021) | argument periastra se při sklonu přesně 0° ztrácí | sklony 0,5–2° |

Ponaučení je prosté, ale důležité: publikované číslo je měřením něčeho a dříve, než se rozdíl
označí za chybu jedné či druhé strany, je třeba zjistit, co přesně bylo spočteno, z jakého vstupu,
s jakou definicí a za jakou dobu.
