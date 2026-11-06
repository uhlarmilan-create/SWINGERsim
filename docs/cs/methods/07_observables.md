# Pozorovatelné veličiny {#observables}

Nikdo nepozoruje vícenásobnou hvězdu milion let. Astronomové mají jen snímky: jasnost soustavy
v průběhu nocí až let, okamžiky, kdy nastávají její zákryty, rychlosti ze spekter. Tato kapitola
popisuje, jak SWINGERsim převádí dráhy a hvězdy z předchozích kapitol na tato měření: zákryty a jejich
detekovatelnost, světelné křivky, radiální rychlosti, okamžiky minim a diagramy O−C, pulzující hvězdy
jako hodiny a prachové disky jako zakrývající tělesa. Všechny pozorovatelné veličiny se počítají
z N-částicového backendu, protože potřebují polohu každé hvězdy na její dráze.

## Geometrie zákrytů {#eclipse-geometry}

Dvojhvězda zakrývá, když při pohledu podél zorného paprsku přechází kotouč jedné hvězdy přes kotouč
druhé. Pro kulové hvězdy o poloměrech $R_1$, $R_2$ na kruhové dráze se vzdáleností složek $a$ je
podmínkou

$$ |\cos i| < \frac{R_1 + R_2}{a} , $$

kde $i$ je sklon k obloze ([Hierarchie a dráhy](#orbital-elements)). U excentrické dráhy se $a$
nahradí vzdáleností složek v konjunkci. Pro izotropní orientaci je pravděpodobnost, že dvojhvězda
zakrývá, rovna $(R_1+R_2)/a$; pro obě dvojhvězdy V994 Her to dává 37 % a 40 %, jak zjistil
Vokrouhlický (2016) (viz *Testování, část I*).

![**Obrázek 12.** Zákrytová dvojhvězda. Vlevo a uprostřed: dvojhvězda na obloze (pozorovatel se dívá
podél $z$); šedá elipsa je promítnutá dráha modré hvězdy kolem oranžové. Při $i = 88°$ (vlevo) je
promítnutá dráha téměř úsečkou a kotouče se v konjunkci překrývají; při $i = 60°$ (uprostřed) se
míjejí bez překrytí. Vpravo: světelná křivka zákrytového případu, jasnost (normovaná na 1 mimo
zákryt) v závislosti na oběžné fázi. Hlubší minimum ve fázi 0 je primární zákryt, kdy je skryta
teplejší hvězda; mělčí ve fázi 0,5 sekundární zákryt.](../../figures/concepts/eclipse_geometry.png)

Obrázek 12 ukazuje oba případy. Hloubka každého minima je podíl celkového světla, který je skryt;
závisí na překrytí kotoučů a na povrchové jasnosti skryté hvězdy: skrytí teplejší hvězdy stojí více
světla, a proto je primární minimum hlubší.

Během běhu udává SWINGERsim pro každou dvojhvězdu ze dvou hvězd podíl času, po který platí výše
uvedená podmínka (kulové hvězdy, tvar pro kruhovou dráhu, jako u Vokrouhlického 2016, odd. 3.3),
s poloměry, které se mění s vývojem hvězd; pro soustavy 2+2 udává také podíl času, kdy zakrývají obě
dvojhvězdy současně – definiční vlastnost dvojitě zákrytových soustav. Ansámbly uvádějí, zda
dvojhvězda v $t = 0$ zakrývá a kdy její zákryty začínají nebo končí.

### Geometrické a detekovatelné zákryty {#detectable-eclipses}

Geometrická podmínka říká jen, že se kotouče překrývají. Zda pozorovatel zákryt *zaznamená*, závisí
na jeho hloubce v celkovém světle. Ve vícenásobné hvězdě světlo všech ostatních hvězd zákryty
dvojhvězdy ředí. Ztratila-li by dvojhvězda sama v minimu podíl $d$ svého vlastního světla a přispívají-li
ostatní hvězdy podílem $l_3$ celkového světla (*třetí světlo*), je pozorovaná hloubka v magnitudách

$$ \Delta m = -2{,}5\,\log_{10}\left[1 - (1 - l_3)\,d\right] . $$

Tečný zákryt (malé $d$) v jasné vícenásobné hvězdě ($l_3$ blízko 1) proto může být hluboko pod
fotometrickou přesností, přestože geometrická podmínka platí.

![**Obrázek 13.** Mez zákrytů a detekovatelnost (nový obrázek pro tento svazek, spočtený z obou výše
uvedených rovnic). (a) Největší odklon od pohledu z boku,
$|90° - i|_\mathrm{max} = \arcsin[(R_1 + R_2)/a]$, při kterém kruhová dvojhvězda ještě zakrývá,
v závislosti na její periodě $P$ ve dnech (obě osy logaritmické). Modře: dvě hvězdy ZAMS 1 M☉
o poloměru 0,89 R☉; oranžově: dvě hvězdy 2 M☉ o poloměru 1,70 R☉; zeleně: dva podobři 1 M☉
o poloměru 5 R☉ (zobrazeno jen tam, kde jsou hvězdy dobře oddělené). (b) Pozorovaná hloubka
$\Delta m$ v magnitudách v závislosti na vlastní poměrné hloubce $d$ dvojhvězdy (obě logaritmicky)
pro třetí světlo $l_3$ = 0 (modře), 0,5 (oranžově), 0,9 (zeleně) a 0,99 (žlutě). Čárkovaná čára
vyznačuje příkladnou mez detekce 1 mmag.](../../figures/methods/eclipse_detectability.png)

Obrázek 13 obě otázky odděluje. Panel (a): dvojice hvězd podobných Slunci s periodou 1 d zakrývá do
odklonu asi 20° od pohledu z boku, při 10 d do asi 4°, při 100 d do méně než 1°; větší hvězdy okno
rozšiřují úměrně $R_1 + R_2$. Proto jsou dlouhoperiodické dvojhvězdy zřídka zákrytové a proto
precedující dvojhvězda své zákryty získává a ztrácí ([Hierarchie a dráhy](#nodal-precession)).
Panel (b): pro malé $d$ je pozorovaná hloubka $\Delta m \approx 1{,}086\,(1 - l_3)\,d$, takže každá čára
je modrou čarou posunutou dolů o faktor $1 - l_3$. Vlastní hloubka 1 % se jeví jako 11 mmag pro
samotnou dvojhvězdu, 5 mmag s $l_3$ = 0,5, 1,1 mmag s $l_3$ = 0,9 a 0,11 mmag s $l_3$ = 0,99 –
poslední je pod příkladnou mezí detekce. Když se SWINGERsim používá k předpovědi, kdy se zákryty stanou
viditelnými nebo zmizí, rozhoduje syntetická světelná křivka (níže), nikoli geometrická podmínka.

## Světelné křivky {#light-curves}

`swingersim lightcurve` počítá světlo, které pozorovatel zaznamená, z N-částicových poloh všech
hvězd. K dispozici jsou tři úrovně věrnosti (tabulka 11).

**Tabulka 11.** Modely světelných křivek. *Model* pojmenovává úroveň; *Hvězdy* tvary; *Fyzika* co je
zahrnuto; *Použití* kde je model vhodný. Všechny úrovně sdílejí polohy, dobu letu světla a normování
(1 = všechny hvězdy mimo zákryt).

| Model | Hvězdy | Fyzika | Použití |
|---|---|---|---|
| koule (výchozí) | koule s okrajovým ztemněním | přesné zákryty libovolného počtu hvězd, analytické členy elipsoidální proměnnosti, odrazu a dopplerovského zesílení | oddělené soustavy, $R/R_L$ hluboko pod 1 |
| sítě | Rocheovy plochy, rotující hvězdy, kontaktní dvojhvězdy; všechny hvězdy jako sítě | gravitační a okrajové ztemnění pro každý plošný element, zákryt podél každého paprsku | těsné, polodotykové a kontaktní dvojhvězdy, rychle rotující hvězdy |
| hybrid s PHOEBE | PHOEBE 2 pro každou vnitřní dvojhvězdu | Rocheova geometrie, modelové atmosféry, ozáření v PHOEBE; vzájemné zákryty různých dvojhvězd jako koule | vybraná řešení, křížové kontroly, trénovací množiny |

### Koule: intenzita, okrajové ztemnění a fotometrická pásma {#limb-darkening}

Povrchová jasnost hvězdy ve fotometrickém pásmu je jasnost černého tělesa integrovaná přes pásmo
s vahami pro počítání fotonů, $w = \int B_\lambda(T_\mathrm{eff})\,S(\lambda)\,\lambda\,d\lambda$, kde
$B_\lambda$ je Planckova funkce a $S(\lambda)$ propustnost pásma. Dostupná pásma jsou U, B, V, R, I
Johnsonova–Cousinsova systému (Bessell 1990), TESS (Ricker et al. 2015) a Gaia G, BP, RP (Riello et
al. 2021). Modelové atmosféry mění poměr povrchových jasností hvězd různých teplot o několik procent;
hybrid s PHOEBE používá modelové atmosféry tam, kde na tom záleží.

Intenzita napříč kotoučem se řídí čtyřkoeficientovým zákonem okrajového ztemnění podle Clareta
(2000),

$$ \frac{I(\mu)}{I(1)} = 1 - \sum_{k=1}^{4} a_k\left(1-\mu^{k/2}\right), $$

kde $\mu$ je kosinus úhlu mezi zorným paprskem a normálou povrchu a $a_k$ jsou koeficienty; lineární
a kvadratický zákon jsou zvláštními případy a koeficienty lze zadat pro každou hvězdu. Pod 15 000 K
pocházejí koeficienty a exponenty gravitačního ztemnění od Clareta & Bloemen (2011), interpolované
v $T_\mathrm{eff}$, $\log g$ a [M/H]; pro TESS a Gaia se použije nejbližší pásmo Johnsonova–Cousinsova
systému z těchto tabulek (TESS a RP: I; G: R; BP: V). Od 15 000 do 55 000 K pocházejí koeficienty
z non-LTE atmosfér TLUSTY (OSTAR2002, BSTAR2006) podle Reeva & Howartha (2016); tam, kde platí obě
mřížky, souhlasí intenzity středované přes kotouč na 0,2–0,7 %. Za hranicí mřížky gravitačního
ztemnění Clareta & Bloemen (hvězdy O nad zhruba 39 000 K) se exponent řídí von Zeipelovým zákonem
v daném pásmu, $y = \tfrac14\,x e^x/(e^x - 1)$, $x = hc/(\lambda k T_\mathrm{eff})$. Hvězdy mimo
tabulky (například chladné hvězdy AGB) přebírají nejbližší tabelovaný uzel s varováním.

**Zákryty.** Každou hvězdu zakrývají všechny hvězdy bližší k pozorovateli, jejichž kotouče se s ní
překrývají. Skrytý tok je $\int_0^R I(r)\,\alpha(r)\,2\pi r\,dr$, kde $\alpha(r)$ je podíl kružnice
o poloměru $r$ kolem středu hvězdy pokrytý sjednocením zakrývajících kotoučů; radiální integrál se
počítá Gaussovou–Legendreovou kvadraturou mezi poloměry, v nichž kružnici protíná okraj některého
kotouče (jako v ellc, Maxted 2016). To je přesné pro libovolný počet koulí, takže vzájemné zákryty
obou dvojhvězd soustavy 2+2 i trojnásobná překrytí nevyžadují zvláštní zacházení. Souhlasí
s analytickým modelem Mandela & Agola (2002), jak je implementován v batman (Kreidberg 2015), na
$10^{-8}$.

**Efekty blízkosti** hvězd dvojhvězd ze dvou hvězd, v aproximaci malých deformací: elipsoidální
proměnnost $-\alpha_\mathrm{ell}\,q\,(R/d)^3(2\cos^2\psi - \sin^2 i)$ s
$\alpha_\mathrm{ell} = 0{,}15(15+u)(1+\tau)/(3-u)$ (Morris 1985; Morris & Naftilan 1993), kde $q$ je
poměr hmotností, $d$ vzdálenost složek, $\psi$ úhel mezi spojnicí středů a zorným paprskem, $u$
lineární koeficient okrajového ztemnění a $\tau$ koeficient gravitačního ztemnění; dopplerovské
zesílení (beaming) $-(3-\alpha_\nu)\,v_z/c$ se spektrálním indexem $\alpha_\nu$ (Loeb & Gaudi 2003);
a odraz světla průvodce lambertovskou koulí s Bondovým albedem 0,5 pro konvektivní a 1 pro zářivé
obálky (Rucinski 1969). Elipsoidální člen podhodnocuje proměnnost o 3, 7, 11 a 17 % při $R/a$ = 0,10,
0,15, 0,20 a 0,25 ve srovnání s Rocheovým modelem v ellc; těsnější dvojice potřebují sítě.

**Expozice a šum.** Konečné expozice se modelují středováním modelu přes dílčí expozice (doba
expozice a počet dílčích bodů na expozici); lze přidat gaussovský bílý šum.

![**Obrázek 14.** Syntetická světelná křivka soustavy 2+2 spočtená s koulemi. Nahoře: celkový tok
(normovaný na 1 mimo zákryty) v závislosti na čase ve dnech. Dole: tok každé hvězdy dělený jejím
vlastním nezakrytým tokem, pro m1 (tmavě modře), m2 (světle modře), m3 (oranžově) a m4 (světle
oranžově). Zákryty dvojhvězdy A (m1 a m2 skryté střídavě) se v celkovém světle jeví zředěné světlem
B; malé pomalé vlnky jsou elipsoidální proměnnost a odraz.](../../figures/concepts/light_curve.png)

Obrázek 14 ukazuje, jak vzniká signál vícenásobné hvězdy. Spodní panel hvězdy odděluje: m2 ztratí při
svém zákrytu asi 12 % vlastního světla, v celkovém světle (nahoře) se z toho však stane pokles o 3 %,
protože ostatní tři hvězdy svítí dál. Toto zředění je třetí světlo z předchozího oddílu a je důvodem,
proč je k interpretaci hloubek zákrytů potřeba model celé soustavy, nikoli jedné dvojhvězdy.

### Doba letu světla {#light-travel-time}

Každou hvězdu vidíme v jejím vlastním *okamžiku vyslání*: světlo pozorované v čase $t$ opustilo hvězdu
v hloubce $z$ (kladné směrem od pozorovatele, barycentrum v $z = 0$) v čase $t_e = t - z/c$. Každá
dvojhvězda ze dvou hvězd se integruje zvlášť do okamžiku vyslání svého těžiště,
$t_e = t - z_\mathrm{com}(t_e)/c$ (řešeno iteračně), a jen odchylky obou jejích hvězd od těžiště se
opraví v prvním řádu, $-\mathbf{v}\,\Delta z/c$ (Rømerovo zpoždění uvnitř dvojhvězdy); ostatní hvězdy
používají $\mathbf{x}(t) - \mathbf{v}\,z/c$. Oprava celého zpoždění v prvním řádu nestačí, když je
$z/c$ znatelným zlomkem periody dvojhvězdy: u QZ Car ($z/c \approx 0{,}2$ d, $P$ = 6 d) by posunula
zákryty asi o 200 s.

### Deformované hvězdy {#distorted-stars}

Hvězdám lze zadat Rocheův tvar (`stars.<m>.shape: roche`, obě hvězdy kruhové synchronizované
dvojhvězdy), zploštělý tvar tuhé rotující hvězdy (`shape: rotating`), nebo může být dráha kontaktní
dvojhvězdou s faktorem přeplnění $f$ (`orbits.<X>.contact_fillout`). Model dvojhvězdy používá
Rocheův potenciál (Kopal 1959; Wilson & Devinney 1971)

$$ \Omega(x,y,z) = \frac{1}{r_1} + q\left(\frac{1}{r_2} - x\right) + \frac{1+q}{2}\,(x^2+y^2) , $$

v jednotkách vzdálenosti složek, s hvězdou 1 v počátku, hvězdou 2 v bodě $(1, 0, 0)$, $r_1$, $r_2$
vzdálenostmi k oběma hvězdám, $q = M_2/M_1$ a osou $z$ podél normály dráhy. Vnitřní kritická plocha
$\Omega_\mathrm{in} = \Omega(\mathrm{L1})$ ohraničuje Rocheovy laloky, vnější $\Omega_\mathrm{out}$
(procházející L2 nebo L3) společnou obálku; faktor přeplnění je
$f = (\Omega_\mathrm{in}-\Omega)/(\Omega_\mathrm{in}-\Omega_\mathrm{out})$, záporný pro oddělené
a mezi 0 a 1 pro kontaktní dvojhvězdy (Lucy & Wilson 1979; Rucinski 1973). Rotující hvězda se řídí
Rocheovým modelem $1/r + (4/27)\,w^2 r^2\sin^2\theta = 1$ v jednotkách svého polárního poloměru,
s $w = \Omega/\Omega_\mathrm{crit}$ a $\theta$ pólovou vzdáleností (Maeder 2009); při $w = 1$ dosáhne
její rovník 1,5násobku polárního poloměru. Zadaný poloměr je objemově ekvivalentní poloměr.

Každá plocha je vzorkována v 64 × 128 směrech. Element vyzařuje
$I = I_\mathrm{ref}\,(g/\bar g)^y\,L(\mu)$, s místním efektivním tíhovým zrychlením $g$, jeho
povrchovým průměrem $\bar g$, exponentem gravitačního ztemnění v pásmu $y$ (von Zeipel 1924; Claret
& Bloemen 2011) a okrajovým ztemněním $L(\mu)$; oba laloky kontaktní dvojhvězdy sdílejí jednu
efektivní teplotu. Element je skryt, když jeho normála míří od pozorovatele, když před ním leží koule
(s vyhlazením hran podle jeho promítnuté velikosti), nebo když jeho zorný paprsek prochází jiným
tělesem či druhým lalokem jeho vlastního Rocheova tělesa. Je-li některá hvězda sítí, jsou sítěmi
všechny hvězdy. Přesnost: $3\times10^{-4}$ toku oproti přesným zákrytům koulí; amplituda kontaktní
dvojhvězdy souhlasí s PHOEBE 2.5 na 0,004 mag (viz *Testování, část I*). Chyba koulí oproti Rocheovým
plochám strmě roste s vyplněním laloku: asi 0,5 ppt rms při $R/R_L$ = 0,3, 2 ppt při 0,5, 5 ppt
při 0,7, 9 ppt při 0,85 a 17 ppt při 0,95 a 15–45 ppt u kontaktních dvojhvězd. Rocheův model
předpokládá kruhovou synchronizovanou dráhu; excentrické deformované dvojhvězdy potřebují hybrid
s PHOEBE.

### Hybrid s PHOEBE {#phoebe}

S výpočetním modulem PHOEBE pochází světlo každé vnitřní dvojhvězdy ze dvou hvězd z PHOEBE 2 (Prša et
al. 2016; Horvat et al. 2018; Conroy et al. 2020): Rocheova geometrie, gravitační ztemnění, ozáření,
modelové atmosféry Castelliho–Kurucze nebo černá tělesa a doba letu světla uvnitř dvojhvězdy.
Vícenásobnou soustavu zajišťuje SWINGERsim:

1. v úsecích délky jedné periody dvojhvězdy dá N-částicová integrace oskulační elementy dvojhvězdy
   ve středu úseku a pohyb jejího těžiště podél zorného paprsku; perioda a okamžik průchodu
   periastrem se určí ze střední délky během jednoho oběhu, protože oskulační keplerovská perioda
   slapově deformované dvojhvězdy je chybná o minuty za oběh;
2. PHOEBE spočte dvojhvězdu v průběhu úseku v okamžicích vyslání;
3. hvězdy mimo dvojhvězdy ze dvou hvězd svítí svým zářivým výkonem v daném pásmu;
4. zákryty mezi hvězdami různých dvojhvězd se převezmou z modelu koulí a odečtou.

Konvence se převádějí ($z$ od pozorovatele; argument periastra v PHOEBE je argumentem primární
složky, `per0` = $\omega$ + 180°). Se stejnou fyzikou reprodukuje hybrid model koulí na několik
$10^{-4}$. PHOEBE potřebuje 0,1–0,2 s na bod a dvojhvězdu, asi $10^4$krát více než model koulí,
a je určen pro vybraná řešení. Rotace je synchronní, dopplerovské zesílení v PHOEBE 2.5 chybí a disky
se ignorují. Nativní zpracování vícenásobných soustav podle Brože et al. (2025) je alternativou,
musí-li veškeré světlo pocházet z jednoho kódu.

## Radiální rychlosti {#radial-velocities}

Radiální rychlost hvězdy $k$ je její rychlost podél zorného paprsku vůči barycentru plus systémová
rychlost $\gamma$: $\mathrm{RV}_k = \gamma + v_{z,k}$, kladná směrem od pozorovatele. Rychlost
naměřená v čase $t_\mathrm{obs}$ patří k okamžiku vyslání $t_e$, pro nějž
$t_e + z_k(t_e)/c = t_\mathrm{obs}$; SWINGERsim jej najde Newtonovými kroky, dokud nesoulad neklesne
pod 1 ms. Doba letu světla od hvězdy (na široké vnější dráze až desítky minut) je tedy zahrnuta, což je
důležité pro strmé křivky radiálních rychlostí těsných dvojhvězd (poloamplituda 100 km s⁻¹ na
třídenní dráze se mění o 0,15 km s⁻¹ za minutu).

**Střed světla a Rossiterův–McLaughlinův efekt.** Spektrograf měří dopplerovský posuv světla hvězdy –
tokem váženou střední rychlost jejího viditelného povrchu podél zorného paprsku –, nikoli rychlost
jejího těžiště. U kulových, pomalu rotujících hvězd se obě shodují. Slapově deformovaná hvězda je na
některých stranách jasnější a větší a její synchronní rotace přidává povrchovou rychlost, takže její
střed světla se pohybuje po jiné dráze (v kontaktních dvojhvězdách o desítky km s⁻¹; Wilson & Sofia
1976); během zákrytu chybí v průměru zakrytá část rotující hvězdy, což měřenou rychlost posune
(Rossiterův–McLaughlinův efekt; Rossiter 1924; McLaughlin 1924). S `fit.rv.centre: light` je modelová
rychlost tokem váženou rychlostí viditelných plošných elementů výše popsaných sítí podél zorného
paprsku – Rocheovy hvězdy korotují se svou dráhou, ostatní hvězdy rotují podle své rotace –
v pásmu `fit.rv.band`. Odchylka od rychlosti těžiště se přičte k přesné rychlosti těžiště, takže doba
letu světla i dynamika zůstávají ty z N-částicového modelu. Ozáření a dopplerovské zesílení nejsou ve
vahách zahrnuty. Pro kontaktní dvojhvězdu ($q$ = 0,42, přeplnění 0,3) dosahují odchylky ±34
a ±52 km s⁻¹ a souhlasí s PHOEBE 2 na 2 % své amplitudy; proložení takových rychlostí rychlostí
těžiště dá hmotnosti o několik procent příliš velké.

## Okamžiky minim a diagramy O−C {#times-of-minima}

Kdyby byla dvojhvězda sama, řídily by se její zákryty lineární *efemeridou* $T_n = T_0 + nP$ s číslem
cyklu $n$. Pozorovatelé porovnávají každý naměřený okamžik minima (O) s okamžikem vypočteným z takové
efemeridy (C) a vynášejí rozdíl O−C v závislosti na čase. Přímka znamená chybnou periodu; cokoli jiného
znamená, že dvojhvězdou něco pohybuje nebo mění její dráhu.

**Syntetické okamžiky minim.** `swingersim timings` integruje soustavu (N-částicově s nakonfigurovanou
fyzikou) a najde každý zákryt každé vnitřní dvojhvězdy: minima promítnuté vzdálenosti jejích dvou
hvězd, vzorkovaná 40krát za periodu dvojhvězdy a zpřesněná dvoučásticovým pohybem mezi vzorky.
Konjunkce je zákrytem, když je promítnutá vzdálenost menší než $R_1 + R_2$ (poměr se uvádí jako
zákrytový parametr). *Primární* zákryt je zákryt první hvězdy dvojhvězdy (m1 v A, m3 v B), tj. první
hvězda je ta vzdálenější. Pozorované okamžiky zahrnují dobu letu světla od těžiště dvojhvězdy,
$t_\mathrm{obs} = t + z_\mathrm{com}/c$; dynamická zpoždění plynou ze samotné integrace. Lineární
efemerida proložená primárními zákryty dává O−C; sekundární zákryty používají poloviční epochy.

### Efekt světelného času {#ltte}

Dvojhvězda obíhá kolem barycentra celé soustavy, takže je někdy dál od pozorovatele a někdy blíž.
Její zákryty přicházejí střídavě později a dříve: *efekt světelného času* (LTTE). Pro vnější dráhu
s periodou $P_3$, excentricitou $e_3$, argumentem periastra $\omega_3$ a pravou anomálií $\nu_3$ je
zpoždění (Irwin 1952)

$$\Delta_\mathrm{LTTE} = \frac{a_{12}\sin i_3}{c}\left[\frac{1-e_3^2}{1+e_3\cos\nu_3}\sin(\nu_3+\omega_3) + e_3\sin\omega_3\right],$$

kde $a_{12}\sin i_3$ je promítnutá velká poloosa dráhy dvojhvězdy kolem barycentra. Jeho poloamplituda
je $K = (a_{12}\sin i_3/c)\sqrt{1-e_3^2\cos^2\omega_3}$ a hmotnostní funkce
$f(m_3) = (a_{12}\sin i_3)^3/P_3^2$ (au, yr, M☉). N-částicové okamžiky minim obsahují tento člen
přesně, spolu se všemi ostatními efekty.

![**Obrázek 15.** Efekt světelného času. Vlevo: dvojhvězda (modrý bod) na své vnější dráze (zelená
elipsa) kolem barycentra (křížek); šipka míří k pozorovateli. Na vzdálené straně dráhy jsou její
zákryty vidět později, na bližší straně dříve. Vpravo: výsledné O−C v minutách v závislosti na čase ve
dnech, periodická křivka s vnější periodou a amplitudou rovnou době letu světla napříč promítnutou
dráhou dvojhvězdy (zde asi 10 min; ilustrativní).](../../figures/concepts/lite_oc.png)

Obrázek 15 spojuje geometrii s diagramem. Křivka O−C se opakuje s vnější periodou; její poloamplituda
je $a_{12}\sin i_3/c$ (asi 8,3 min na au) a pro excentrickou vnější dráhu se její tvar odchyluje od
sinusoidy. Ve čtyřhvězdě 2+2 obíhají obě dvojhvězdy kolem sebe, takže jejich křivky LTTE jsou
v protifázi s amplitudami v obráceném poměru jejich hmotností; v kompaktních soustavách s vnějšími
periodami do jednoho roku přidávají dynamické poruchy periody složku, která v protifázi není a může
být větší než LTTE.

### Apsidální stáčení v diagramu O−C {#apsidal-motion-in-o-c}

V excentrické dvojhvězdě závisí doba mezi primárním a sekundárním minimem na orientaci dráhy. Když se
periastrum stáčí (apsidální stáčení), posouvají se primární a sekundární minima v diagramu O−C
opačnými směry, s apsidální periodou $U$.

![**Obrázek 16.** Apsidální stáčení (ilustrativní). Vlevo: excentrická dráha (modré elipsy, odstíny
pro po sobě jdoucí okamžiky) kolem oranžové hvězdy; body vyznačují periastrum, které se pomalu stáčí.
Vpravo: O−C v minutách v závislosti na čase v letech pro primární (modře) a sekundární (oranžově)
minima: dvě křivky v protifázi s apsidální periodou, zde asi 25 let.](../../figures/concepts/apsidal.png)

Obrázek 16 ukazuje znak, podle něhož se měří $\dot\omega$. Amplituda obou křivek roste
s excentricitou a jejich protifáze odlišuje apsidální stáčení od LTTE, které posouvá obě minima
společně. Analytický tvar použitý v proložení v kroku 1 je uveden v kapitole
[Prokládání](#stage-1).

### Fenomenologická minima {#phenomenological-minima}

Pozorovatelé zákryt neprokládají fyzikálním modelem; měří jeho okamžik minima, hloubku a šířku
fenomenologickým profilem. `swingersim minima` redukuje syntetickou světelnou křivku stejně,
s profilem Mikuláška (2015, rov. 13),

$$\Delta m(t) = A\,(1 + C x^2)\left\{1 - \left[1 - e^{1 - \cosh x}\right]^{\Gamma}\right\},
\qquad x = \frac{t - t_\mathrm{min}}{D},$$

s centrální hloubkou $A$ (mag), šířkou $D$, špičatostí $\Gamma$ a korekcí $C$, na kvadratické základní
úrovni (světlo ostatních hvězd, elipsoidální proměnnost a odraz, které jsou do druhého řádu symetrické
kolem minima, a pomalý drift). Zákryty se lokalizují pomocí syntetických okamžiků minim; okno kolem
každého je 2,5násobkem jeho trvání (nejvýše polovina periody). *Tvar* ($A$, $D$, $\Gamma$, $C$) se
prokládá pro každou dvojhvězdu a druh (primární, sekundární) všemi jejími zákryty společně, volitelně
pro každou sezónu (což sleduje pomalou změnu hloubky a šířky při precesi dráhy); potom se okamžik
každého zákrytu proloží tímto tvarem s volnou hloubkou, šířkou a základní úrovní a jeho chyba se určí
z kovarianční matice škálované rezidui. Standardně se každá dvojhvězda měří na svém vlastním světle
plus středním světle ostatních (ideální rozdělení světla); s volbou smíšeného světla se použije
celkové světlo tak, jak je pozorováno, zákryty jiných dvojhvězd v okně se zamaskují a jejich proměnnost
z blízkosti posouvá minima (v příkladu podobném CzeV4315 vzroste rozptyl minim B z 0,2 s na 8 s).

Rozdíl mezi minimem a barycentrickou konjunkcí je fyzikální. Doba letu světla uvnitř dvojhvězdy
zpožďuje zákryt vzdálenější hvězdy: sekundární minimum kruhové dráhy přichází
$(P/\pi c)(K_2 - K_1) = 2a\sin i\,(m_1 - m_2)/(Mc)$ po polovině periody (Kaplan 2010; $K_1$, $K_2$ jsou
poloamplitudy radiálních rychlostí), což je reprodukováno na 0,4 %. Pro křivku bez šumu prokládá
Mikuláškův profil oddělené zákryty asi na 1 % jejich hloubky.

## Pulzující hvězdy jako hodiny {#pulsating-stars}

Mnoho složek vícenásobných soustav pulzuje: hvězdy typu δ Scuti a γ Doradus kolem 1,5–2,5 M☉, hvězdy
typu β Cephei a pomalu pulzující hvězdy B (SPB), hvězdy RR Lyrae, cefeidy a dlouhoperiodické
proměnné na AGB. S `stars.<m>.pulsation` se vlastní tok hvězdy mění jako

$$\Delta m(t) = \sum_k A_k \sin\left[2\pi\left(\frac{t_e}{P_k} + \phi_k\right)\right],$$

s amplitudami $A_k$, periodami $P_k$ a fázemi $\phi_k$ módů, vyhodnocenými v *okamžiku vyslání* $t_e$
hvězdy. Obíhající pulzující hvězda je proto hodinami: její fáze se zpožďují o dobu letu světla napříč
jejími drahami, fázová modulace podle Murphyho et al. (2014). Pro kruhovou dráhu jsou zpoždění
$\tau(t) = (a_1\sin i/c)\sin(2\pi t/P_\mathrm{orb} + \varphi)$ ($a_1$ je velká poloosa dráhy pulzující
hvězdy kolem barycentra); ve vícenásobné soustavě se zpoždění vnitřní dráhy sčítají se zpožděními
vnější. SWINGERsim reprodukuje zpoždění naměřená ze světelné křivky vzorkované jako data Kepleru na
$3\times10^{-7}$.

![**Obrázek 17.** Pulzující složka jako hodiny (ilustrativní). Nahoře: jasnost pulzující hvězdy
(libovolné jednotky) v závislosti na čase v hodinách, na bližší straně její dráhy (šedě) a o půl
oběhu později (oranžově): tytéž pulzy přicházejí později. Dole: zpoždění fáze pulzace v minutách
v závislosti na čase v letech (oranžová čára) a O−C zákrytů téže dvojhvězdy (modré body): obě sledují
tutéž křivku doby letu světla na dráze.](../../figures/concepts/pulsator_clock.png)

Obrázek 17 ukazuje, proč pulzující hvězda v zákrytové dvojhvězdě nese dvoje hodiny, které se musejí
shodovat: zákryty vzorkují křivku doby letu světla jednou nebo dvakrát za periodu dvojhvězdy, pulzace
spojitě.

Periody plynou z modelu hvězdy (současná hmotnost, poloměr, zářivý výkon, $T_\mathrm{eff}$), pokud
není zadáno `period_d` (tabulka 12).

**Tabulka 12.** Typy pulzací a zdroj jejich period. *Typ* je hodnota klíče `pulsation.type`; *Vztah pro
periodu* použitý vzorec; *Zdroj* reference a její rozsah.

| Typ | Vztah pro periodu | Zdroj |
|---|---|---|
| `lpv` | nelineární vztah perioda–hmotnost–poloměr pro základní mód, lomený mocninný zákon s nasycením | Trabucchi et al. (2021), rov. 1–2; AGB bohaté na kyslík, 0,6–7 M☉ |
| `delta_scuti` | $Q = P\sqrt{\bar\rho/\bar\rho_\odot}$ = 0,0329 d (F), 0,0251, 0,0203, 0,0172 d (1. až 3. overton) | Breger (1979); v souladu s Poro et al. (2024) |
| `rr_lyrae` | $\log P_0 = 11{,}627 + 0{,}823\log L - 3{,}506\log T_\mathrm{eff} - 0{,}582\log M$; $P_1/P_0 = 0{,}744$ | Bono et al. (1996), jak je cituje Alcock et al. (1997) |
| `cepheid` | $\log P = a + b\log T_\mathrm{eff} + c\log M + d\log L$ (F), první overton zvlášť | Bono, Castellani & Marconi (2000), Z = 0,004–0,02 |
| `beta_cep` | radiální základní mód, $Q$ = 0,033 d | Handler (2013) |
| `spb`, `gamma_dor` | řada gravitačních módů $P_0 + k\,\Delta P$; $P_0$, $\Delta P$ zadává uživatel | Wu, Li & Deng (2018) |
| `manual` | zadané `period_d` nebo $Q$ | – |

Každý typ varuje, leží-li $T_\mathrm{eff}$ hvězdy mimo její třídu (kontrola věrohodnosti).
`swingersim pulsation` udává také $dP/dt$ z vyvíjejícího se modelu hvězdy, tedy vývojovou změnu
periody, kterou měří diagramy O−C pulzujících hvězd. Hvězdy MIST uprostřed pásu nestability δ Scuti
leží na pozorovaném vztahu perioda–zářivý výkon podle Pora et al. (2024) do 0,1 mag a cefeidy MIST na
svých modrých smyčkách souhlasí se vztahem perioda–poloměr podle Andersona et al. (2016) na 0,01 dex.
Amplitudy jsou pro každé pásmo pevné, tvary módů jsou sinusoidy (bez kombinačních frekvencí
a nesinusových tvarů RR Lyrae a mirid) a pulzace nemění poloměr použitý pro zákryty.

## Disky jako zakrývající tělesa {#discs}

Některé zákryty způsobuje temný prachový disk: dvouleté zákryty ε Aur diskem kolem jejího průvodce,
zákryty EE Cep, jejichž hloubka se mění z cyklu na cyklus, jak její disk preceduje, nebo KH 15D, mladá
dvojhvězda mizející za okrajem precedujícího cirkumbinárního prstence. Disk je ve SWINGERsimu deskou:
válcem o vnějším poloměru $r_\mathrm{out}$, s volitelným otvorem $r_\mathrm{in}$ a celkovou tloušťkou
$H$, se středem ve hvězdě nebo v těžišti dráhy. Pro každý zorný paprsek z bodu hvězdy k pozorovateli
se přesně spočte délka $L$ uvnitř desky a propustnost je $e^{-\tau L/H}$, kde $\tau$ je optická
tloušťka kolmo k disku (bez $\tau$ je disk neprůhledný). Disk viděný z boku je tmavým pásem výšky $H$,
skloněný disk elipsou s okrajem. Hvězdný kotouč je vzorkován ve 4000 bodech stejné plochy vážených
okrajovým ztemněním; ostrý okraj přecházející přes hvězdu je rozlišen asi na $10^{-3}$ jejího toku.
Disk může být pevný na obloze, skloněný vůči dráze, precedovat zadanou rychlostí nebo precedovat jako
tuhé těleso pod vlivem momentu síly své dvojhvězdy (`dynamics: rigid`; Zanazzi & Lai 2018, rov. 3;
Larwood et al. 1996 pro cirkumstelární disky), s polárním srovnáním kolem excentrických dvojhvězd
(Farago & Laskar 2010; Martin & Lubow 2017). Disky mohou svítit jako ploché disky přezařující světlo
hvězdy, $T = (2/3\pi)^{1/4}(R_*/r)^{3/4}T_*$ (Chiang & Goldreich 1997), a rozptylovat světlo hvězdy
s albedem a Henyeyovou–Greensteinovou fázovou funkcí (Henyey & Greenstein 1941). Nemodeluje se:
hmotnost a gravitace disku působící na dráhy, dopředný rozptyl přes opticky tenké části, zborcení
disku.

![**Obrázek 18.** Světelné křivky soustav se zákryty disky spočtené SWINGERsimem. Vlevo: ε Aur
2009–2011, zeslabení v magnitudách v závislosti na roce; čárkované svislé čáry: pozorované kontakty I
a IV; tečkovaná vodorovná čára: pozorovaná hloubka. Uprostřed: EE Cep, hloubka každého zákrytu
v magnitudách v závislosti na čase v letech během jedné precese jejího disku (ilustrativní model;
stínováno: pozorovaný rozsah 0,5–2 mag). Vpravo: KH 15D, zeslabení v magnitudách v závislosti na čase
v letech, jak okraj jejího neprůhledného prstence přejíždí přes dráhu dvojhvězdy (ilustrativní model;
fáze, v nichž je dvojhvězda zcela skryta, leží mimo stupnici).](../../figures/discs.png)

Obrázek 18 ukazuje tři režimy. ε Aur: neprůhledná deska přecházející před hvězdou dává zákryt
s plochým dnem, jehož kontakty a hloubku model vystihuje. EE Cep: precedující disk mění hloubku od
zákrytu k zákrytu, takže hloubky sledují precesi. KH 15D: okraj prstence postupně zakrývá dráhu
dvojhvězdy, takže podíl každého oběhu, během něhož je hvězda skryta, se v průběhu let mění. Pro
ilustrativní disk EE Cep (0,8 au) dává moment síly hvězdy B periodu precese 183 let oproti pozorovaným
62–67 letům; souhlasil by disk o velikosti asi 1,6 au – precese tedy omezuje velikost disku.

## Co se nemodeluje {#observables-not-modelled}

Hvězdné skvrny, kombinační frekvence pulzací, světlo a gravitace materiálu disku nad rámec modelu
desky a astrometrie: SWINGERsim nepočítá pohyb fotocentra ani astrometrický signál pro Gaiu. Relativní
i absolutní polohy hvězd na obloze jsou ve výstupu běhu k dispozici a mohly by takovému modelu
posloužit, žádná astrometrická pozorovatelná veličina se však nepočítá ani neprokládá.

## Implementace {#implementation-observables}

**Externí knihovny a data.** Tabulky okrajového a gravitačního ztemnění Clareta & Bloemen (2011)
a koeficienty TLUSTY podle Reeva & Howartha (2016); křivky propustnosti pásem podle Bessella (1990),
TESS (Ricker et al. 2015) a Gaia DR3 (Riello et al. 2021); PHOEBE 2 (Prša et al. 2016; Conroy et al.
2020) ve volitelném hybridu. batman, ellc a PHOEBE slouží jen jako reference v testech. **Vlastní
kód.** Integrace černého tělesa přes fotometrická pásma a interpolace v tabulkách; přesný zákryt
libovolného počtu koulí s okrajovým ztemněním; řešení doby letu světla pro každou dvojhvězdu;
analytické členy efektů blízkosti; Rocheovy plochy a plochy rotujících hvězd, jejich sítě, gravitační
ztemnění a zákryty sledováním paprsků včetně vlastního zákrytu kontaktních dvojhvězd; napojení na
PHOEBE (úseky, převod elementů, vzájemné zákryty); radiální rychlosti v okamžiku vyslání a rychlosti
středu světla na sítích; syntetické okamžiky zákrytů a O−C; fenomenologická minima s Mikuláškovým
profilem; periody pulzací a fázová modulace; a geometrie, propustnost, světlo a tuhá precese disků.
