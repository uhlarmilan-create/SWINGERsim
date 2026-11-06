# Data a metody

## Tři soustavy

**Tabulka 1.** Testované soustavy a použitá data. *Minima*: publikované okamžiky minim (počet, roky,
zdroje). *TESS*: sektory srovnání světelných křivek (SPOC: zpracování TESS Science Processing Operations
Center; QLP: Quick-Look Pipeline).

| | BU CMi | TIC 120362137 | V994 Her |
|---|---|---|---|
| hierarchie | 2+2: (Aa+Ab)+(Ba+Bb) | 3+1: ((Aa+Ab)+B)+C | (2+2)+2: ((A)+(B))+(C) |
| zákrytové páry | A 2,94 d, B 3,26 d | Aa–Ab 3,28 d | A 2,08 d, B 1,42 d, C 1,96 d |
| vnější dráhy | A–B 121,8 d | A–B 51,3 d, AB–C 1046 d | A–B 1062 d; C kolem AB (neměřeno) |
| zdroj řešení | Pribulla et al. (2023), tab. 7 | Borkovits et al. (2026), tab. 1–2 | Zasche et al. (2023), tab. 3 |
| epocha elementů (BJD) | 2458491,5 | 2458683,3 | 2459010,0 |
| minima | 183, 2005–2023; HAT, KELT, MASCARA, TESS, pozemní dopozorování, Ogloza et al. (2008), Volkov et al. (2021) | 147, 2019–2024; doplňková tabulka S3 / Supplementary Data 6 článku | 343, 2002–2022; TESS, pozorovatelé RU, SD & UM, FRAM, MM, GP |
| radiální rychlosti | 305 všech čtyř hvězd (Pribulla et al. 2023) | 380 všech čtyř hvězd (Borkovits et al. 2026) | 128 hvězd párů A a B (Lee et al. 2008) |
| sektory TESS | 7 (SPOC 30 min), 34, 61 (SPOC 2 min) | 14, 26, 40–41, 53–54, 74 (SPOC) | 26 (SPOC 2 min), 40, 53 (QLP 10 min) |

Publikovaná řešení všech tří soustav byla získána programem Lightcurvefactory (Borkovits et al. 2019,
2020, 2021), který fituje světelné křivky, radiální rychlosti, změny okamžiků zákrytů a spektrální
rozdělení energie společně, s poloměry a teplotami hvězd převzatými z hvězdných izochron (PARSEC). Dráhy se
počítají dvěma různými způsoby. U BU CMi a TIC 120362137 program integruje pohyb hvězd numericky
(*fotodynamický* model); jeho integrátor umí zahrnout slapové a rotační efekty (Borkovits et al. 2019,
příloha A), ale ani jeden článek neuvádí, zda byly tyto členy použity. U V994 Her jsou poruchy mezi páry
malé a Zasche et al. (2023, odd. 3) místo N-částicové integrace použili analytické keplerovské dráhy
s konstantní, volně fitovanou rychlostí apsidálního stáčení každého páru. Dráhové elementy v článcích jsou
*okamžité oskulační elementy* v epoše z tabulky 1, tj. keplerovská dráha, po které by se hvězdy od toho
okamžiku pohybovaly, kdyby ostatní hvězdy zmizely. SWINGERsim používá stejný druh elementů jako počáteční
podmínky, s argumentem periastra ω sekundární složky (relativní dráhy), podle konvence všech tří článků.
U V994 Her články neuvádějí uzly; předpokládali jsme Ω = 0 pro všechny dráhy (autoři považují soustavu za
téměř koplanární) a dráhu C kolem A+B s P = 3000 let a e = 0,3 (promítnutá vzdálenost je ~300 au); ansámbl
tyto předpoklady obměňuje.

## Pojmy použité ve zprávě

- **Okamžik minima (okamžik zákrytu).** Okamžik největšího poklesu jasnosti během zákrytu, změřený
  pozorovatelem ze světelné křivky: to jsou **data**. Články je uvádějí s chybami; *primární* minimum je
  zákryt první hvězdy páru.
- **Modelový okamžik zákrytu.** Pro každé naměřené minimum SWINGERsim spočítá okamžik nejmenší promítnuté
  vzdálenosti obou hvězd páru v N-částicové integraci plus dobu, kterou světlo potřebuje k pozorovateli
  napříč soustavou (Irwin 1952). Je to modelová podoba téže události, ne předpověď.
- **O − C** („pozorováno minus počítáno“). Naměřený okamžik minus okamžik přísně periodické (lineární)
  efemeridy T = T0 + P E, kde E je číslo cyklu. Přísně periodická dvojhvězda dává vodorovnou čáru. Ve
  vícenásobné hvězdě se čára ohýbá: *apsidální pohyb* posouvá primární a sekundární minima opačnými směry
  (dvě zrcadlové vlny); *efekt světelného času* (LTTE) vnější dráhy posouvá všechna minima páru
  společně; a v kompaktních soustavách třetí hvězda přímo mění periodu (*dynamické* změny okamžiků zákrytů,
  ETV). V této zprávě je lineární efemerida každého páru střední efemeridou modelových primárních minim,
  takže naměřené body i modelové křivky mají společný základ.
- **Předpověď.** Modelový okamžik zákrytu nebo světelná křivka k datu bez dat.
- **Oskulační a střední perioda.** Oskulační perioda je perioda keplerovské dráhy v jednom okamžiku;
  střední perioda je průměrný interval mezi zákryty. V kompaktní vícenásobné hvězdě se liší o 10⁻⁴–10⁻³
  periody, protože dráha páru je neustále rušena.
- **Ukotvení.** Úprava jen periody a fáze (střední anomálie v epoše) každého zákrytového páru tak, aby
  jeho střední perioda zákrytů a fáze souhlasily s minimy.
- **Apsidální pohyb.** Stáčení přímky apsid (dlouhé osy) excentrické dráhy; jeho perioda U je doba jednoho
  úplného otočení. Pohánějí jej ostatní hvězdy (*dynamický*), slapové a rotační výdutě samotných hvězd
  (*klasický*) a obecná relativita.
- **Poměr stability.** Pro pár uvnitř vnější dráhy poměr vnější pericentrové vzdálenosti ke kritické
  vzdálenosti podle Mardlinga & Aarsetha (2001, „MA01“) nebo Vynatheyi et al. (2022, „V22“); nad 1 je trojice
  stabilní. Pro čtyřhvězdy SWINGERsim používá také klasifikátor strojového učení Vynatheyi et al. (2023) a
  *test duchem*: dva N-částicové běhy lišící se o 10⁻⁶ au v jedné velké poloose musí zůstat blízko sebe.
- **Světlo navíc.** Světlo ve fotometrické cloně, které se zákrytů neúčastní: další hvězdy soustavy,
  sousední hvězdy nebo světlo ubrané či přidané korekcí zpracování na přeplněné pole. Mění všechny hloubky
  zákrytů stejným faktorem.
- **ppt**: promile celkového toku. **rms**: střední kvadratická hodnota reziduí. **Medián |O − C|**:
  typické reziduum, méně citlivé na několik špatných bodů.

## Tři konfigurace modelu

Každá soustava byla počítána se třemi konfiguracemi:

- **publikovaná**: přesně elementy, hmotnosti, poloměry a teploty z článku;
- **ukotvená**: publikovaná konfigurace s periodou a fází každého zákrytového páru přizpůsobenými minimům;
  všechna ostatní čísla jako v článku;
- **nafitovaná**: výsledek fitu SWINGERsimu (níže).

Proč je ukotvení nutné, je jedním z výsledků testu: oskulační periody nejsou přenositelné. Dva programy,
které popisují týž pohyb v mírně odlišných souřadnicích nebo s mírně odlišnými silami, přiřadí stejným
oskulačním číslům mírně odlišné střední periody, a za několik tisíc cyklů se rozdíl 10⁻⁴ v periodě změní
v hodiny. Ukotvení tento rozdíl odstraní dvěma čísly na pár a dynamiku nechá nedotčenou.

## Model

SWINGERsim integruje pohyb všech hvězd s plnou gravitací (REBOUND s integrátorem IAS15, Rein & Liu 2012, Rein
& Spiegel 2015), včetně obecné relativity a slapové a rotační deformace hvězd (REBOUNDx, Tamayo et al.
2020), s apsidálními konstantami z modelů hvězd (Claret 2023) ve stáří podle stop MIST (Choi et al. 2016;
Dotter 2016). Světelná křivka sčítá světlo všech hvězd: koule (nebo Rocheovy plochy) s okrajovým ztemněním
podle čtyřkoeficientového zákona Clareta (2000) s koeficienty Clareta & Bloemen (2011), se zákryty každé
hvězdy každou jinou hvězdou, která je před ní, s elipsoidální proměnností, odrazem a dopplerovským
zjasněním těsných párů a s dobou šíření světla napříč soustavou.

## Srovnání světelných křivek

Světelné křivky družice TESS (Ricker et al. 2015) byly staženy programem lightkurve (Lightkurve
Collaboration 2018) ze zpracování SPOC (Jenkins et al. 2016) a QLP (Huang et al. 2020). Model byl spočítán
v pozorovaných časech a zprůměrován přes expozici (2, 10 nebo 30 min). Dvě věci ve fotometrii TESS nejsou
vlastnostmi hvězd a je třeba s nimi počítat: světlo navíc v cloně mění všechny hloubky stejným faktorem
a pomalé instrumentální trendy ohýbají základní úroveň. Pozorování jsme proto srovnávali
s `a × model + baseline` (základní úroveň), kde `a` je jeden škálovací faktor na sektor a základní úroveň je
kvadratický polynom v čase pro každý souvislý úsek. To jsou jediná upravovaná čísla; zákryty přicházejí
z modelu nezměněné. Měřítkem shody je rms rozdílu pozorování a modelu v ppt, po odstranění bodů vzdálených
více než pět robustních směrodatných odchylek.

## Fit

Fit mění elementy N-částicového modelu, dokud nesouhlasí s daty. Pro každé naměřené minimum se soustava
zintegruje těsně před očekávaný zákryt a modelový zákryt se najde z pohybu páru. Fit minimalizuje χ² =
Σ ((data − model)/σ)² přes všechna data najednou:

- **okamžiky minim**, každý zdroj (každá přehlídka, každý pozorovatel, každý sektor TESS) s vlastní
  dodatečnou chybou („jitter“), přičtenou kvadraticky k uváděným chybám a určenou fitem tak, aby každý zdroj
  měl redukované χ² rovné 1; minima vzdálená od modelu více než 5 robustních směrodatných odchylek jsou
  vynechána (a vrácena, pokud je pozdější iterace přiblíží);
- **radiální rychlosti**, s volnými hmotnostmi hvězd, systémovou rychlostí a nulovým bodem každého
  přístroje;
- **světelné křivky**: okna ±0,3 d kolem zákrytů jednoho nebo dvou sektorů TESS s volnými sklony a poloměry
  a lineární kalibrací každého okna (škála pro světlo navíc, nulový bod, trend). Sousední body světelné
  křivky nejsou nezávislé (červený šum); jejich chyby se násobí faktorem β podle Ponta et al. (2006),
  změřeným z reziduí;
- **apriorní rozdělení**: gaussovská apriorní rozdělení z článků na elementy, které data omezují jen slabě
  (excentricity a ω párů, vnější dráhy), nikdy na periody párů (oskulační periody nejsou přenositelné).

Minimalizace používá metodu nejmenších čtverců s oblastí důvěry (SciPy, Virtanen et al. 2020).
Aposteriorní rozdělení bylo poté vzorkováno programem emcee (Foreman-Mackey et al. 2013; 72 chodců × 1000
kroků, první polovina zahozena jako zahoření). Autokorelační čas řetězců je τ ≈ 52–53 kroků, takže celý
řetězec pokrývá asi 19 τ a ponechaná druhá polovina asi deset τ na chodce. To splňuje praktické pravidlo
Foreman-Mackeyho et al. (2013), kteří doporučují běh o délce asi deseti τ, a stačí pro zde uváděné rozsahy
16–84 % (asi 700 nezávislých vzorků), ale nedosahuje 50 τ, které pro spolehlivý odhad samotného τ doporučuje
online dokumentace emcee (návod „Autocorrelation analysis & convergence“). Běh o 400 krocích dal chyby,
které se od chyb z 1000 kroků liší v mediánu o 3–5 % (nejvýše dvakrát, u vektoru excentricity téměř
kruhového páru TIC 120362137), a mediány, které se liší o méně než 0,4σ. Uváděné chyby jsou proto
spolehlivé na několik desítek procent.

Aby bylo vidět, co který druh dat přináší, fit byl spuštěn i po krocích: jen okamžiky minim, pak
s dodatečnou chybou pro každý zdroj, s apriorními rozděleními, automaticky z publikovaných elementů,
s radiálními rychlostmi a se světelnými křivkami. Kroky jsou shrnuty u každé soustavy.

## Ansámbly

Ansámbl je množina soustav vylosovaných náhodně v rámci publikovaných nejistot (hmotnosti, excentricity,
sklony, uzly, vnější dráhy; u V994 Her i neznámá dráha C). Každá soustava se zintegruje dopředu i dozadu
a SWINGERsim zaznamená její osud (neporušená, nestabilní, splynulá), její periody precese a zda její páry
zakrývají. Soustavy se pak váží podle toho, jak dobře reprodukují publikované apsidální periody a podmínku,
že páry dnes zakrývají. Efektivní počet členů n_eff udává, kolik soustav nese váhu.
