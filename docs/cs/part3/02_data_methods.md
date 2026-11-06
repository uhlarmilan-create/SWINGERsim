# Data a metody

## Jak SWINGERsim modeluje soustavu

SWINGERsim popisuje soustavu jako hierarchii drah: zde těsný pár hvězd (zákrytová dvojhvězda) a třetí hvězda
obíhající kolem páru. Z hmotností, drah a poloměrů a teplot hvězd integruje pohyb všech hvězd s plnou
gravitací (N-částicově; REBOUND s integrátorem IAS15, Rein & Liu 2012, Rein & Spiegel 2015). Pohyb
zahrnuje obecnou relativitu a slapovou a rotační deformaci hvězd (REBOUNDx, Tamayo et al. 2020), které
pomalu stáčejí excentrické dráhy (apsidální pohyb). Ze zintegrovaného pohybu SWINGERsim počítá okamžiky
zákrytů, radiální rychlosti hvězd, sklon dráhy páru k libovolnému datu a syntetické světelné křivky v daném
fotometrickém pásmu (zde TESS, Johnson B a V). Hloubky zákrytů plynou z geometrie a povrchů hvězd (koule,
nebo Rocheovy plochy u těsných hvězd).

## Pojmy použité ve zprávě

Několik slov používáme v přesném smyslu. Zde je jednou definujeme.

- **Okamžik minima (okamžik zákrytu).** Okamžik nejmenší jasnosti během zákrytu, změřený pozorovatelem
  ze světelné křivky. To jsou **data**. Pocházejí z literatury, z O − C brány České astronomické
  společnosti (var.astro.cz), z O − C souborů R. H. Nelsona (AAVSO) a u TESS z našeho vlastního měření:
  pro každý sektor TESS se každým zákrytem proloží střední profil zákrytu s volným časem.
- **Modelový okamžik zákrytu.** Pro každé pozorované minimum model spočítá okamžik, kdy obě hvězdy páru
  stojí při pohledu pozorovatele v zákrytu za sebou (konjunkce), nalezený v N-částicové integraci, plus
  dobu, kterou světlo potřebuje k pozorovateli z polohy páru na jeho dráze kolem třetí hvězdy. Počítá se
  v epoše každého pozorovaného minima a srovnává se s ním; **není to předpověď**, ale modelová podoba
  téže události. Geometrická konjunkce a fotometrické minimum se u těchto soustav liší o sekundy (asi 1 s
  u AH Cep, ze syntetických světelných křivek), hluboko pod chybami dat.
- **O − C** („pozorováno minus počítáno“). Rozdíl mezi okamžikem zákrytu a přísně periodickými hodinami
  (lineární efemeridou). Ukazuje, o kolik minut zákryty přicházejí dříve nebo později. O − C modelu
  používá tytéž hodiny, takže data i model lze nakreslit do jednoho diagramu.
- **Předpověď.** Modelový okamžik zákrytu, sklon nebo světelná křivka k datu **bez dat**. V této zprávě
  předpovídáme jen mimo pozorované intervaly a takto to i označujeme.
- **Mez zákrytů.** Pár s poloměry R₁, R₂ zakrývá jen tehdy, je-li promítnutá vzdálenost obou hvězd
  v konjunkci menší než R₁ + R₂, tj. když cos i < (R₁ + R₂)/r, kde r je jejich vzdálenost v konjunkci.
  U V907 Sco to je i > 72,8°, u QX Cas i > 70,7° (pomalu se mění, protože se excentrická dráha páru
  stáčí), u AH Cep i > 50,4°. Nad mezí hvězdy zakrývat **mohou**; **neznamená** to, že zákryt byl dost
  hluboký, aby byl pozorován. Blízko meze jsou zákryty tečné a velmi mělké. V celé zprávě proto
  rozlišujeme **geometrický zákryt** (sklon je nad mezí) a **zjistitelný zákryt** (jeho hloubka je větší,
  než jakou tehdejší fotometrie dokázala zachytit: zhruba 0,1 mag u fotografických desek, několik setin
  magnitudy u fotoelektrické fotometrie, několik desetin procenta u TESS).
- **Třetí světlo.** Jakékoli světlo, které se zákrytu neúčastní: třetí hvězda soustavy, hvězda ve stejné
  fotometrické cloně nebo rozptýlené světlo. Tvoří-li třetí světlo podíl l₃ celkového světla, zákryt,
  který by ubral podíl d světla dvojhvězdy, ubere jen (1 − l₃) d celkového světla a je v magnitudách
  mělčí. Hloubka zákrytu proto vždy závisí na třech věcech zároveň: na sklonu, na velikostech
  a povrchových jasnostech obou hvězd a na třetím světle. Kdykoli tato zpráva uvádí hloubku zákrytu,
  říká, jaké třetí světlo bylo předpokládáno.
- **Uzlová perioda.** Doba, za kterou se dráha páru jednou otočí kolem celkového momentu hybnosti
  soustavy. Během jednoho otočení se sklon páru kývá mezi dvěma krajními hodnotami.

## Sklony nad 90°

Světelná křivka měří jen sin i: dráha se sklonem i a dráha se sklonem 180° − i dávají stejné zákryty.
Publikované sklony jsou proto vždy ≤ 90° a hlavní panel každého obrázku se sklonem ukazuje model v této
**pozorovatelné** podobě (0–90°). Když dráha během precese prochází polohou, v níž ji vidíme přesně z boku
(90°), pozorovatelný sklon stoupne na 90° a klesne zpět, a poté, co dráha mine svůj bod obratu a vrací se
přes polohu z boku, stoupne znovu. Během jednoho takového průchodu má proto pozorovatelná křivka dva vrcholy
v 90° s prohlubní mezi nimi, tvar „M“.

N-částicový model ví i to, na které straně oblohy dráha je. Druhý panel ukazuje tutéž křivku **rozloženou**
do rozsahu 0–180°: tam je to jediná hladká vlna a prohlubeň „M“ se ukáže jako skutečné maximum sklonu nad
90°. Oba panely obsahují stejnou informaci. Pozorovatelný panel je ten, se kterým lze srovnávat světelnou
křivku; rozložený panel ukazuje precesi.

## Jak SWINGERsim fituje data

Fit mění hmotnosti a dráhy N-částicového modelu, dokud model nesouhlasí se všemi daty najednou:

- okamžiky minim. Každý zdroj (vizuální, fotografický, fotoelektrický, CCD, každý sektor TESS) má
  vlastní dodatečnou chybu, kterou určí fit, a odlehlé body za 4–5σ jsou vyřazeny;
- radiální rychlosti (každá sada dat s vlastním nulovým bodem a vlastní dodatečnou chybou);
- sklony z publikovaných řešení světelných křivek v jejich epochách. Používá se i jednostranná informace:
  „sklon byl pod 75°“ nebo „v tomto roce nebyly zákryty“ (model pak musí zůstat pod svou vlastní mezí
  zákrytů);
- světelné křivky, kde je to uvedeno.

Bez astrometrie je měřitelný jen rozdíl uzlů obou drah na obloze, takže uzel páru je fixován.

**Chyby.** Není-li uvedeno jinak, pocházejí chyby ze zakřivení χ² kolem nejlepšího fitu (lineární chyby).
U V907 Sco bylo aposteriorní rozdělení navíc vzorkováno programem emcee (Foreman-Mackey et al. 2013;
48 chodců × 1000 kroků, prvních 500 kroků zahozeno jako zahoření). Autokorelační čas řetězce je τ ≈ 54
kroků, takže celý řetězec pokrývá asi 19 τ a 500 ponechaných kroků asi 9 τ, tj. asi devět nezávislých
vzorků na chodce a celkem asi 430. Foreman-Mackey et al. (2013) doporučují běh o délce asi deseti τ, čehož
ponechaná část téměř dosahuje. To stačí pro zde uváděné rozsahy 16–84 %, ale ne pro 50 τ, které pro
spolehlivý odhad samotného τ doporučuje online dokumentace emcee (návod „Autocorrelation analysis
& convergence“). Řetězec o 400 krocích dal chyby, které se od chyb z 1000 kroků liší v mediánu o 10 %
(nejvýše 1,5krát, u m₃ a fáze vnější dráhy), a mediány, které se liší o méně než 0,8σ. Uváděné chyby jsou
proto spolehlivé na několik desítek procent. 100 vzorků tohoto aposteriorního rozdělení bylo zintegrováno
přes roky 1890–2070, aby daly rozsahy předpovězeného sklonu a přechodů zákrytů.

Tam, kde jsou data na jednoznačný fit příliš chudá (QX Cas), a pro otázku čtvrtého tělesa (AH Cep)
sestaví SWINGERsim **ansámbl**: tisíce možných soustav vylosovaných ze širokých rozsahů se zintegrují,
soustavy nestabilní podle kritéria Mardlinga & Aarsetha (2001) se vyřadí a zbytek se váží podle toho, jak
dobře reprodukuje měření. Jeden fit je jedna možná soustava; ansámbl ukazuje rozsah soustav, které data
připouštějí.

## Soustavy a jejich data

**Tabulka 1.** Tři soustavy a použitá data.

| | V907 Sco | QX Cas | AH Cep |
|---|---|---|---|
| zákrytový pár | 3,78 d, e ≈ 0,01, 2,74 + 2,56 M☉, pozdní B (Zasche et al. 2023) | 6,00 d, e ≈ 0,21, ≈ 5,5 + 6,5 M☉, B1.5 V + B3 IV (Guinan et al. 2012a) | 1,775 d, kruhová, 16,1 + 13,7 M☉, B0.5 V (Pavlovski et al. 2018) |
| třetí těleso | 142 d, e = 0,47, 1,06 M☉ (Zasche et al. 2023) | nezjištěno | dráha z efektu světelného času, ~68 let (Kim et al. 2005) |
| sklony | 14 epoch 1900–2022 (Zasche et al. 2023, tab. 2) | ~1960, ~1980, bez zákrytů od 2003 (Guinan et al. 2012a) | ~1983 (Pavlovski et al. 2018) |
| zákryty pozorovány | 1899–1918, 1963–1986 (Lacy et al. 1999), od 2019 (TESS) | fotografická minima 1933–1974; fotoelektricky ~1960 a ~1980; žádné v polovině 90. let, 2003–2008 (Guinan et al. 2012a,b) a TESS 2019–2024 | vždy |
| okamžiky minim | 58 z TESS 2021–2025 | 96 z let 1933–1974 (74 fotografických) | 325 z let 1901–2026 (216 TESS) |
| radiální rychlosti | 172: CHIRON 2022 (Zasche et al. 2023), 1983–1994 (Lacy et al. 1999) | – | – |
| světelné křivky | TESS 2021, 2023 | TESS 2019–2024 (bez zákrytů) | TESS 2019, 2024 |

Data TESS (Ricker et al. 2015) jsou světelné křivky SPOC, není-li uvedeno jinak.

## Jak číst obrázky soustav

Každá soustava je zobrazena ve třech panelech (obr. 1, 8, 12), publikovaná soustava čárkovaně a oranžově,
výsledek SWINGERsimu plnou čarou a modře:

- **Shora.** Dráhy viděné z bodu nad soustavou, se Zemí směrem ke spodnímu okraji panelu. Pár (popisek
  „pair m1+m2“) a třetí hvězda se pohybují kolem společného těžiště.
- **Zákrytový pár na obloze.** Co by viděl dalekohled s neomezenou rozlišovací schopností: dráhy obou
  hvězd páru promítnuté na oblohu, s hvězdami nakreslenými v měřítku. Zákryt nastane, když jeden disk
  přejde před druhým. Čím je promítnutá dráha plošší, tím blíže se díváme podél roviny dráhy a tím hlubší
  jsou zákryty.
- **Kam míří dráha páru.** Osa dráhy páru (jednotkový vektor) promítnutá na oblohu. Bod ve vzdálenosti
  sin i od středu znamená sklon i; okraj (vzdálenost 1) znamená, že se díváme přesně podél roviny dráhy.
  Zákryty jsou možné, jen když osa leží ve žlutém prstenci mezi okrajem a tečkovanou kružnicí o poloměru
  sin(mez zákrytů). Jak dráha precesuje, osa se pohybuje kolem celkového momentu hybnosti soustavy (×);
  popisky jsou roky a velká tečka označuje rok 2025. Panel je přiblížen na dráhu osy.
