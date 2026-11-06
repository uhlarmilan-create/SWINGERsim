# Diskuse

## Publikované elementy nejsou bez ukotvení přenositelné

Nejviditelnějším výsledkem je, že publikované elementy, zadané do SWINGERsimu doslovně, dávají zákryty, které
jsou po několika letech mimo o hodiny až dny. Dráhy nejsou špatně: stejná oskulační čísla znamenají mírně
odlišné pohyby ve dvou programech, které pro hierarchii používají odlišné souřadnice a pro slapové a rotační
síly hvězd odlišné modely. Relativní rozdíl několika 10⁻⁴ ve střední periodě je ve světelné křivce jednoho
sektoru TESS neviditelný, ale za tisíce cyklů se nasčítá do hodin. Dvě čísla na pár (perioda a fáze) efekt
úplně odstraní. Kdo přenáší fotodynamické řešení do N-částicového programu, má s tím počítat a řešení nejprve
ukotvit k minimům; SWINGERsim to dělá automaticky, když fit začíná z publikovaných elementů.

## Dynamika je správně; klasický apsidální pohyb testuje stavbu hvězd

Tam, kde apsidální pohyb pohánějí ostatní hvězdy, SWINGERsim reprodukuje publikovanou integraci: BU CMi
s hmotnými body na 1 % a celé ETV TIC 120362137 stejně dobře jako vlastní model autorů. Tam, kde přispívají
samotné hvězdy, se modely liší v tom, co předpokládají. U BU CMi SWINGERsim přidává klasické stáčení, které
předpovídají modely hvězd (2–4 % celku), a fit jej kompenzuje vnější excentricitou (2,6σ). Publikovaná
integrace se chová jako model hmotných bodů SWINGERsimu, ačkoli vlastní analytický odhad slapového stáčení
v článku (0,65 a 0,39°/rok) je stejně velký jako ve SWINGERsimu (0,56 a 0,33°/rok) a publikovaný program
slapy a rotaci zahrnout umí; zda byl jeho slapový člen menší, nebo vypnutý, článek neuvádí, a necháváme
to otevřené. S volnými apsidálními konstantami chtějí minima ještě slabší klasické stáčení
($k_2 \le 0{,}0023$ pro pár A při Δχ² = 4 proti 0,0042 z modelů), ale vnější excentricitu ponechávají
na 0,266, takže konstanty rozdíl nevysvětlují. U V994 Her článek rychlosti stáčení fitoval volně, kdežto SWINGERsim je předpovídá; B a A by
potřebovaly apsidální konstanty změněné asi o 20 % a 10 %, což je pro modely hvězd přijatelné, C o faktor
1,76, který žádný model nedává. Místo toho by to dokázal průvodce C o 0,13–0,31 M☉ na dráze
s periodou 80–182 d; současná minima C jej neukazují ani nevylučují. U B zjistili Zasche et al. (2023) totéž sami: jejich změřená k₂ je
0,83násobkem teoretické hodnoty, což dávají do souvislosti s možnou vyšší metalicitou; faktor SWINGERsimu
je 0,81. V obou soustavách tedy rozdíly nejsou chybami žádného z programů, ale měřením klasického
apsidálního pohybu: časování příštího desetiletí ukáže, jak je velký.

## Co fit dokáže a co ne

- **Dodatečné chyby pro každý zdroj** jsou nejúčinnějším samostatným prvkem fitu. Ve všech třech soustavách
  přivedou redukované χ² na 1, změří, nakolik lze minimům každého pozorovatele věřit (od 1 min u TESS po
  10–18 min u přehlídek a raných pozemních minim), a zabrání několika starým minimům strhávat řešení.
- **Samotné okamžiky minim** určují zákrytové páry dobře (V994 Her: v rámci 2σ od publikovaných hodnot), ale
  ne dráhy viděné jen nepřímo: vnější dráhu TIC 120362137, viditelnou přes jeden zákrytový pár, a dráhu A–B
  soustavy V994 Her, viditelnou přes vlny doby šíření světla. Tam fit zabloudí k řešením, která časování
  vyhovují stejně dobře, ale odporují jiným datům.
- **Apriorní rozdělení** dělají, co mají, u dráhy, kterou minima neurčují (TIC 120362137), a odhalují dráhy,
  které minima určují jinak než články (BU CMi, V994 Her).
- **Radiální rychlosti** dávají hmotnosti, ale testují rychlosti stejně jako model: u TIC 120362137 souhlasí
  na své chyby, u BU CMi mají rozptyl daleko nad uváděnými chybami, u V994 Her samy o sobě souhlasí
  s publikovanými hmotnostmi, ale ve společném fitu je zčásti pohltí dráha A–B.
- **Světelné křivky** určují sklony a poloměry tam, kde jsou zákryty čisté, poměr poloměrů částečně
  zakrývajících se hvězd bez spektroskopického poměru světla neurčí a vyžadují nejprve odstranit jakýkoli
  další signál (oscilaci V994 Her). Je třeba počítat s jejich korelovaným šumem; bez něj byly chyby několikrát
  příliš malé.
- **MCMC s 1000 kroky** dává chyby, které souhlasí s lineárními asi na 10 % (medián), a chyby běhu o 400
  krocích změnilo o několik procent. Řetězce pokrývají celkem asi 19 autokorelačních časů a asi deset po
  zahození poloviny na zahoření, což stačí pro uváděné rozsahy, ale není to formální důkaz konvergence.

## Role pozorovatelů

Minima TESS pokrývají nejvýše pět let. Apsidální pohyb BU CMi (25–29 let) a V994 Her (století a více), 2,9letá
dráha A–B soustavy V994 Her a vnější dráha TIC 120362137 jsou viditelné jen díky pozemním minimům, která
prodlužují základnu na 15–20 let. U V994 Her souhlasí 64 minim R. Uhlaře a minima FRAM s modelem na 2–5 minut,
blízko TESS. Nejstarší minima (2002–2004 u V994 Her, přehlídky z let 2008–2016 u BU CMi) jsou nejméně přesná
a pro dlouhodobé členy nejrozhodnější. Světelné křivky pozorovatelů jsme neměli k dispozici, jen z nich
odvozené okamžiky minim, takže srovnání světelných křivek používá jen TESS; pozemní světelné křivky, např. od
pozorovatelů SQUADRA zpracované pipeline VYVAR, by umožnily totéž srovnání pro zákryty mezi sektory TESS.

## Meze tohoto testu

- Teploty byly převzaty z článků; u V994 Her jsou uzly a dráha C předpoklady.
- Světelné křivky SWINGERsimu používají koule nebo Rocheovy plochy s okrajovým ztemněním; páry v těchto
  soustavách jsou oddělené (Rocheovy laloky vyplněné nejvýše na 57 %) a Rocheovy plochy změnily výsledky
  o méně, než jsou chyby.
- Oscilace V994 Her byla odstraněna jako součet sinusovek, ne modelována.
