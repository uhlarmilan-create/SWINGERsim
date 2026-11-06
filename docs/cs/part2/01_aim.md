---
title: "Testování SWINGERsimu"
subtitle: "Část II: Tři publikované vícenásobné hvězdy – BU CMi, TIC 120362137 a V994 Her"
---

**Abstrakt.** Testujeme, zda N-částicový simulátor SWINGERsim reprodukuje tři dobře prostudované kompaktní
vícenásobné hvězdy, jejichž řešení publikovali členové a spolupracovníci skupiny SQUADRA: čtyřhvězdu 2+2
BU CMi (Pribulla et al. 2023), čtyřhvězdu 3+1 TIC 120362137 (Borkovits et al. 2026) a trojnásobně
zákrytovou šestinásobnou hvězdu V994 Her (Zasche et al. 2023). Pro každou soustavu SWINGERsim integroval
publikované řešení dopředu i dozadu v čase, otestoval její stabilitu, spočítal syntetické světelné křivky
TESS a okamžiky všech zákrytů, nafitoval vlastní N-částicový model na naměřené okamžiky minim, radiální
rychlosti a světelné křivky TESS a zintegroval ansámbl 400–600 soustav vylosovaných v rámci publikovaných
nejistot.

Všechny tři soustavy jsou stabilní ve všech testech. Publikované elementy, převzaté doslovně, dávají
zákryty, které se během několika let posunou o hodiny až dny, protože oskulační periody nejsou periodami,
se kterými se zákryty opakují; po přizpůsobení periody a fáze každého zákrytového páru minimům („ukotvení“)
souhlasí syntetické světelné křivky TESS s pozorovanými v nejlepších sektorech na 2–4 promile toku (ppt),
zhruba na jedno procento hloubky zákrytu, a okamžiky zákrytů na rozptyl dat. Společný fit všech dat,
s realistickou chybou pro každý zdroj minim a aposteriorním rozdělením vzorkovaným metodou MCMC,
reprodukuje publikované hmotnosti na 1–5 % (BU CMi, TIC 120362137) a sklony na 0,15°. Ukazuje také tři
rozdíly, z nichž každý jsme sledovali tak daleko, jak to publikované informace dovolují. (1) U BU CMi
souhlasí publikované apsidální periody na 1 % s integrací SWINGERsimu pro hmotné body. Když SWINGERsim přidá
stáčení přímky apsid způsobené slapovými a rotačními výdutěmi hvězd, fit reprodukuje pozorovaný apsidální
pohyb tak, že sníží vnější excentricitu o 2,6σ; s hmotnými body rozpor zmizí. Publikovaný kód umí slapové
a rotační efekty zahrnout a vlastní analytický odhad slapového stáčení v článku je dokonce o něco větší než
odhad SWINGERsimu; zda publikovaná integrace použila menší slapový člen, nebo žádný, článek neuvádí, a tato
otázka zůstává otevřená. (2) U TIC 120362137 je publikovaná apsidální perioda 51denní dráhy teoretická
hodnota, ne měření; 32,4 roku ze SWINGERsimu tedy s daty v rozporu není. (3) U V994 Her je apsidální pohyb
převážně klasický: SWINGERsim jej předpovídá ze stavby hvězd, kdežto článek rychlosti volně fitoval. Pro pár B
nachází sám článek apsidální konstantu asi 0,83krát menší než teoretickou, blízko faktoru 0,81, který
potřebuje SWINGERsim; pár C by potřeboval apsidální konstanty 1,8krát větší, než dávají hvězdné modely.

# Cíl testu

SWINGERsim je dopředný model: z hmotností, poloměrů a drah počítá, jak se hierarchická vícenásobná hvězda
pohybuje, jak vypadá její světelná křivka a kdy nastávají její zákryty. SQUADRA pracuje opačným směrem:
vychází z pozorování, světelných křivek a okamžiků minim naměřených družicemi a pozorovateli, a ptá se,
jaká soustava je vytvořila. Než může SWINGERsim s touto úlohou pomoci, musí projít jednoduchým testem:
**reprodukuje SWINGERsim to, co bylo pozorováno, u soustavy, kterou už jiní pečlivě vyřešili?**

Vybrali jsme tři soustavy s úplnými publikovanými řešeními, dlouhými řadami okamžiků minim a světelnými
křivkami TESS, po jedné z každé hierarchie, na které SQUADRA záleží:

- **BU CMi**: dvě zákrytové dvojhvězdy obíhající jedna kolem druhé (2+2), nejkompaktnější známá plochá
  dvojitě zákrytová čtyřhvězda (vnější perioda 121,8 d);
- **TIC 120362137**: zákrytová dvojhvězda se třetí hvězdou na 51denní dráze a čtvrtou na 1046denní dráze
  (3+1), nejkompaktnější známá čtyřhvězda tohoto typu;
- **V994 Her**: tři zákrytové dvojhvězdy, A a B v uspořádání 2+2 a vizuální průvodce C ((2+2)+2), první
  známá trojnásobně zákrytová šestinásobná hvězda.

U každé soustavy se ptáme na pět otázek:

1. **Stabilita.** Je publikovaná soustava dynamicky stabilní podle analytických kritérií a podle přímého
   N-částicového testu?
2. **Dynamika v čase.** Co se děje s drahami během století dopředu a dozadu a souhlasí periody precese
   spočítané SWINGERsimem s publikovanými?
3. **Světelné křivky.** Souhlasí syntetické světelné křivky SWINGERsimu, celé soustavy i jednotlivých párů, se
   světelnými křivkami TESS?
4. **Okamžiky minim.** Nastávají zákryty SWINGERsimu tehdy, kdy je pozorovatelé viděli, po celé době pokryté
   publikovanými minimy?
5. **Fit.** Dokáže vlastní fit SWINGERsimu na všechna tato data najít konfiguraci, která je reprodukuje, a jak
   se srovná s publikovaným řešením? Ansámbl soustav v rámci publikovaných chyb ukazuje, jak robustní jsou
   závěry.

Úspěch znamená: publikovaná soustava je v SWINGERsimu stabilní; syntetické světelné křivky souhlasí s TESS na
úrovni reziduí publikovaných modelů (několik ppt); okamžiky zákrytů souhlasí s naměřenými minimy na
minuty, jakmile jsou periody ukotveny; a fit reprodukuje data, aniž by opustil rozsah fyzikálně rozumných
parametrů. Kde se SWINGERsim s článkem neshoduje, hledáme důvod, ne shodu.
