---
title: "Testování SWINGERsimu"
subtitle: "Část III: Zákrytové dvojhvězdy s proměnným sklonem – V907 Sco, QX Cas a AH Cep"
---

**Abstrakt.** Má-li zákrytová dvojhvězda třetí hvězdu na skloněné dráze, dráha dvojhvězdy se pomalu stáčí
(precesuje) jako roztočená káča. Její sklon vůči našemu zornému paprsku se mění, takže zákryty se během
desetiletí prohlubují, změlčují nebo mizí a znovu se objevují. Testujeme, zda SWINGERsim, N-částicový
simulátor vázaných soustav tří až deseti hvězd, dokáže takovou soustavu převzít z publikovaných dat,
najít model, který reprodukuje století pozorování, a říci, co přijde dál. Tři dvojhvězdy představují tři
různé situace.

**V907 Sco** má úplné publikované řešení (Zasche et al. 2023). Z 58 okamžiků minim z TESS, 172 radiálních
rychlostí a 14 sklonů měřených od roku 1900 nachází SWINGERsim třetí hvězdu o 0,98 (+0,05/−0,04) M☉ na
dráze s periodou 142,09 ± 0,04 d, skloněné o 25,8 ± 0,8° k dráze dvojhvězdy, s uzlovou periodou 65,5
(+0,9/−1,0) roku. Vše souhlasí s publikovanými hodnotami zhruba v rámci jedné směrodatné odchylky
(publikované rozsahy m₃, vzájemného sklonu a uzlové periody jsou 90% meze). Model vypíná a zapíná zákryty
ve stejných letech, jaké ukazují desky, stará fotometrie a TESS. Předpovídá hloubky zákrytů pozorovaných
TESS v letech 2021 a 2023 s přesností 5 % a klade konec současných zákrytů do let 2052–2057 (publikováno:
kolem roku 2049). Radiální rychlosti Lacyho et al. (1999) odpovídají téže 142denní dráze; fit nastartovaný
na jejich 99denní periodě z ní odejde a skončí ve zřetelně horším řešení.

**QX Cas** přestala zakrývat v 90. letech a její třetí hvězda nebyla nikdy pozorována (Guinan et al.
2012a,b). SWINGERsim ukazuje, které třetí hvězdy dokážou dráhu stáčet dostatečně rychle (nejpravděpodobněji
0,5–3 roku a 0,5–4 M☉), a že sklon musel klesat nejméně dvakrát rychleji, než udává publikovaná lineární
rychlost. Dvojhvězdu také nachází ve světelných křivkách TESS z let 2019–2024 podle změn její jasnosti na
oběžné periodě, s amplitudou, kterou model předpovídá. Výpočet **hloubky** zákrytů, nejen toho, zda je
zákryt geometricky možný, odhaluje slabé místo: v našem nejlepším modelu jsou zákryty z 30. let příliš
mělké pro fotografické desky, které je zaznamenaly. Model s hlubokými zákryty ve 30. letech existuje, ale
zákryty si podrží až do konce 90. let, kdy žádné pozorovány nebyly. Staré hloubky zákrytů navíc vyžadují
13–47 % třetího světla, víc, než dává kterákoli třetí hvězda modelu. Obojí ukazuje na převzaté parametry
hvězd, ne na dynamiku.

**AH Cep** je hmotná dvojhvězda, u níž byla kdysi ohlášena změna hloubky zákrytů. SWINGERsim reprodukuje
125 let okamžiků zákrytů třetím tělesem na 70leté dráze o hmotnosti nejméně 6,6 M☉. Tak vzdálený průvodce
mění sklon jen o 0,02° za 400 let, takže ohlášenou změnu způsobit nemůže. Čtvrté těleso, které by pár
znatelně naklánělo, by se projevilo v okamžicích zákrytů nebo v hloubkách zákrytů v TESS, a kdysi navržené
čtvrté těleso s periodou 9 let by bylo dynamicky nestabilní. Světelné křivky TESS nejlépe popisuje sklon
asi 72,50 ± 0,21° (publikováno 69,8°) s asi 5 % světla navíc; to je v rozporu s třetím světlem < 2 %
připouštěným spektry a závisí to na modelu světelné křivky.

**Gaia DR3** neukazuje u žádné ze tří soustav RUWE nad 1,4. U V907 Sco, kde náš model posouvá světelný střed
asi o 0,3 mas, je to rozpor, který stojí za to sledovat v Gaia DR4.

# Cíl

Předchozí části testu SQUADRA ověřovaly SWINGERsim na kompaktních vícenásobných hvězdách s úplnými
publikovanými řešeními. Zde je otázka těžší: **když máme to, co bylo naměřeno za století – několik
sklonů, seznam okamžiků zákrytů velmi různé kvality, nějaké radiální rychlosti –, dokáže SWINGERsim najít
soustavu, která to vytvořila, ověřit publikované závěry a říci něco nového?**

U každé soustavy se ptáme na čtyři věci:

1. K jakým závěrům dospívají publikované práce a z jakých dat?
2. Co dostane SWINGERsim ze stejných dat?
3. Práce SWINGERsim **potvrzuje**, něco k nim **přidává**, nebo s nimi **nesouhlasí**, a proč?
4. Co SWINGERsim předpovídá, co lze ověřit na již existujících datech (zde hlavně TESS a Gaia)?
