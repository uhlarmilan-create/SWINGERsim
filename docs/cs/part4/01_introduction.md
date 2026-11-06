---
title: "Testování SWINGERsimu"
subtitle: "Část IV: Zvláštní případy – budoucnost CzeV343, hierarchie CzeV4315 a konvergence aposteriorních rozdělení"
---

# Cíl {#aim}

Části I–III této série testů se ptaly, zda SWINGERsim reprodukuje to, co je již známo: publikované
srovnávací úlohy fyziky (část I), tři dobře pozorované vícenásobné hvězdy (část II) a tři zákrytové
dvojhvězdy, jejichž sklon se mění působením třetího tělesa (část III). Část IV klade otázky, na které
publikovaná odpověď neexistuje. Právě pro takové otázky byl program vytvořen.

1. **Jak a kdy se CzeV343 rozpadne?** CzeV343 je dvojitě zákrytová čtyřhvězda typu 2+2.
   Pejcha et al. (2022) zjistili, že se čtyřhvězda stane dynamicky nestabilní, jestliže se jedna z jejích
   dvojhvězd přenosem hmoty rozšíří nad asi 50 R☉. Neuvedli, kdy k tomu dojde ani co bude následovat
   („dynamical instability caused by orbital expansion … has not been explored so far in 2+2
   quadruples“). Sledujeme vývoj hvězd a přenos hmoty v obou dvojhvězdách až do počátku nestability.
   Od tohoto okamžiku pak integrujeme pohyb všech čtyř těles
   ([oddíl 2](#czev343)).
2. **Je CzeV4315 čtyřhvězda 2+2, nebo pětihvězda 2+(2+1)?** V CzeV4315 se sklon dvojhvězdy B mění
   o 1,4° za rok (Zasche et al. 2025). Okamžiky zákrytů B navíc vykazují 119denní cyklus. Srovnáváme
   dvě hierarchie se stejnými daty, vážíme je bayesovskou evidencí a ptáme se, co předpovídá ta
   preferovaná ([oddíl 3](#czev4315)).
3. **Jsou aposteriorní rozdělení fitů zkonvergovaná?** Fity metodou Monte Carlo s Markovovými řetězci
   (MCMC) v částech II a III běžely 1000 kroků. Ověřujeme, jak daleko mají ke konvergenci a co to
   znamená pro uváděné nejistoty ([oddíl 4](#mcmc)).

Ke dvěma druhům výsledků v této části. **Dopředný model** vychází ze současného stavu soustavy
a integruje ji do budoucnosti; jeho odpovědí je rozdělení přes neznámé parametry, nikoli jediné datum.
**Ansámbl** je množina takových modelů, jejichž neznámé parametry jsou taženy z uvedených apriorních
rozsahů. Každý člen dostane **váhu pozorování**: věrohodnost pozorovaných dat pro daného člena.
Průměrná váha přes ansámbl odhaduje **bayesovskou evidenci** $Z$ hypotézy, na níž je ansámbl
postaven; poměr dvou evidencí, **Bayesův faktor**, udává, o kolik lépe jedna hypotéza vysvětluje
data než druhá. **Efektivní velikost vzorku**
$n_\mathrm{eff} = (\sum w)^2/\sum w^2$ udává, kolik členů nese výsledek. Výsledek opřený
o hrstku členů je náznak, nikoli měření. Metody popisuje dokument *Fyzika a metody*, kapitoly
Ansámbly a inference a Fitování pozorovaných dat.

Všechna čísla v této části pocházejí z výstupních souborů uvedených v oddílu
[Data a reprodukovatelnost](#data). Obrázky vytváří `tools/part4_figures.py`.
