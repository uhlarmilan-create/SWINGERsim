# Jsou aposteriorní rozdělení zkonvergovaná? MCMC s 1000 kroky {#mcmc}

## Otázka {#mcmc-question}

Nejistoty uváděné v částech II a III pocházejí ze vzorkování aposteriorní hustoty pravděpodobnosti
fitovaných parametrů metodou Monte Carlo s Markovovými řetězci (MCMC). Používáme afinně invariantní
ansámblový vzorkovač emcee (Foreman-Mackey et al. 2013); řetězec fitu popisuje dokument *Fyzika
a metody*, kapitola Fitování pozorovaných dat. Výpočet MCMC je množina **chodců**, kteří se pohybují
prostorem parametrů. Každý **krok** stojí jednu úplnou N-částicovou integraci na chodce. Po
**zahoření** (burn-in) vzorkují polohy chodců aposteriorní rozdělení.

Polohy jednoho chodce jsou mezi kroky korelované. **Integrovaný autokorelační čas** $\tau$ je počet
kroků, po nichž chodec „zapomene“, kde byl. Řetězec o $N$ krocích obsahuje zhruba $N/\tau$
nezávislých vzorků na chodce. Foreman-Mackey et al. (2013) doporučují běh „a few (say 10)
autocorrelation times“; dokumentace emcee doporučuje nejméně 50 τ, než je spolehlivý i samotný odhad
$\tau$. SWINGERsim zahazuje první polovinu řetězce jako zahoření, měří $\tau$ na druhé polovině
a varuje, je-li ponechaná polovina kratší než 10 τ.

První fity používaly 400 kroků, pro něž emcee hlásil $\tau \approx 21$ (19 $\tau$). Abychom ověřili,
zda jsou uváděné nejistoty spolehlivé, přepočítali jsme tři fity s 1000 kroky: jednu soustavu 2+2
(BU CMi), jednu soustavu 3+1 (TIC 120362137) a jednu soustavu 2+1 s proměnným sklonem (V907 Sco).

## Výsledky {#mcmc-results}

**Tabulka 7.** Výpočty s 1000 kroky. Sloupce:

- *doba výpočtu*: skutečný čas na 12 jádrech;
- *$\tau$*: průměrný integrovaný autokorelační čas přes parametry, v krocích;
- *kroky/$\tau$*: délka řetězce v jednotkách $\tau$;
- *poměr chyb*: nejistoty 1σ z 1000 kroků dělené nejistotami ze 400 kroků, jako medián přes všechny
  parametry a rozsah;
- *největší posun*: největší změna mediánu mezi oběma výpočty, v jednotkách jeho nejistoty.

| soustava | doba výpočtu | $\tau$ [kroky] | kroky/$\tau$ | poměr chyb (medián; rozsah) | největší posun |
|:--|--:|--:|--:|:--|:--|
| TIC 120362137 | 10,5 h | 53,3 | 18,7 | 1,05; 0,86–2,02 | < 0,4σ |
| V907 Sco | 19,3 h | 53,9 | 18,6 | 1,10; 0,86–1,53 | < 0,8σ |
| BU CMi | 19,4 h | 52,3 | 19,1 | 1,03; 0,87–1,36 | < 0,25σ |

Vystupují tři věci.

1. **Autokorelační čas byl podceněn.** Při prodloužení řetězce na 2,5násobek vzrostl z asi 21 na asi
   53 kroků. To je známé chování krátkých řetězců: $\tau$ nelze změřit z řetězce kratšího než několik
   desítek $\tau$ a odhad s délkou řetězce roste (Foreman-Mackey et al. 2013; Sokal 1997). Řetězce
   jsou proto celkem asi 19 $\tau$ dlouhé a polovina ponechaná po zahoření má asi 9,4 $\tau$. To je
   těsně pod 10 τ podle Foreman-Mackeyho et al. a hluboko pod 50 τ podle dokumentace emcee.
2. **Výsledky se nepohnuly.** Žádný medián se neposunul o více než 0,8σ a většina o mnohem méně.
   Typická nejistota vzrostla o 3–10 %. U několika slabě omezených parametrů se širokými,
   negaussovskými aposteriorními rozděleními vzrostla až na dvojnásobek (TIC 120362137), resp.
   1,5násobek (V907 Sco).
3. **Nejistoty z MCMC a z linearizace souhlasí.** Kovarianční matice fitu metodou nejmenších čtverců
   („lineární“ chyby) souhlasí se šířkami z MCMC v mediánu asi na 10 %. Obojí popisuje totéž dobře
   definované minimum.

## Co to znamená pro uváděná čísla {#mcmc-meaning}

- **Mediány** uváděné v částech II a III jsou robustní. Zdvojnásobení délky řetězce by je
  neposunulo o více než zlomek jejich nejistoty.
- **Nejistoty 1σ** jsou spolehlivé na několik desítek procent. Jsou spíše mírně podhodnocené než
  nadhodnocené, typicky asi o 5–10 % a více u slabě omezených elementů.
- **Okraje rozdělení** (rozsahy 95 %) jsou méně jisté než šířky 1σ. Při asi 9–10 nezávislých
  vzorcích na chodce stojí na menším počtu vzorků. Nejvíce to má vliv na data budoucích událostí
  uváděná jako rozsahy 95 %, například konec zákrytů V907 Sco v letech 2052,5–2056,8. Šířka takových
  rozsahů je sama nejistá, více než šířka intervalů 1σ.
- Deset τ v ponechané polovině by vyžadovalo asi 1100 kroků, těsně nad použitými 1000. Padesát τ
  podle dokumentace emcee by vyžadovalo asi 5300 kroků, tj. 55–100 h na soustavu na 12 jádrech.
  Srovnání 400 a 1000 kroků ukazuje, že delší řetězce zpřesňují šířky, ale nemění závěry. Uvádíme
  výsledky s 1000 kroky.

Hodnoty v částech II a III jsou hodnotami z výpočtů s 1000 kroky.
