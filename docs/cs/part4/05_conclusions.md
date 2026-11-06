# Závěry {#conclusions}

Část IV použila SWINGERsim na tři otázky, které literatura ponechává otevřené.

1. **CzeV343.**
   - Čtyřhvězda se stane dynamicky nestabilní jen tehdy, je-li donor dvojhvězdy A obnažen
     a dvojhvězda se rozšíří na asi 57–62 R☉. S předpisy pro přenos hmoty ve SWINGERsimu se to
     stane ve 14 % budoucností základního řešení (1,3–1,6 Gyr od současnosti) a v 1 % budoucností
     alternativního řešení; ve většině budoucností dvojhvězda A dříve splyne nebo se stane kontaktní
     dvojhvězdou. Jak často se A rozšíří, závisí na vývojových stopách: tytéž tahy se SSE místo MIST
     ji rozšíří v 82 % (COMPAS, postavený na týchž vzorcích: ≳ 85 %; Pejcha et al. 2022) a 75 %
     z nich se stane nestabilními. Vzdálenost složek v okamžiku počátku nestability a to, co následuje, na
     tom nezávisejí.
   - Blízko počátku nestability je pro takové soustavy bezpečným kritériem kritérium Mardlinga
     a Aarsetha: v 1200 N-částicových realizacích najde 98 % rozpadů do $4 \times 10^4$ let, kritérium
     Vynatheyi et al. (2022) 51 % a jejich neuronová síť 25 %.
   - CzeV343 není v rezonanci 3:2 mezi svými dvojhvězdami; leží asi 80 šířek rezonance od ní.
   - Po počátku nestability není soustava roztržena. S apsidální konstantou spočítanou pro
     obnaženého obra ($k_2 \approx 10^{-3}$) zvýší cykly ZLK buzené dvojhvězdou B excentricitu A
     na 0,7 a obr se se svým průvodcem srazí během 50–550 let ve všech
     600 integracích.
   - Jen nepravděpodobně hmotná obálka ($k_2 \gtrsim 0{,}01$) by ZLK potlačila a udržela soustavu
     neporušenou po dobu $10^4$ let.
   - Odpověď na otázku „kdy a jak“ tedy zní: rozšíří-li se dvojhvězda A (v našich modelech zřídka,
     v COMPAS často), nestabilita začne při 57–62 R☉, v našich modelech asi za 1,3–1,6 Gyr,
     a během staletí skončí splynutím uvnitř čtyřhvězdy způsobeným ZLK, nikoli vymrštěním.
2. **CzeV4315.**
   - Rychlá změna sklonu dvojhvězdy B nemůže pocházet od dvojhvězdy A. V uspořádání 2+2 porucha,
     která stáčí B, stáčí i A o 2–50°/rok, nebo se hvězdy soustavy srazí.
   - Třetí hvězda na 119denní dráze kolem B vysvětluje všechna pozorování s evidencí větší nejméně
     dvacetkrát (formálně $10^{13}$krát). Zasche et al. (2025) ponechali uspořádání otevřené se
     slabým náznakem stejným směrem.
   - Její hmotnost je degenerovaná s excentricitou a sklonem její dráhy. 74 % váhy leží nad 2 M☉,
     což lze ověřit spektroskopií nebo amplitudou efektu světelného času 119denního signálu.
   - Další fyzický průvodce ve vzdálenosti 1,84″ činí ze soustavy nejméně šestihvězdu.
3. **Konvergence fitů.**
   - Výpočty MCMC s 1000 kroky jsou asi 19 autokorelačních časů dlouhé; polovina ponechaná po
     zahoření má asi 9,4 τ.
   - Jejich mediány souhlasí s výpočty se 400 kroky v rámci 0,8σ a jejich nejistoty jsou typicky
     o 3–10 % větší.
   - Hodnoty uváděné v částech II a III jsou robustní. Jejich nejistoty 1σ jsou spolehlivé na několik
     desítek procent, jejich rozsahy 95 % o něco méně.

**Co tyto testy ukazují o programu.** Každá otázka vyžadovala souhru několika částí
SWINGERsimu:

- vývoj hvězd a přenos hmoty se sledováním stability (CzeV343, fáze 1);
- N-částicovou integraci s relativitou, slapy a rotací, sestavenou automaticky z vyvinutého
  stavu (fáze 2);
- ansámbly s vahami pozorování a evidencemi pro dvě různé hierarchie (CzeV4315);
- řetězec fitu s MCMC (oddíl 4).

Odpovědi závisejí na fyzikálních vstupech, které pozorování neurčují: na apsidální konstantě
obnaženého obra a na hmotnosti neviděné hvězdy. SWINGERsim tyto závislosti zviditelňuje
a kvantifikuje, místo aby je skrýval v jediném výpočtu.

**Otevřené body**, shrnuté v dokumentu *Otevřené otázky a další práce*:

- podrobný kód vývoje dvojhvězd se vstupy pro CzeV343, který rozhodne mezi chováním donoru
  dvojhvězdy A podle MIST a podle SSE;
- spektroskopie CzeV4315, která odhalí nebo vyloučí hmotnou třetí hvězdu kolem B, a změření
  amplitudy jejího 119denního signálu v okamžicích zákrytů;
- delší řetězce MCMC (≈ 1100 kroků pro 10 τ po zahoření, ≈ 5300 pro 50 τ) tam, kde záleží na
  okrajích aposteriorního rozdělení.

# Data a reprodukovatelnost {#data}

**Tabulka 8.** Vstupy a výstupy části IV. Každý řádek uvádí konfiguraci nebo nástroj (vstup)
a výstupní složku, kterou vytváří; cesty jsou relativní vůči adresáři SWINGERsimu.

| položka | vstup | výstup |
|:--|:--|:--|
| CzeV343, fáze 1 a nastavení fáze 2 | `tools/czev343_future.py`, `configs/special/czev343.yaml` | `output/special/czev343/stage1.csv`, `stage2/*.yaml` |
| CzeV343, fáze 2 | `swingersim ensemble output/special/czev343/stage2/<onset>_k2_<k2>.yaml` | `output/special/czev343/s2tides/` |
| apsidální konstanta obnaženého obra | `tools/k2_condensed.py` | tabulka 2 |
| CzeV4315, 2+2 | `configs/special/czev4315_2p2_nbody_ensemble.yaml` | `output/special/czev4315_2p2/` |
| CzeV4315, 2+(2+1) | `configs/squadra/czev4315_2p2p1_wide_ensemble.yaml` | `output/special/czev4315_2p2p1_wide/` |
| MCMC s 1000 kroky | konfigurace fitů z částí II a III s 1000 kroky | složky fitů z částí II a III |
| obrázky a čísla | `tools/part4_figures.py` | `docs/figures/part4/` |

Příkazy jsou:

```bash
python tools/czev343_future.py -o output/special/czev343          # fáze 1 + nastavení fáze 2
for f in output/special/czev343/stage2/*_k2_*.yaml; do
  swingersim ensemble "$f" -o "output/special/czev343/s2tides/$(basename "$f" .yaml)" -j 12
done
swingersim ensemble configs/special/czev4315_2p2_nbody_ensemble.yaml -o output/special/czev4315_2p2 -j 12
swingersim ensemble configs/squadra/czev4315_2p2p1_wide_ensemble.yaml -o output/special/czev4315_2p2p1_wide -j 12
python tools/k2_condensed.py
python tools/part4_figures.py
```
