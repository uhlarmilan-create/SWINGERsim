# Diskuse

## Co činí tyto výsledky spolehlivými a kde jsou slabší

- **V907 Sco** spočívá na nejlepších datech: sklonech ze čtyř nezávislých zdrojů, rychlostech ze dvou epoch
  vzdálených 40 let, přesných minimech z TESS a aposteriorním rozdělení vzorkovaném metodou MCMC. Výsledek
  SWINGERsimu souhlasí s publikovaným v rámci chyb a SWINGERsim reprodukuje publikované předpovědi, dostane-li
  publikované parametry. Orientace drah na obloze (uzly) v práci uvedena není; do pozorovatelných veličin
  vstupuje jen rozdíl uzlů.
- **QX Cas** spočívá na dvou sklonech, jedné horní mezi a nepřítomnosti zákrytů. Hlavním výsledkem je proto
  ansámbl a fit je příkladem. Výpočet hloubek zákrytů ukázal, že fit splňuje 30. léta jen geometricky, ne
  pozorovatelně, a že staré hloubky vyžadují třetí světlo a jiné parametry hvězd. Detekce dvojhvězdy v TESS
  je nezávislou kontrolou amplitudy modelu, zatím ne jeho dráhy.
- **AH Cep** spočívá, pokud jde o třetí těleso, jen na efektu světelného času. Výsledek o třetím světle
  závisí na modelu světelné křivky.

## Poučení: „model zakrývá“ neznamená „zákryty byly vidět“

U všech tří soustav byla data přechodů nejprve srovnávána jako geometrická překročení meze zákrytů. U V907
Sco to stačilo, protože její zákryty se po překročení rychle prohlubují a modelové hloubky v pozorovaných
datech přechodů byly ověřeny (tab. 2). U QX Cas totéž srovnání skrylo skutečný problém: geometrický zákryt
v roce 1933 hluboký 0,0006 mag není pro fotografickou desku žádný zákryt. Kdykoli pozorování říká „zákryty
byly vidět“, musí model dát hloubku nad mezí zjistitelnosti tohoto pozorování, s uvedeným třetím světlem.
Taková pozorování nyní zadáváme v této podobě.

## Model světelné křivky

U všech tří soustav byly syntetické světelné křivky SWINGERsimu srovnány s TESS. Bez fitování reprodukují
zákryty V907 Sco s přesností 5 % v hloubce a amplitudu tlukotu srdce QX Cas. Tam, kde se hloubky fitují,
jsou zákryty SWINGERsimu pro dané parametry systematicky asi o 5 % mělčí než zákryty z PHOEBE (V907 Sco,
AH Cep). To posouvá nafitované sklony nahoru o zlomek stupně u V907 Sco a pravděpodobně o více u mělkých
částečných zákrytů AH Cep. Výsledky, které spočívají na hloubkách zákrytů (fit světelné křivky V907 Sco,
sklon a třetí světlo AH Cep, hloubky a zjistitelnost zákrytů QX Cas, tvar signálu QX Cas), je třeba číst
s tímto vědomím. Dynamické výsledky (okamžiky zákrytů, rychlosti, geometrická data přechodů) na něm
nezávisejí.

## Chyby

Chyby jsou lineární chyby ze zakřivení χ², kromě V907 Sco, kde se uvádějí chyby z aposteriorního rozdělení
MCMC; u V907 Sco se obě shodují v rámci 50 %. Řetězec MCMC pokrývá celkem asi 19 autokorelačních časů
a asi 9 po zahození poloviny na zahoření (viz Data a metody); řetězec o 400 krocích dal chyby, které se od
těchto liší v mediánu o 10 %, takže je považujeme za správné na několik desítek procent. U silně
korelovaných parametrů (například hmotnost a sklon dráhy třetího tělesa AH Cep) lineární chyby skutečný
rozsah podhodnocují.
