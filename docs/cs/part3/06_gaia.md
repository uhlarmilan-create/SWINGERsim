# Nezávislé kontroly pomocí Gaia

Gaia měří polohu světelného středu hvězdy mnohokrát během tří let. Pohybuje-li se světelný střed po
dráze, řešení pro jednu hvězdu sedí hůře: renormalizovaná chyba jednotkové váhy (RUWE) roste nad asi 1,4
(práh navržený Lindegrenem 2018), astrometrický nadbytečný šum je významný, paralaxy z DR2 a DR3 (různá
časová období) se mohou lišit a v jasných případech dává katalog Non-Single Star astrometrickou dráhu. To
vše je nezávislé na všem, co bylo použito výše. Třetí hvězda o hmotnosti m₃ a podílu světla l₃ na dráze
o velikosti a kolem páru posouvá světelný střed asi o a (m₃/M − l₃), kde M je celková hmotnost.

**Tabulka 12.** Očekávaný pohyb světelného středu podle modelů této zprávy. Vzdálenosti jsou hrubé hodnoty
z jasnosti hvězd; nahradí je paralaxa z Gaia.

| | V907 Sco | QX Cas | AH Cep |
|---|---|---|---|
| vnější dráha | 142 d, a ≈ 0,98 au | 190–1030 d, a ≈ 1,6–4,8 au (fit: 288 d, 2,2 au) | 70 let, a ≈ 56 au |
| m₃/M − l₃ | 0,16 − 0,02 = 0,14 | 0,04–0,17 (fit: 0,28 − 0,13 = 0,15) | ≥ 0,19 − 0,07 ≈ 0,12 |
| světelný střed se pohybuje o | ≈ 0,14 au | ≈ 0,06–0,8 au (fit: 0,3 au) | ≈ 7 au za 70 let |
| ve vzdálenosti | 0,3–0,5 kpc: ≈ 0,3–0,5 mas | ~3,3 kpc (NGC 7790): 0,02–0,25 mas (fit: 0,1 mas) | ~1 kpc: ≈ 7 mas; vlastní pohyb se za 25 let změní o ~1 mas/rok |
| očekáváno v Gaia | RUWE > 1,4 pravděpodobně; dráha v NSS možná | na mezi | rozdíl vlastních pohybů Hipparcos–Gaia |

V907 Sco by měla třetí hvězdu ukázat zřetelně; QX Cas je na mezi toho, co Gaia vidí, takže normální RUWE
by ansámbl nevyloučila; AH Cep, dost jasná pro Hipparcos, by měla 70letou dráhu ukázat jako změnu svého
vlastního pohybu mezi Hipparcosem (1991) a Gaia (2016).

## Co ukazuje Gaia DR3 {#gaia-results}

V archivu Gaia (DR3, Gaia Collaboration et al. 2023, a pro srovnání DR2) jsme na polohách tří soustav
provedli dotaz nástrojem `tools/gaia_check.py`. Astrometrie DR3 je řešení popsané Lindegrenem et al.
(2021).

**Tabulka 13.** Astrometrie tří soustav z Gaia DR3. Sloupce:

- *G*: hvězdná velikost v pásmu Gaia;
- *paralaxa*: paralaxa z DR3 s chybou a odpovídající vzdálenost (1/paralaxa);
- *paralaxa DR2*: hodnota z DR2 a rozdíl DR3 − DR2 v jednotkách sloučené chyby;
- *RUWE*: hodnoty nad asi 1,4 ukazují na špatný fit pro jednu hvězdu;
- *nadb. šum*: astrometrický nadbytečný šum (AEN, dodatečný rozptyl potřebný k proložení poloh)
  a jeho významnost;
- *NSS*: typ řešení v katalogu Non-Single Star, pokud existuje.

| soustava | G | paralaxa [mas] (vzdálenost) | paralaxa DR2 (rozdíl) | RUWE | nadb. šum (význ.) | NSS |
|---|---|---|---|---|---|---|
| V907 Sco | 8,55 | 2,279 ± 0,033 (439 pc) | 2,337 ± 0,064 (−0,8σ) | 0,81 | 0,16 mas (11) | SB2, P = 3,978 d |
| QX Cas | 10,12 | 0,268 ± 0,019 (3,7 kpc) | 0,115 ± 0,031 (+4,2σ) | 1,18 | 0,17 mas (39) | – |
| AH Cep | 6,84 | 1,139 ± 0,024 (878 pc) | 1,102 ± 0,034 (+0,9σ) | 1,18 | 0,23 mas (52) | zákrytová dvojhvězda, P = 1,774737 d |

**Jak číst tabulku 13.** Každý řádek je jedna soustava; sloupce jsou vysvětleny výše. Žádná ze tří soustav
nemá RUWE nad 1,4. Nadbytečný šum je formálně významný u všech tří, ale u jasných hvězd je to spíše
pravidlo než výjimka: mezi pětiparametrovými řešeními Gaia EDR3 (stejná astrometrie jako DR3) s G = 9–12
má asi 97 % významný nadbytečný šum, s mediánem asi 0,11 mas (online dokumentace Gaia EDR3, odd. 4.5.2).
Velká část z něj pravděpodobně pochází z nedokonalé kalibrace jasných hvězd, ne z oběžného pohybu.
Významný nadbytečný šum 0,16–0,23 mas, jako zde, je proto sám o sobě slabým důkazem; robustnějším
ukazatelem je RUWE.

- **V907 Sco: rozpor.** Ve vzdálenosti 439 pc z Gaia posouvá náš model světelný střed na 142denní dráze
  asi o 0,14 au, tj. asi o 0,32 mas. Gaia během časového rozsahu DR3 pokryla asi sedm takových oběhů.
  Pohyb o 0,3 mas při G = 8,5 obvykle zvýší RUWE nad 1,4, přesto je 0,81. Nadbytečný šum (0,16 mas) má
  správný řád, sám o sobě však dráhu neprokazuje. Katalog NSS dává řešení dvojčarové spektroskopické
  dvojhvězdy (SB2) s P = 3,978 d. To je o 5 % více než perioda zákrytů 3,776 d, takže pravděpodobně
  odráží vnitřní pár viděný se špatně určenou periodou, ne 142denní dráhu. Možnými důvody nízké RUWE jsou
  menší pohyb světelného středu, než jaký modelujeme (více světla od třetí hvězdy nebo menší m₃/M), nebo
  vyrušení geometrií skenování. Samotná třetí hvězda je prokázána okamžiky zákrytů a radiálními
  rychlostmi (Zasche et al. 2023), takže ji to nezpochybňuje. Astrometrická dráha v příštím vydání Gaia
  (DR4, s epochovou astrometrií) je však přímým testem našeho m₃.
- **QX Cas: konzistentní, ne rozhodující.** RUWE je normální. Nadbytečný šum 0,17 mas má zhruba velikost
  pohybu světelného středu podle našeho fitu (0,1 mas ve vzdálenosti 3,3 kpc), ale u hvězdy s G = 10 ve
  hvězdokupě to důkazem není. Paralaxa se mezi DR2 a DR3 změnila o 4,2σ. Hodnota z DR2 (0,115 mas, tj.
  8,7 kpc) je neslučitelná s hvězdokupou NGC 7790 ve vzdálenosti asi 3,3 kpc, kdežto DR3 s ní souhlasí.
  Oběžný pohyb světelného středu může zkreslit paralaxu měřenou během kratšího časového úseku, takže je to
  náznak, ale stejně pravděpodobnou příčinou jsou systematické chyby DR2.
- **AH Cep: podle očekávání.** RUWE je normální, jak se čeká u vnější dráhy o periodě 70 let, kterou Gaia
  vidí jako téměř přímočarý pohyb. Katalog NSS rozpoznává zákrytovou dvojhvězdu s přesně periodou
  vnitřního páru (1,774737 d). Výše navržený test, změna vlastního pohybu mezi Hipparcosem a Gaia,
  využívá katalog zrychlení Hipparcos–Gaia (Brandt 2021). Zatím jsme jej neprovedli; je uveden
  v dokumentu *Otevřené otázky a další práce*.
