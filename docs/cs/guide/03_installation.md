# Instalace {#installation}

SWINGERsim je balík Pythonu, který běží v Linuxu a macOS (ve Windows prostřednictvím subsystému
Windows Subsystem for Linux, WSL). Tato kapitola uvádí, co potřebuje, jak ho nainstalovat a instalaci
ověřit, jaké má volitelné části, jak funguje zpracovatel úloh pro dlouhé výpočty a jak sestavit tuto
dokumentaci.

## Požadavky {#requirements}

| Položka | Požadavek | Proč |
| --- | --- | --- |
| Python | 3.11 nebo novější | jazyk balíku |
| překladač C | `gcc` nebo `clang` | REBOUND a REBOUNDx, N-částicové knihovny, se při instalaci překládají, pokud pip nenajde hotový balík |
| operační systém | Linux nebo macOS; Windows přes WSL | zpracovatel úloh používá skupiny procesů a signály systémů POSIX |
| disk | asi 100 MB pro vývojové dráhy hvězd a tabulky balíku, k tomu volitelné knihovny (samotná skupina `sbi` asi 1 GB) | data jsou součástí balíku |

**Tabulka 7.** Požadavky. Levý sloupec uvádí položku, prostřední požadavek, pravý sloupec, proč je
potřeba. Paměť závisí na úloze: paralelní ansámbl potřebuje paměť pro každý proces, a proto
zpracovatel úloh hlídá volnou paměť (viz níže).

## Postup {#install-steps}

```bash
git clone https://github.com/uhlarmilan-create/SWINGERsim.git
cd SWINGERsim
./install.sh
source .venv/bin/activate
```

Instalační skript `install.sh` najde Python 3.11 nebo novější, ověří, že jsou přítomny vývojové
dráhy hvězd a tabulky balíku, vytvoří ve složce `.venv` izolované prostředí Pythonu (aby se balíky
SWINGERsimu nemísily se zbytkem systému), nainstaluje SWINGERsim se skupinami `all` a `dev`
z tabulky 8, provede zkušební běh a spustí rychlé srovnávací testy. Nikdy nepoužívá `sudo`: když
chybí systémový balík (Python, jeho modul `venv`, překladač C), zastaví se a vypíše příkaz, který ho
nainstaluje. Jeho volby:

| Volba | Účinek |
| --- | --- |
| `--minimal` | jen jádro a sada testů (skupina `dev`) |
| `--with-sbi` | navíc skupina `sbi` (PyTorch, asi 1 GB) |
| `--with-phoebe` | navíc skupina `phoebe` |
| `--venv DIR` | prostředí v `DIR` místo `.venv` |
| `--python PY` | tento interpret Pythonu |
| `--no-check` | na konci žádný zkušební běh ani srovnávací testy |

**Tabulka 7a.** Volby skriptu `install.sh` (vypíše je `./install.sh --help`). Lze je kombinovat,
např. `./install.sh --with-sbi --no-check`.

Prostředí aktivujte v každém novém terminálu příkazem `source .venv/bin/activate`; poté je k
dispozici příkaz `swingersim`. **Windows:** jednou nainstalujte WSL (PowerShell jako správce:
`wsl --install`, restart, jednou otevřít *Ubuntu* a vytvořit uživatele) a výše uvedené kroky
proveďte v terminálu Ubuntu; případně `powershell -ExecutionPolicy Bypass -File install.ps1` ověří
WSL a spustí v něm `install.sh`. Webová aplikace spuštěná ve WSL se otevře v prohlížeči Windows na
vypsané adrese.

Totéž ručně:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[all,dev]"
swingersim --help
```

`pip install -e` nainstaluje balík v upravitelném režimu: přímo se používá zdrojový kód
v `src/swingersim/`, takže se najdou konfigurace, data a výstupy ve složce repozitáře. Názvy
v hranatých závorkách jsou volitelné skupiny závislostí (tabulka 8). Poslední příkaz vypíše
podpříkazy; pokud je vypíše, instalace funguje.

| Skupina | Balíky | Co umožňuje |
| --- | --- | --- |
| (vždy) | numpy, scipy, pydantic, pyyaml, h5py, rebound, reboundx | běhy, veškerou fyziku, výstup HDF5 |
| `dev` | pytest, ruff | sadu testů |
| `all` | skupiny `viz`, `ensemble`, `fast`, `emulator`, `fit` | vše kromě `sbi` a `phoebe`; instaluje ji `install.sh` |
| `viz` | plotly, kaleido, matplotlib | interaktivní 3D prohlížeč, grafy ve zprávách HTML, obrázky PNG (matplotlib; kaleido s prohlížečem Chrome pro soubory PNG z grafů plotly); webová aplikace plotly použije, je-li k dispozici, jinak načte plotly.js z internetu |
| `fit` | emcee, corner, astropy | chyby proložení metodou MCMC, korelační grafy (corner plots), převod HJD/JD (UTC) na BJD_TDB |
| `fast` | numba | přeložené sekulární rovnice, asi stokrát rychlejší sekulární běhy a ansámbly |
| `ensemble` | pandas, pyarrow | souhrny ansámblů také jako soubory Parquet |
| `emulator` | scikit-learn | neuronové emulátory (`emulator-train`, `emulator-report`) a klasifikátory stability popsané níže |
| `sbi` | sbi, torch | simulační inferenci (`swingersim sbi`); asi 1 GB |
| `phoebe` | phoebe | PHOEBE 2 jako výpočetní modul světelných křivek vnitřních dvojhvězd (`--engine phoebe`) |

**Tabulka 8.** Volitelné skupiny závislostí. Skupina se instaluje příkazem
`pip install -e ".[<group>]"`; více skupin se odděluje čárkami. Pravý sloupec říká, které příkazy
nebo volby skupinu potřebují; bez ní se zastaví se zprávou, která jmenuje chybějící balík.

## Ověření instalace {#install-check}

```bash
pytest -m "not slow and not validation"     # unit tests, a few minutes
swingersim validate                          # fast benchmarks against published results
swingersim run configs/examples/example_2p2.yaml -o output/check/example.h5
```

Jednotkové testy kontrolují převody mezi elementy a polohami, hierarchickou inicializaci, zachování
energie, relativistické stáčení periastra, schéma konfigurace, webovou aplikaci a dokumentaci.
`swingersim validate` reprodukuje sadu publikovaných výsledků a pro každý vypíše PASS (vyhovuje),
EXPLAINED (rozdíl s doloženou příčinou) nebo FAIL (nevyhovuje); srovnávací testy a jejich verdikty
rozebírá dokument *Testování, část I*. Poslední příkaz je prvním během z kapitoly
[Rychlý start](#quick-start).

## Volitelně: klasifikátory stability Vynatheya et al. (2023) {#install-vynatheya}

Klasifikátory stability čtyřhvězd 2+2 a 3+1 založené na strojovém učení (Vynatheya et al. 2023)
publikovali jejich autoři a nejsou součástí SWINGERsimu. Chcete-li je použít, stáhněte jejich
repozitář do `data/`:

```bash
git clone https://github.com/pavanvyn/quadruple-stability data/quadruple-stability
pip install -e ".[emulator]"
```

`swingersim stability config.yaml` pak uvede jejich verdikt pro každou dílčí čtyřhvězdu. Sítě byly
uloženy ve scikit-learn 1.2.2; novější verze je načtou a varování o verzi je potlačeno. Jiné umístění
lze zadat proměnnou prostředí `SWINGERSIM_VYNATHEYA_DIR`.

## Volitelně: MSE pro nezávislá srovnání {#install-mse}

Kód Multiple Stellar Evolution, MSE (Hamers et al. 2021), je nezávislý kód pro hierarchické
vícenásobné hvězdy, který posloužil ke křížové kontrole SWINGERsimu (*Testování, část I*). Pro běh
SWINGERsimu potřeba není. Chcete-li provádět vlastní srovnání, nainstalujte ho z jeho repozitáře;
potřebuje překladače C, C++ a Fortranu:

```bash
sudo apt install gfortran
git clone https://github.com/hamers/mse ~/mse
cd ~/mse && make
```

## Volitelně: MESA pro vlastní mřížky modelů hvězd {#install-mesa}

SWINGERsim umí spočítat vlastní mřížky modelů hvězd programem MESA (Modules for Experiments in Stellar
Astrophysics; Paxton et al. 2011; Jermyn et al. 2023): stavbu ve fázi před hlavní posloupností
a stavbu obálky, $k_2$ a disipaci inerciálních vln z úplných profilů hvězd. To je potřeba jen tehdy,
když vestavěné vývojové dráhy MIST a Claretovy tabulky danou hvězdu nepokrývají (neobvyklá
metalicita, zvláštní nastavení fyziky). MESA se jednou nainstaluje ručně; SWINGERsim ji pak spouští
(i prostřednictvím zpracovatele úloh).

Požadavky: Linux x86_64, nejméně 8 GB paměti, asi 30 GB volného místa na disku (MESA asi 20 GB, její
sada vývojových nástrojů asi 2 GB, mřížka několik GB) a jedna až dvě hodiny času.

1. Systémové balíky (sada vývojových nástrojů MESA, SDK, si přináší vlastní překladače):

   ```bash
   sudo apt install binutils make perl libx11-6 libx11-dev zlib1g zlib1g-dev tcsh unzip python3
   ```

2. Stáhněte MESA SDK pro Linux x86_64 z
   <http://user.astro.wisc.edu/~townsend/static.php?ref=mesasdk> a rozbalte ho do `~/mesa`.
3. Stáhněte vydání MESA (záznam na Zenodo <https://doi.org/10.5281/zenodo.2602941> nebo
   <https://github.com/MESAHub/mesa/releases>) a rozbalte ho do `~/mesa`. SDK musí mít verzi, kterou
   dané vydání MESA vyžaduje.
4. Přidejte do `~/.bashrc` následující řádky a otevřete nový terminál:

   ```bash
   export MESASDK_ROOT=$HOME/mesa/mesasdk
   source $MESASDK_ROOT/bin/mesasdk_init.sh
   export MESA_DIR=$HOME/mesa/mesa-<release>
   export OMP_NUM_THREADS=2
   ```

5. Přeložte MESA (30–60 minut): `cd $MESA_DIR && ./install`.
6. Otestujte se SWINGERsimem – jedna vývojová dráha hvězdy 1 M☉ od fáze před hlavní posloupností do
   konce hlavní posloupnosti:

   ```bash
   swingersim mesa-grid configs/mesa/lowmass_pms_ms.yaml --test
   ```

7. Spočítejte celou mřížku a zredukujte ji (hodiny; oba kroky mohou běžet jako úlohy zpracovatele):

   ```bash
   swingersim mesa-grid configs/mesa/lowmass_pms_ms.yaml -j 6
   swingersim mesa-reduce configs/mesa/lowmass_pms_ms.yaml
   ```

Zredukovaná mřížka je `mesa/lowmass_pms_ms/lowmass_pms_ms.npz`; konfigurace ji použije zápisem
`stellar_model.structure: mesa:lowmass_pms_ms`.

## Zpracovatel úloh pro dlouhé výpočty {#install-worker}

Velké ansámbly (desítky tisíc členů, N-částicové ansámbly), proložení mnoha volných parametrů
a mřížky MESA trvají hodiny až dny. Nejlépe se počítají na počítači s mnoha jádry, který zůstává
zapnutý. `swingersim worker` z takového počítače udělá spouštěč úloh: sleduje složku, spustí každý
soubor úlohy, který se v ní objeví, a výsledky zapíše do `output/`.

```bash
cd swingersim
source .venv/bin/activate
swingersim worker                  # watches jobs/, Ctrl+C stops it
swingersim worker --cores 8        # use at most 8 cores for all jobs together
```

Zpracovatel vytvoří složky z tabulky 9 a každých několik sekund zapisuje soubor `jobs/worker.json`
s údaji o své činnosti (název počítače, počet jader, verze, s nimiž byl spuštěn, aktuální úloha, čas
poslední kontroly). Úlohy, které potřebují jedno jádro (běh, světelná křivka, okamžiky zákrytů),
běží vedle sebe až do počtu jader; úloha, která potřebuje všechna jádra (ansámbl, proložení nebo
mapa chaosu s `processes: 0`, `validate --all`), počká na ostatní a běží sama.

| Složka nebo soubor | Obsah |
| --- | --- |
| `jobs/drafts/` | konfigurace napsané ve webové aplikaci (zpracovatel je nespouští) |
| `jobs/submitted/` | zmrazené kopie konfigurací odeslaných úloh |
| `jobs/inbox/` | `<id>.yaml`: úlohy ve frontě, brané v pořadí jejich názvů; `<id>.yaml.hold`: úlohy čekající na zařazení do fronty |
| `jobs/running/` | právě počítané úlohy, každá s `<id>.yaml`, svým výpisem na obrazovku `<id>.log` a svým průběhem `<id>.status.json` |
| `jobs/done/`, `jobs/failed/` | dokončené a neúspěšné úlohy se stejnými třemi soubory |
| `jobs/cancel/` | prázdný soubor `<id>` danou úlohu zastaví; její výstup v `output/` se smaže a její soubory se uchovají zde |
| `jobs/pause/` | prázdný soubor `<id>` pozastaví běžící úlohu; odstranění souboru ji nechá pokračovat přesně tam, kde se zastavila |
| `jobs/stop` | prázdný soubor zastaví zpracovatele; běžící úlohy se vrátí do fronty |
| `jobs/worker.json`, `jobs/worker.lock` | údaje o činnosti (heartbeat); zámek, který brání druhému zpracovateli na téže složce |

**Tabulka 9.** Složka úloh zpracovatele. Levý sloupec je složka nebo soubor uvnitř složky
SWINGERsimu, pravý sloupec to, co obsahuje. Webová aplikace používá stejné uspořádání.

Aby zpracovatel běžel i po odhlášení, spusťte ho v `tmux` nebo jako uživatelskou službu, např.
`systemd-run --user --unit=swingersim-worker swingersim worker --root ~/swingersim`. Klesne-li volná
paměť počítače pod 8 % celkové (nejméně 1,5 GB), zpracovatel zastaví úlohu, která používá nejvíc
paměti, a označí ji jako neúspěšnou s vysvětlením, aby počítač zůstal použitelný. Po aktualizaci
SWINGERsimu zpracovatele restartujte; `worker.json` ukazuje verzi, kterou spouští. Soubory úloh
a pravidla pro ně popisuje kapitola [Pracovní postupy](#wf-worker).

## Sestavení této dokumentace {#install-docs}

Dokumentace je napsána v Markdownu v `docs/en/` (anglicky) a `docs/cs/` (česky), jedna složka pro
každý dokument. Soubory PDF se sestavují programy pandoc a XeLaTeX:

```bash
sudo apt install pandoc texlive-xetex texlive-fonts-recommended texlive-lang-czechslovak
python tools/build_docs_pdf.py guide --lang en      # this guide: docs/pdf/SWINGERsim-guide-en.pdf
python tools/build_docs_pdf.py                      # every document in both languages
```
