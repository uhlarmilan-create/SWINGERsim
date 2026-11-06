"""User documentation is kept in English and Czech; these tests keep the two in step."""

import re
from pathlib import Path

import pytest

from swingersim.config import (
    ActivitySpec,
    CommonEnvelopeSpec,
    Config,
    ContactSpec,
    CoordinatesSpec,
    DiscSpec,
    Distribution,
    EncountersSpec,
    EnsembleSpec,
    EphemerisSpec,
    ExtraModeSpec,
    FitExcludeSpec,
    FitOutputSpec,
    FitSpec,
    GalacticTideSpec,
    GasDamping,
    GasPlane,
    GasSpec,
    MagneticBrakingSpec,
    MassLossSpec,
    MassTransferSpec,
    ObservedValue,
    OrbitSpec,
    Physics,
    PulsationSpec,
    RunSpec,
    Stage1Spec,
    Stage2Spec,
    StarSpec,
    TertiaryTransferSpec,
)
from swingersim.scatter import ScatterBody, ScatterConfig, ScatterSpec

DOCS = Path(__file__).resolve().parents[1] / "docs"
LANGS = ("en", "cs")
DOCUMENTS = ("guide", "methods", "part1", "part2", "part3", "part4", "outlook")
CONFIG_CHAPTER = "guide/09_configuration.md"

pytestmark = pytest.mark.skipif(not (DOCS / "en" / "guide").is_dir(), reason="docs not present")


def _pages(lang):
    return sorted(
        str(p.relative_to(DOCS / lang)) for d in DOCUMENTS for p in (DOCS / lang / d).glob("*.md")
    )


def test_same_pages_in_both_languages():
    assert _pages("en") == _pages("cs")
    assert _pages("en"), "no documentation pages found"


@pytest.mark.parametrize("page", _pages("en") if (DOCS / "en").is_dir() else [])
def test_same_structure(page):
    """Same number of headings, tables, code blocks, figures and display formulas in both
    languages."""

    def shape(text):
        return (
            len(re.findall(r"^#+ ", text, re.MULTILINE)),
            len(re.findall(r"^\| ?(-|:-)", text, re.MULTILINE)),
            text.count("```"),
            text.count("!["),
            text.count("$$"),
        )

    en = (DOCS / "en" / page).read_text(encoding="utf-8")
    cs = (DOCS / "cs" / page).read_text(encoding="utf-8")
    assert shape(en) == shape(cs)


@pytest.mark.parametrize("lang", LANGS)
def test_every_config_key_documented(lang):
    text = (DOCS / lang / CONFIG_CHAPTER).read_text(encoding="utf-8")
    code = set()  # every word inside `...` and inside the YAML examples
    for span in re.findall(r"`([^`]+)`", text) + re.findall(r"```(.*?)```", text, re.DOTALL):
        code |= set(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", span))
    missing = []
    for model in (
        Config,
        OrbitSpec,
        StarSpec,
        PulsationSpec,
        ExtraModeSpec,
        Physics,
        RunSpec,
        EnsembleSpec,
        Distribution,
        ObservedValue,
        GasSpec,
        GasDamping,
        GasPlane,
        MassLossSpec,
        MassTransferSpec,
        CommonEnvelopeSpec,
        MagneticBrakingSpec,
        GalacticTideSpec,
        EncountersSpec,
        ActivitySpec,
        FitSpec,
        EphemerisSpec,
        FitExcludeSpec,
        CoordinatesSpec,
        Stage1Spec,
        Stage2Spec,
        FitOutputSpec,
        ScatterConfig,
        ScatterSpec,
        ScatterBody,
        DiscSpec,
        TertiaryTransferSpec,
        ContactSpec,
    ):
        for name in model.model_fields:
            if name not in code:
                missing.append(f"{model.__name__}.{name}")
    assert not missing, f"missing in {lang}/{CONFIG_CHAPTER}: {missing}"


@pytest.mark.parametrize("lang", LANGS)
def test_figures_exist(lang):
    for d in DOCUMENTS:
        for page in (DOCS / lang / d).glob("*.md"):
            for ref in re.findall(r"!\[[^\]]*\]\(([^)\s]+)", page.read_text(encoding="utf-8")):
                if ref.startswith("#"):
                    continue
                assert (page.parent / ref).resolve().exists(), f"{d}/{page.name}: missing {ref}"


@pytest.mark.parametrize("lang", LANGS)
def test_no_version_numbers(lang):
    for page in _pages(lang):
        text = (DOCS / lang / page).read_text(encoding="utf-8")
        assert not re.search(r"SWINGERsim \d+\.\d+|\bv0\.\d+|quadsim", text), page
