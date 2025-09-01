"""Unit tests for normalization utilities."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

import pytest

from mylib.transforms import normalize_name


def test_isotope_flag() -> None:
    res = normalize_name("[3H] 8 - oh dpat")
    assert res["flag_isotope"] is True
    assert res["search_name"] == "8-oh dpat"


def test_biotin_and_peptide() -> None:
    res = normalize_name("biotinylated peptide")
    assert res["flag_biotin"] is True
    assert res["category"] == "peptide"
    assert res["search_name"] == "peptide"


def test_sequence_detection() -> None:
    res = normalize_name("Ala-Gly-Ser")
    assert res["category"] == "peptide"
    assert res["peptide_info"]["type"] == "sequence_like"



def test_protective_group_sequence() -> None:
    """Sequences with protective termini are still classified as peptides."""

    res = normalize_name("H-Ala-Gly-OH")
    assert res["category"] == "peptide"
    assert res["peptide_info"]["type"] == "sequence_like"



def test_noise_and_concentration_removal() -> None:
    res = normalize_name("Sample solution 10 mM")
    assert "solution" in res["flags"].get("noise", [])
    assert res["search_name"] == "sample"


def test_parenthetical_noise_two_pass() -> None:
    res = normalize_name("histamine (solution, 10 mM in PBS)")
    noise_flags = res["flags"].get("noise", [])
    assert res["search_name"] == "histamine"
    assert "solution" in noise_flags
    assert "PBS" in noise_flags
    assert any("10 mM" in t for t in noise_flags)


def test_empty_after_clean_status_flag() -> None:
    """Names that empty out after cleaning fallback and set status."""

    res = normalize_name("solution 10 mM")
    assert res["status"] == "empty_after_clean"
    assert res["flag_empty_after_clean"] is True
    assert res["normalized_name"] == "solution 10 mm"
    assert res["search_name"] == "solution 10 mm"


@pytest.mark.parametrize("connector", ["-", "/", "+", ":"])
def test_spacing_compaction_after_flag_removal(connector: str) -> None:
    """Removal of flagged tokens should not leave spaces around connectors."""

    res = normalize_name(f"a FITC {connector} b")
    assert res["search_name"] == f"a{connector}b"


def test_repeated_connector_collapses() -> None:
    res = normalize_name("a - biotin - b")
    assert res["search_name"] == "a-b"




def test_spacing_for_comma_and_decimal() -> None:
    """Spaces around commas and decimals are compacted."""

    res = normalize_name("N , N-dimethyl 1 . 5")
    assert res["search_name"] == "n,n-dimethyl 1.5"





def test_salt_tokens_removed_and_logged() -> None:
    res = normalize_name("histamine hydrochloride")
    assert res["search_name"] == "histamine"
    assert res["flags"].get("salt") == ["hydrochloride"]


def test_mineral_acid_tokens() -> None:
    res = normalize_name("dextromethorphan HBr")
    assert res["search_name"] == "dextromethorphan"
    assert res["flags"].get("salt") == ["HBr"]


def test_salt_and_hydrate_combination() -> None:
    res = normalize_name("metformin hydrochloride monohydrate")
    assert res["search_name"] == "metformin"
    assert res["flags"].get("salt") == ["hydrochloride"]
    assert res["flags"].get("hydrate") == ["monohydrate"]


@pytest.mark.parametrize(
    "token",
    [
        "dihydrate",
        "trihydrate",
        "tetrahydrate",
        "pentahydrate",
        "anhydrous",
    ],
)
def test_various_hydrate_tokens(token: str) -> None:
    res = normalize_name(f"glucose {token}")
    assert res["search_name"] == "glucose"
    assert res["flags"].get("hydrate") == [token]


def test_removed_tokens_flat() -> None:
    text = "Alexa Fluor 488 [3H] histamine hydrochloride"
    res = normalize_name(text)
    assert (
        res["removed_tokens_flat"]
        == "fluorophore:Alexa Fluor 488|isotope:[3H]|salt:hydrochloride"
    )


def test_removed_tokens_flat_empty() -> None:
    res = normalize_name("aspirin")
    assert res["removed_tokens_flat"] == ""


def test_search_name_defaults_to_normalized() -> None:
    """search_name should equal normalized_name when no override reason."""

    res = normalize_name("Histamine")
    assert res["search_name"] == res["normalized_name"]
    assert res["search_override_reason"] == ""


@pytest.mark.parametrize(
    "text, expected, tokens",
    [
        ("[3H]-histamine", "histamine", ["[3H]"]),
        ("d5-amphetamine", "amphetamine", ["d5"]),
        ("U-13C glucose", "glucose", ["U-13C"]),
        ("d5 [3H] amphetamine", "amphetamine", ["d5", "[3H]"]),
    ],
)
def test_isotope_variants(text: str, expected: str, tokens: list[str]) -> None:
    """Isotope markers are stripped and logged across representations."""

    res = normalize_name(text)
    assert res["search_name"] == expected
    assert res["flags"].get("isotope") == tokens


@pytest.mark.parametrize(
    "text, expected, token",
    [
        ("poly-Glu:Tyr Alexa Fluor 488", "poly-glu:tyr", "Alexa Fluor 488"),
        ("HiLyte Fluor 555 peptide", "peptide", "HiLyte Fluor 555"),
        ("DyLight-650 antibody", "antibody", "DyLight-650"),
        ("peptide CF568", "peptide", "CF568"),
        ("Janelia Fluor 549 ligand", "ligand", "Janelia Fluor 549"),
        ("BODIPY-581/591 conjugate", "conjugate", "BODIPY-581/591"),
    ],
)
def test_expanded_fluorophore_tokens(text: str, expected: str, token: str) -> None:
    """Expanded fluorophore patterns are removed prior to other processing."""

    res = normalize_name(text)
    assert res["search_name"] == expected
    assert res["flags"].get("fluorophore") == [token]


@pytest.mark.parametrize(
    "text, composition",
    [
        ("poly-Glu:Tyr", "glu:tyr"),
        ("poly (Glu, Tyr)", "glu:tyr"),
        ("poly Glu Tyr", "glu:tyr"),
    ],
)
def test_poly_peptide_detection(text: str, composition: str) -> None:
    res = normalize_name(text)
    assert res["category"] == "peptide"
    assert res["peptide_info"] == {"type": "polymer", "composition": composition}


def test_polymer_non_peptide() -> None:
    res = normalize_name("polymer support resin")
    assert res["category"] == "small_molecule"

