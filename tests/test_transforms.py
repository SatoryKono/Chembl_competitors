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


def test_noise_and_concentration_removal() -> None:
    res = normalize_name("Sample solution 10 mM")
    assert "solution" in res["flags"].get("noise", [])
    assert res["search_name"] == "sample"


@pytest.mark.parametrize("connector", ["-", "/", "+", ":"])
def test_spacing_compaction_after_flag_removal(connector: str) -> None:
    """Removal of flagged tokens should not leave spaces around connectors."""

    res = normalize_name(f"a FITC {connector} b")
    assert res["search_name"] == f"a{connector}b"


def test_repeated_connector_collapses() -> None:
    res = normalize_name("a - biotin - b")
    assert res["search_name"] == "a-b"


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
