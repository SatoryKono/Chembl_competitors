"""Tests for build_compound_name_dictionary."""

from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))

from pubchem_library import build_compound_name_dictionary


def test_build_compound_name_dictionary_basic() -> None:
    df = pd.DataFrame(
        {
            "search_name": ["Aspirin", "Foo", "Ignored"],
            "pubchem_cid": ["1", "2", "3"],
            "canonical_smiles": ["", "", ""],
            "inchi": ["", "", ""],
            "inchi_key": ["AAAA-BBBB-N", "CCCC-DDDD-N", "EEEE-FFFF-N"],
            "molecular_formula": ["", "", ""],
            "molecular_weight": ["", "", ""],
            "iupac_name": ["Acetylsalicylic acid", "Foo", "Bar"],
            "synonyms": ["aspirin|asa", "", "unknown"],
        }
    )

    result = build_compound_name_dictionary(df)

    assert set(result.columns) == {
        "synonyms",
        "preferent_name",
        "pubchem_cid",
        "inchi_key",
        "merge_index",
        "reference_synonyms",
    }

    aspirin_rows = result[result["preferent_name"] == "Aspirin"]
    assert set(aspirin_rows["synonyms"]) == {
        "aspirin",
        "asa",
        "acetylsalicylic acid",
    }
    assert aspirin_rows["reference_synonyms"].iloc[0] == "aspirin|asa|acetylsalicylic acid"

    foo_rows = result[result["preferent_name"] == "Foo"]
    assert foo_rows["synonyms"].tolist() == ["foo"]
    assert set(result["preferent_name"]) == {"Aspirin", "Foo"}


def test_synonyms_are_unique() -> None:
    """Duplicate synonym strings should appear only once."""
    df = pd.DataFrame(
        {
            "search_name": ["A", "B"],
            "pubchem_cid": ["1", "2"],
            "canonical_smiles": ["", ""],
            "inchi": ["", ""],
            "inchi_key": ["AAAA-BBBB-N", "CCCC-DDDD-N"],
            "molecular_formula": ["", ""],
            "molecular_weight": ["", ""],
            "iupac_name": ["A", "B"],
            "synonyms": ["foo|shared", "bar|shared"],
        }
    )

    result = build_compound_name_dictionary(df)

    # Ensure that each synonym string appears only once across the dictionary
    assert not result["synonyms"].duplicated().any()
