"""Tests for PubChem lookup utilities."""

from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd
import pytest
import requests

sys.path.append(str(Path(__file__).resolve().parents[1]))

from mylib import annotate_pubchem_info
from mylib.pubchem import fetch_pubchem_cid, fetch_pubchem_record


class DummyResponse:
    """Minimal response stub for :mod:`requests` Session.get."""

    def __init__(self, *, text: str = "", json_data: dict | None = None, status_code: int = 200) -> None:
        self.text = text
        self._json = json_data or {}
        self.status_code = status_code

    def raise_for_status(self) -> None:
        if self.status_code >= 400 and self.status_code not in {400, 404}:
            raise requests.HTTPError(f"HTTP {self.status_code}")

    def json(self) -> dict:
        return self._json


# ---------------------------------------------------------------------------
# fetch_pubchem_cid tests
# ---------------------------------------------------------------------------

def test_fetch_pubchem_cid_short_name(monkeypatch: pytest.MonkeyPatch) -> None:
    sess = requests.Session()

    def fake_get(*args, **kwargs):  # pragma: no cover - should not be called
        raise AssertionError("get should not be called for short names")

    monkeypatch.setattr(sess, "get", fake_get)
    assert fetch_pubchem_cid("abcd", session=sess) == "compound name is too short"


def test_fetch_pubchem_cid_single(monkeypatch: pytest.MonkeyPatch) -> None:
    sess = requests.Session()
    monkeypatch.setattr(sess, "get", lambda *args, **kwargs: DummyResponse(text="123\n"))
    assert fetch_pubchem_cid("aspirin", session=sess) == "123"


def test_fetch_pubchem_cid_multiple(monkeypatch: pytest.MonkeyPatch) -> None:
    sess = requests.Session()

    monkeypatch.setattr(sess, "get", lambda *args, **kwargs: DummyResponse(text="1\n2\n"))

    assert fetch_pubchem_cid("foo bar", session=sess) == "multiply"


def test_fetch_pubchem_cid_unknown(monkeypatch: pytest.MonkeyPatch) -> None:
    sess = requests.Session()
    monkeypatch.setattr(sess, "get", lambda *args, **kwargs: DummyResponse(text="", status_code=404))
    assert fetch_pubchem_cid("unknowncompound", session=sess) == "unknown"


def test_fetch_pubchem_cid_bad_request(monkeypatch: pytest.MonkeyPatch) -> None:
    sess = requests.Session()
    monkeypatch.setattr(sess, "get", lambda *args, **kwargs: DummyResponse(text="", status_code=400))
    assert fetch_pubchem_cid("badname", session=sess) == "unknown"


def test_fetch_pubchem_cid_fallback(monkeypatch: pytest.MonkeyPatch) -> None:
    """The lookup falls back to a broader query when exact matching fails."""
    sess = requests.Session()
    responses = [DummyResponse(text="", status_code=404), DummyResponse(text="789\n")]

    def fake_get(*args, **kwargs):
        return responses.pop(0)

    monkeypatch.setattr(sess, "get", fake_get)
    assert fetch_pubchem_cid("almorexant", session=sess) == "789"


def test_fetch_pubchem_cid_connection_error(monkeypatch: pytest.MonkeyPatch) -> None:
    """Network issues return 'unknown' rather than raising."""
    sess = requests.Session()

    def raise_error(*args, **kwargs):
        raise requests.ConnectionError("boom")

    monkeypatch.setattr(sess, "get", raise_error)
    assert fetch_pubchem_cid("aspirin", session=sess) == "unknown"


# ---------------------------------------------------------------------------
# fetch_pubchem_record tests
# ---------------------------------------------------------------------------

def test_fetch_pubchem_record(monkeypatch: pytest.MonkeyPatch) -> None:
    sess = requests.Session()

    # Stub CID resolution
    monkeypatch.setattr("mylib.pubchem.fetch_pubchem_cid", lambda *a, **k: "2244")

    prop_json = {
        "PropertyTable": {
            "Properties": [
                {
                    "CID": 2244,
                    "CanonicalSMILES": "SMILES",
                    "InChI": "InChI",
                    "InChIKey": "KEY",
                    "MolecularFormula": "C9H8O4",
                    "MolecularWeight": 180.16,
                    "IUPACName": "Name",
                }
            ]
        }
    }
    syn_json = {
        "InformationList": {
            "Information": [
                {"CID": 2244, "Synonym": ["aspirin", "acetylsalicylic acid"]}
            ]
        }
    }

    responses = [
        DummyResponse(json_data=prop_json),
        DummyResponse(json_data=syn_json),
    ]

    def fake_get(url: str, *a, **k):
        return responses.pop(0)

    monkeypatch.setattr(sess, "get", fake_get)

    rec = fetch_pubchem_record("aspirin", session=sess)
    assert rec["pubchem_cid"] == "2244"
    assert rec["canonical_smiles"] == "SMILES"
    assert rec["inchi"] == "InChI"
    assert rec["inchi_key"] == "KEY"
    assert rec["molecular_formula"] == "C9H8O4"
    assert rec["molecular_weight"] == "180.16"
    assert rec["iupac_name"] == "Name"
    assert rec["synonyms"] == "aspirin|acetylsalicylic acid"


def test_fetch_pubchem_record_handles_400(monkeypatch: pytest.MonkeyPatch) -> None:
    """Missing properties or synonyms return empty strings rather than crash."""

    sess = requests.Session()

    # CID resolution succeeds
    monkeypatch.setattr("mylib.pubchem.fetch_pubchem_cid", lambda *a, **k: "42")

    responses = [
        DummyResponse(status_code=400),  # property request returns 400
        DummyResponse(
            json_data={
                "InformationList": {
                    "Information": [
                        {"CID": 42, "Synonym": ["name1", "name2"]}
                    ]
                }
            }
        ),
    ]

    def fake_get(url: str, *a, **k):
        return responses.pop(0)

    monkeypatch.setattr(sess, "get", fake_get)

    rec = fetch_pubchem_record("foo", session=sess)
    assert rec["canonical_smiles"] == ""
    assert rec["synonyms"] == "name1|name2"


def test_fetch_pubchem_record_handles_synonym_400(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    sess = requests.Session()

    monkeypatch.setattr("mylib.pubchem.fetch_pubchem_cid", lambda *a, **k: "42")

    prop_json = {
        "PropertyTable": {
            "Properties": [
                {
                    "CID": 42,
                    "CanonicalSMILES": "S",
                    "InChI": "I",
                    "InChIKey": "K",
                    "MolecularFormula": "F",
                    "MolecularWeight": "W",
                    "IUPACName": "U",
                }
            ]
        }
    }

    responses = [
        DummyResponse(json_data=prop_json),
        DummyResponse(status_code=400),
    ]

    def fake_get(url: str, *a, **k):
        return responses.pop(0)

    monkeypatch.setattr(sess, "get", fake_get)

    rec = fetch_pubchem_record("foo", session=sess)
    assert rec["synonyms"] == ""


def test_fetch_pubchem_record_property_exception(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    sess = requests.Session()

    monkeypatch.setattr("mylib.pubchem.fetch_pubchem_cid", lambda *a, **k: "42")

    syn_json = {
        "InformationList": {
            "Information": [{"CID": 42, "Synonym": ["a", "b"]}]
        }
    }

    responses = [
        requests.ConnectionError("boom"),
        DummyResponse(json_data=syn_json),
    ]

    def fake_get(url: str, *a, **k):
        resp = responses.pop(0)
        if isinstance(resp, Exception):
            raise resp
        return resp

    monkeypatch.setattr(sess, "get", fake_get)

    rec = fetch_pubchem_record("foo", session=sess)
    assert rec["canonical_smiles"] == ""
    assert rec["synonyms"] == "a|b"


def test_fetch_pubchem_record_synonym_exception(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    sess = requests.Session()

    monkeypatch.setattr("mylib.pubchem.fetch_pubchem_cid", lambda *a, **k: "42")

    prop_json = {
        "PropertyTable": {
            "Properties": [
                {
                    "CID": 42,
                    "CanonicalSMILES": "S",
                    "InChI": "I",
                    "InChIKey": "K",
                    "MolecularFormula": "F",
                    "MolecularWeight": "W",
                    "IUPACName": "U",
                }
            ]
        }
    }

    responses = [
        DummyResponse(json_data=prop_json),
        requests.ConnectionError("boom"),
    ]

    def fake_get(url: str, *a, **k):
        resp = responses.pop(0)
        if isinstance(resp, Exception):
            raise resp
        return resp

    monkeypatch.setattr(sess, "get", fake_get)

    rec = fetch_pubchem_record("foo", session=sess)
    assert rec["synonyms"] == ""


# ---------------------------------------------------------------------------
# annotate_pubchem_info tests
# ---------------------------------------------------------------------------

def test_annotate_pubchem_info(monkeypatch: pytest.MonkeyPatch) -> None:
    df = pd.DataFrame({"search_name": ["aspirin", "abcde"]})

    def fake_fetch(name: str, session: requests.Session) -> dict[str, str]:
        base = {
            "pubchem_cid": "111" if name == "aspirin" else "222",
            "canonical_smiles": "S" + name,
            "inchi": "I" + name,
            "inchi_key": "K" + name,
            "molecular_formula": "F" + name,
            "molecular_weight": "W" + name,
            "iupac_name": "U" + name,
            "synonyms": f"{name}|{name}2",
        }
        return base

    monkeypatch.setattr("mylib.pubchem.fetch_pubchem_record", fake_fetch)
    out_df = annotate_pubchem_info(df, session=requests.Session())
    assert out_df["pubchem_cid"].tolist() == ["111", "222"]
    assert out_df["synonyms"].iloc[0] == "aspirin|aspirin2"
