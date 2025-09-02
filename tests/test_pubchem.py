"""Tests for PubChem lookup utilities."""

from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd
import pytest
import requests

sys.path.append(str(Path(__file__).resolve().parents[1]))

from mylib import annotate_pubchem_cids
from mylib.pubchem import fetch_pubchem_cid


class DummyResponse:
    """Minimal response stub for requests Session.get."""

    def __init__(self, text: str, status_code: int = 200) -> None:
        self.text = text
        self.status_code = status_code

    def raise_for_status(self) -> None:
        if self.status_code >= 400 and self.status_code != 404:
            raise requests.HTTPError(f"HTTP {self.status_code}")


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
    monkeypatch.setattr(sess, "get", lambda *args, **kwargs: DummyResponse("123\n"))
    assert fetch_pubchem_cid("aspirin", session=sess) == "123"


def test_fetch_pubchem_cid_multiple(monkeypatch: pytest.MonkeyPatch) -> None:
    sess = requests.Session()
    monkeypatch.setattr(sess, "get", lambda *args, **kwargs: DummyResponse("1\n2\n"))
    assert fetch_pubchem_cid("foo bar", session=sess) == "multiply"


def test_fetch_pubchem_cid_unknown(monkeypatch: pytest.MonkeyPatch) -> None:
    sess = requests.Session()
    monkeypatch.setattr(sess, "get", lambda *args, **kwargs: DummyResponse("", 404))
    assert fetch_pubchem_cid("unknowncompound", session=sess) == "unknown"


# ---------------------------------------------------------------------------
# annotate_pubchem_cids tests
# ---------------------------------------------------------------------------

def test_annotate_pubchem_cids(monkeypatch: pytest.MonkeyPatch) -> None:
    df = pd.DataFrame({"search_name": ["aspirin", "abcde"]})

    def fake_fetch(name: str, session: requests.Session) -> str:
        return "111" if name == "aspirin" else "222"

    monkeypatch.setattr("mylib.pubchem.fetch_pubchem_cid", fake_fetch)
    out_df = annotate_pubchem_cids(df, session=requests.Session())
    assert out_df["pubchem_cid"].tolist() == ["111", "222"]
