"""Tests for issue-checking helpers."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))

from library.cleanup_competitor_io_validation import check_issues


def test_check_issues() -> None:
    df = pd.DataFrame(
        [
            {"input_name": "oligo ACGTACGT", "category": "small_molecule", "flags": {}, "oligo_info": {}},
            {
                "input_name": "oligo",
                "category": "oligonucleotide",
                "flags": {"oligo_mods": []},
                "oligo_info": {"sequences": [], "mods": {"five_prime": [], "three_prime": [], "internal": []}},
            },
            {
                "input_name": "oligo",
                "category": "oligonucleotide",
                "flags": {"oligo_mods": []},
                "oligo_info": {
                    "sequences": [{"seq": "AC", "length": 2}],
                    "mods": {"five_prime": [], "three_prime": [], "internal": []},
                },
            },
            {
                "input_name": "oligo",
                "category": "oligonucleotide",
                "flags": {"oligo_mods": []},
                "oligo_info": {
                    "sequences": [{"seq": "ACGTACGT", "length": 8}],
                    "mods": {"five_prime": ["Phos"], "three_prime": [], "internal": []},
                },
            },
        ]
    )
    res = check_issues(df)
    assert "oligo_missed" in res.loc[0, "issues"]
    assert "oligo_parse_failed" in res.loc[1, "issues"]
    assert "oligo_len_suspect" in res.loc[2, "issues"]
    assert "oligo_mod_unparsed" in res.loc[3, "issues"]
