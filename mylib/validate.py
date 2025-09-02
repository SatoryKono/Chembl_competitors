"""Validation helpers for chemical name normalization."""

from __future__ import annotations

import logging
from typing import Iterable

from .transforms import _has_oligo_signal


import pandas as pd

logger = logging.getLogger(__name__)


def validate_input(df: pd.DataFrame, required: Iterable[str] | None = None) -> None:
    """Validate input DataFrame columns.

    Parameters
    ----------
    df:
        Input DataFrame to validate.
    required:
        Additional required column names.

    Raises
    ------
    ValueError
        If required columns are missing.
    """

    required_cols = set(required or []) | {"input_name"}
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        msg = f"Missing required columns: {missing}"
        logger.error(msg)
        raise ValueError(msg)


def check_issues(df: pd.DataFrame) -> pd.DataFrame:
    """Run basic issue checks on normalized output.

    The function adds an ``issues`` column where each row lists pipe-separated
    problem codes. Currently supported codes are ``oligo_missed``,
    ``oligo_parse_failed``, ``oligo_mod_unparsed`` and ``oligo_len_suspect``.
    """

    issue_list = []
    for _, row in df.iterrows():
        row_issues: list[str] = []
        flags = row.get("flags", {})
        info = row.get("oligo_info", {})
        name = row.get("input_name", "")

        if _has_oligo_signal(str(name)) and row.get("category") != "oligonucleotide":
            row_issues.append("oligo_missed")

        if row.get("category") == "oligonucleotide":
            seqs = info.get("sequences", [])
            if not seqs:
                row_issues.append("oligo_parse_failed")
            mods_recorded = set(flags.get("oligo_mods", []))
            all_mods = set(
                info.get("mods", {}).get("five_prime", [])
                + info.get("mods", {}).get("three_prime", [])
                + info.get("mods", {}).get("internal", [])
            )
            if all_mods and not all_mods.issubset(mods_recorded):
                row_issues.append("oligo_mod_unparsed")
            total_len = sum(seq.get("length", 0) for seq in seqs)
            if total_len and (total_len < 8 or total_len > 200):
                row_issues.append("oligo_len_suspect")

        issue_list.append("|".join(row_issues))

    result = df.copy()
    result["issues"] = issue_list
    return result
