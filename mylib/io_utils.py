"""Utility functions for input/output operations.

This module provides helper routines to load CSV files with automatic
encoding and delimiter detection.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable, Tuple

import pandas as pd


DEFAULT_ENCODINGS: Iterable[str] = (
    "utf-8-sig",
    "utf-16",
    "cp1251",
    "cp1252",
    "latin1",
)
DEFAULT_DELIMITERS: Iterable[str] = (",", ";", "\t", "|")
NA_VALUES = {"", "na", "n/a", "none", "null"}


def smart_read_csv(
    path: str | Path,
    encodings: Iterable[str] = DEFAULT_ENCODINGS,
    seps: Iterable[str] = DEFAULT_DELIMITERS,
) -> Tuple[pd.DataFrame, str, str]:
    """Read a CSV file trying multiple encodings and delimiters.

    Parameters
    ----------
    path:
        Path to the CSV file on disk.
    encodings:
        Candidate text encodings to try. The first successful one is used.
    seps:
        Candidate field delimiters.

    Returns
    -------
    tuple
        DataFrame along with the encoding and delimiter that succeeded.

    Raises
    ------
    UnicodeDecodeError
        If none of the encodings can decode the sample.
    ValueError
        If the file cannot be parsed as CSV with any of the tested combinations.
    """
    path = Path(path)
    sample = path.read_bytes()[:2048]
    last_err: Exception | None = None
    for enc in encodings:
        try:
            snippet = sample.decode(enc)
        except UnicodeDecodeError as err:
            last_err = err
            continue

        # Use csv.Sniffer to guess the delimiter on the decoded sample
        try:
            dialect = csv.Sniffer().sniff(snippet, delimiters="".join(seps))
            sep = dialect.delimiter
        except Exception:
            # Fallback to testing all delimiters if sniffing fails
            for sep in seps:
                try:
                    df = pd.read_csv(  # type: ignore[call-overload]
                        path,
                        dtype=str,
                        encoding=enc,
                        sep=sep,
                        na_values=NA_VALUES,
                        keep_default_na=False,
                    ).fillna("")
                except Exception as err:
                    last_err = err
                    continue
                if (df.iloc[:, 0].astype(str).str.strip() != "").any():
                    return df, enc, sep
            continue

        try:
            df = pd.read_csv(  # type: ignore[call-overload]
                path,
                dtype=str,
                encoding=enc,
                sep=sep,
                na_values=NA_VALUES,
                keep_default_na=False,
            ).fillna("")
        except Exception as err:
            last_err = err
            continue

        if (df.iloc[:, 0].astype(str).str.strip() != "").any():
            return df, enc, sep
    if last_err:
        raise last_err
    raise ValueError(f"Failed to parse CSV file: {path}")


__all__ = ["smart_read_csv"]
