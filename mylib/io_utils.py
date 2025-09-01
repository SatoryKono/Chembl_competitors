"""Input/output helpers for chemical normalization.

This module provides convenience wrappers around pandas CSV readers and
writers with additional validation.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pandas as pd

logger = logging.getLogger(__name__)


REQUIRED_COLUMNS = ["input_name"]


def read_input_csv(
    path: str | Path,
    *,
    sep: str = ",",
    encoding: str = "utf-8",
) -> pd.DataFrame:
    """Load input CSV and validate required columns.

    Parameters
    ----------
    path:
        Path to the input CSV file.
    sep:
        Field separator used in the CSV.
    encoding:
        Encoding of the CSV file.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing at least the required columns.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If required columns are missing.
    """

    path = Path(path)
    logger.debug("Reading input CSV from %s", path)

    try:
        df = pd.read_csv(path, sep=sep, encoding=encoding)
        # Heuristic: if the index is not a simple RangeIndex or extra columns
        # appear, treat the file as malformed and trigger fallback parsing.
        if not isinstance(df.index, pd.RangeIndex) or df.shape[1] != len(REQUIRED_COLUMNS):
            raise pd.errors.ParserError("irregular column structure")
    except pd.errors.ParserError as exc:
        # If names contain unescaped delimiters, fall back to line-based parsing
        logger.warning("Standard CSV parsing failed: %s", exc)
        with Path(path).open(encoding=encoding) as handle:
            lines = handle.read().splitlines()
        if not lines:
            msg = "Input file is empty"
            logger.error(msg)
            raise ValueError(msg) from exc
        header = lines[0].strip()
        if header.lower() != "input_name":
            msg = "Missing required column 'input_name' in header"
            logger.error(msg)
            raise ValueError(msg) from exc
        df = pd.DataFrame({"input_name": [line.strip() for line in lines[1:]]})


    missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing:
        msg = f"Missing required columns: {missing}"
        logger.error(msg)
        raise ValueError(msg)
    return df


def write_output_csv(
    df: pd.DataFrame,
    path: str | Path,
    *,
    sep: str = ",",
    encoding: str = "utf-8",
    index: bool = False,
    **to_csv_kwargs: Any,
) -> None:
    """Write DataFrame to CSV.

    Parameters
    ----------
    df:
        DataFrame to save.
    path:
        Destination file path.
    sep:
        Field separator.
    encoding:
        Text encoding for the output file.
    index:
        Whether to write DataFrame index.
    to_csv_kwargs:
        Additional keyword arguments forwarded to :meth:`pandas.DataFrame.to_csv`.
    """

    path = Path(path)
    logger.debug("Writing output CSV to %s", path)
    df.to_csv(path, sep=sep, encoding=encoding, index=index, **to_csv_kwargs)
