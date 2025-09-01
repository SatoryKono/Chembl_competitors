"""Validation helpers for chemical name normalization."""

from __future__ import annotations

import logging
from typing import Iterable

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
