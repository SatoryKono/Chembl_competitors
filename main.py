"""Command-line interface for name normalisation.

Usage
-----
Run ``python main.py --input examples.csv --output normalized.csv`` to read and
normalise compound names. The resulting DataFrame is written to ``normalized.csv``.
If ``--output`` is omitted, the DataFrame is printed to stdout.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

from mylib.io_utils import smart_read_csv
from mylib.transforms import normalize_name


log = logging.getLogger(__name__)


def configure_logging(level: str) -> None:
    """Configure application logging.

    Parameters
    ----------
    level:
        Logging level name (e.g. ``"INFO"``).
    """

    logging.basicConfig(level=getattr(logging, level.upper(), "INFO"))


def normalise_file(path: Path, output: Path | None = None) -> pd.DataFrame:
    """Read a CSV file, normalise names and optionally write a CSV output.

    Parameters
    ----------
    path:
        Input CSV file containing an ``input_name`` column.
    output:
        Path to write the normalised CSV. If ``None`` the result is returned but
        not written to disk.

    Returns
    -------
    pandas.DataFrame
        DataFrame with additional ``normalized_name`` and ``flags`` columns.
    """

    df, enc, sep = smart_read_csv(path)
    log.info("Loaded %s with encoding %s and delimiter '%s'", path, enc, sep)
    df[["normalized_name", "flags"]] = df["input_name"].apply(
        lambda s: pd.Series(normalize_name(s))
    )
    if output:
        df.to_csv(output, index=False, encoding=enc, sep=sep)
        log.info("Wrote normalised data to %s", output)
    return df


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="CSV file")
    parser.add_argument("--output", type=Path, help="Where to write output CSV")
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    )
    return parser


if __name__ == "__main__":
    args = build_parser().parse_args()
    configure_logging(args.log_level)
    df = normalise_file(args.input, args.output)
    if args.output is None:
        print(df)
