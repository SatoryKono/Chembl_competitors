"""Command-line interface for name normalisation.

Usage
-----
Run ``python main.py --input examples.csv`` to read and normalise compound
names. The resulting DataFrame is printed to stdout.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from mylib.io_utils import smart_read_csv
from mylib.transforms import normalize_name


log = logging.getLogger(__name__)


def configure_logging(level: str) -> None:
    logging.basicConfig(level=getattr(logging, level.upper(), "INFO"))


def normalise_file(path: Path) -> None:
    df, enc, sep = smart_read_csv(path)
    log.info("Loaded %s with encoding %s and delimiter '%s'", path, enc, sep)
    df[["normalized_name", "flags"]] = df["input_name"].apply(
        lambda s: pd.Series(normalize_name(s))
    )
    print(df)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="CSV file")
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    )
    return parser


if __name__ == "__main__":
    import pandas as pd

    args = build_parser().parse_args()
    configure_logging(args.log_level)
    normalise_file(args.input)
