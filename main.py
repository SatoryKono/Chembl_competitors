"""Command-line interface for chemical name normalization."""

from __future__ import annotations

import argparse
import json
import logging
from typing import Any, Dict

import pandas as pd

from mylib import normalize_name, read_input_csv, validate_input, write_output_csv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Path to input CSV with column 'input_name'.")
    parser.add_argument("--output", required=True, help="Destination path for normalized CSV.")
    parser.add_argument("--sep", default=",", help="CSV delimiter (default ',').")
    parser.add_argument("--encoding", default="utf-8", help="File encoding (default utf-8).")
    parser.add_argument("--log-level", default="INFO", help="Logging level (default INFO).")
    return parser.parse_args()


def setup_logging(level: str) -> None:
    logging.basicConfig(level=getattr(logging, level.upper(), logging.INFO), format="%(levelname)s:%(name)s:%(message)s")


def main() -> None:
    args = parse_args()
    setup_logging(args.log_level)

    df = read_input_csv(args.input, sep=args.sep, encoding=args.encoding)
    validate_input(df)

    results: Dict[str, Any] = df["input_name"].apply(normalize_name).apply(pd.Series)
    out_df = pd.concat([df, results], axis=1)

    # Serialize complex columns as JSON strings for CSV compatibility
    out_df["flags"] = out_df["flags"].apply(json.dumps)
    out_df["peptide_info"] = out_df["peptide_info"].apply(json.dumps)
    out_df["oligo_info"] = out_df["oligo_info"].apply(json.dumps)
    out_df["small_molecule_info"] = out_df["small_molecule_info"].apply(json.dumps)

    write_output_csv(out_df, args.output, sep=args.sep, encoding=args.encoding)


if __name__ == "__main__":
    main()
