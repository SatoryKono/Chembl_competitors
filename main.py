"""Command-line interface for PubChem CID lookup by compound name."""

from __future__ import annotations

import argparse
import logging

import pandas as pd

from mylib import annotate_pubchem_cids


# ---------------------------------------------------------------------------
# CLI helpers
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Path to input CSV with column 'search_name'.")
    parser.add_argument("--output", required=True, help="Destination path for CSV with PubChem CIDs.")
    parser.add_argument("--sep", default=",", help="CSV delimiter (default ',').")
    parser.add_argument("--encoding", default="utf-8", help="File encoding (default utf-8).")
    parser.add_argument("--log-level", default="INFO", help="Logging level (default INFO).")
    return parser.parse_args()


def setup_logging(level: str) -> None:
    """Configure basic logging."""

    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(levelname)s:%(name)s:%(message)s",
    )


# ---------------------------------------------------------------------------
# Main execution
# ---------------------------------------------------------------------------

def main() -> None:
    """Entry point for PubChem CID lookup."""

    args = parse_args()
    setup_logging(args.log_level)

    df = pd.read_csv(args.input, sep=args.sep, encoding=args.encoding)
    out_df = annotate_pubchem_cids(df, name_column="search_name")
    out_df.to_csv(args.output, sep=args.sep, encoding=args.encoding, index=False)


if __name__ == "__main__":
    main()
