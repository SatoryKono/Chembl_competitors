"""Name normalisation utilities for chemical compounds."""

from __future__ import annotations

import re
from typing import Dict, List, Tuple

# Known isotopic element symbols. Using uppercase for comparison
ISOTOPE_ELEMENTS = {
    "H",
    "C",
    "N",
    "O",
    "P",
    "S",
    "I",
    "F",
    "CL",
    "BR",
}

# Regular expression for stereochemistry descriptors like (4S,5R)
RS_PATTERN = re.compile(r"\((?:\d+[RSrs](?:[,/-]\d+[RSrs])*)\)")
# Regular expression for isotope notations such as [3H], (14C) or 125I-
ISOTOPE_PATTERN = re.compile(r"\[?(\d{1,3})([A-Za-z]{1,2})\]?", re.IGNORECASE)


def normalize_name(name: str) -> Tuple[str, Dict[str, List[str]]]:
    """Normalise a compound name and capture annotation flags.

    Parameters
    ----------
    name:
        Raw compound name.

    Returns
    -------
    tuple
        Normalised name and dictionary with detected flags. The dictionary
        contains keys ``isotope``, ``fluorophore`` and ``noise`` amongst
        others, each mapping to a list of strings found in the input.
    """
    flags: Dict[str, List[str]] = {
        "fluorophore": [],
        "isotope": [],
        "biotin": [],
        "salt": [],
        "hydrate": [],
        "noise": [],
        "parenthetical": [],
    }

    if not name:
        return "", flags

    out = name.strip()

    # Remove common "labeled" noise but keep track
    if "labeled" in out.lower():
        flags["noise"].append("labeled")
        out = re.sub(r"\b[- ]*labeled\b", "", out, flags=re.IGNORECASE)

    # Detect and strip fluorophore tags such as FITC or FAM
    for fl in ("fitc", "fam"):
        if fl in out.lower():
            flags["fluorophore"].append(fl)
            out = re.sub(fl, "", out, flags=re.IGNORECASE)

    # Remove stereochemistry descriptors before isotope detection
    rs_matches = RS_PATTERN.findall(out)
    if rs_matches:
        flags["parenthetical"].extend(rs_matches)
        out = RS_PATTERN.sub("", out)

    # Detect isotopes
    def isotope_repl(match: re.Match[str]) -> str:
        mass, elem = match.groups()
        elem_up = elem.upper()
        if elem_up in ISOTOPE_ELEMENTS:
            flags["isotope"].append(f"{mass}{elem_up}")
            return ""
        return match.group(0)

    out = ISOTOPE_PATTERN.sub(isotope_repl, out)
    out = re.sub(r"\(\s*\)", "", out)  # drop empty parentheses left by isotope removal

    # Collapse extra whitespace and punctuation
    out = re.sub(r"\s+", " ", out)
    out = out.strip().strip(",;")
    return out, flags


__all__ = ["normalize_name"]
