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
# Bracketed blocks that may contain isotope tokens, e.g. [1-14C]
BRACKET_BLOCK_RE = re.compile(r"\[([^\]]+)\]")
# Generic isotope token matcher used both inside brackets and in free text
ISOTOPE_TOKEN_RE = re.compile(r"(\d{1,3})([A-Za-z]{0,2})", re.IGNORECASE)


def _normalise_element(elem: str) -> str:
    """Normalise element symbols for isotope detection.

    Parameters
    ----------
    elem:
        Raw element symbol extracted from the name.

    Returns
    -------
    str
        Upper-cased element symbol with common typos corrected. ``L`` is
        interpreted as iodine (``I``).
    """

    elem_up = elem.upper()
    return "I" if elem_up == "L" else elem_up


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

    # Strip bracketed blocks containing isotope labels
    def _bracket_repl(match: re.Match[str]) -> str:
        inner = match.group(1)
        tokens = ISOTOPE_TOKEN_RE.findall(inner)
        isotope_found = False
        if len(tokens) == 1 and tokens[0][1] == "":
            # Numeric-only placeholder such as [125]
            mass, _ = tokens[0]
            flags["isotope"].append(f"{mass}I")
            isotope_found = True
        else:
            for mass, elem in tokens:
                elem_up = _normalise_element(elem)
                if elem_up and elem_up in ISOTOPE_ELEMENTS:
                    flags["isotope"].append(f"{mass}{elem_up}")
                    isotope_found = True
        return "" if isotope_found else match.group(0)

    out = BRACKET_BLOCK_RE.sub(_bracket_repl, out)

    # Detect isotopes appearing outside brackets, e.g. 33P or 125I-
    def isotope_repl(match: re.Match[str]) -> str:
        mass, elem = match.groups()
        elem_up = _normalise_element(elem)
        if elem_up in ISOTOPE_ELEMENTS:
            flags["isotope"].append(f"{mass}{elem_up}")
            return ""
        return match.group(0)

    out = ISOTOPE_TOKEN_RE.sub(isotope_repl, out)
    out = re.sub(r"\(\s*\)", "", out)  # drop empty parentheses
    out = re.sub(r"\[\s*\]", "", out)  # drop empty brackets
    out = re.sub(r"-\s*-", "-", out)  # collapse hyphens

    # Collapse extra whitespace and punctuation
    out = re.sub(r"\s+", " ", out)
    out = out.strip().strip(",;-")
    return out, flags


__all__ = ["normalize_name"]
