"""Text normalization utilities for chemical names."""

from __future__ import annotations

import logging
import re
import unicodedata
from typing import Dict, List, Tuple

logger = logging.getLogger(__name__)

# Regex patterns for various flags
PATTERNS: Dict[str, re.Pattern[str]] = {
    "fluorophore": re.compile(
        r"\b(FITC|Alexa\s?\d+|Cy\d+|Rhodamine|TRITC|DAPI|Texas\sRed|PE|PerCP|APC)\b",
        re.IGNORECASE,
    ),
    "isotope": re.compile(
        r"\b(\[\d+[A-Za-z]\]|\d+[A-Za-z]|[dD]\d+|U-?13C|tritiated|deuterated)\b",
        re.IGNORECASE,
    ),
    "biotin": re.compile(r"\bbiotin(?:ylated)?\b", re.IGNORECASE),
    "salt": re.compile(
        r"\b(hydrochloride|phosphate|mesylate|citrate|tartrate|acetate|sulfate|nitrate|maleate|fumarate|oxalate|sodium|potassium|calcium|lithium)\b",
        re.IGNORECASE,
    ),
    "hydrate": re.compile(r"\b(?:mono|di|tri|tetra|penta)?hydrate|anhydrous\b", re.IGNORECASE),
    "noise": re.compile(
        r"\b(solution|aqueous|buffer|USP|EP|ACS|reagent|analytical|grade|powder|crystalline|purity|lyophilized)\b",
        re.IGNORECASE,
    ),
}

AA1 = set("ACDEFGHIKLMNPQRSTVWY")
AA3 = {
    "Ala",
    "Cys",
    "Asp",
    "Glu",
    "Phe",
    "Gly",
    "His",
    "Ile",
    "Lys",
    "Leu",
    "Met",
    "Asn",
    "Pro",
    "Gln",
    "Arg",
    "Ser",
    "Thr",
    "Val",
    "Trp",
    "Tyr",
}


def _unicode_normalize(text: str) -> str:
    """Apply Unicode normalization and basic whitespace fixes."""

    text = unicodedata.normalize("NFKC", text)
    text = text.replace("µ", "u")
    # Remove NBSP and zero-width spaces
    text = re.sub(r"[\u00A0\u200B\u200C\u200D\uFEFF]", "", text)
    text = re.sub(r"\s+", " ", text).strip()
    return text


def _fix_spacing(text: str) -> str:
    """Remove spaces around punctuation like '-', '/', ':', '+'."""

    return re.sub(r"\s*([-/:+])\s*", r"\1", text)


def _remove_concentrations(text: str, flags: Dict[str, List[str]]) -> str:
    pattern = re.compile(
        r"\b\d+(?:\.\d+)?\s*(?:mM|M|uM|µM|nM|pM|%|mg/mL|mg\/mL|g/mL|mg|g|mL)\b",
        re.IGNORECASE,
    )
    matches = pattern.findall(text)
    if matches:
        flags.setdefault("noise", []).extend(matches)
        text = pattern.sub(" ", text)
    return text


def _remove_parenthetical(text: str, flags: Dict[str, List[str]]) -> str:
    pattern = re.compile(r"(\([^)]*\)|\[[^]]*\])")
    matches = pattern.findall(text)
    if matches:
        flags.setdefault("parenthetical", []).extend(matches)
        text = pattern.sub(" ", text)
    return text


def _detect_and_remove(text: str, key: str, flags: Dict[str, List[str]]) -> str:
    pattern = PATTERNS[key]
    matches = pattern.findall(text)
    if matches:
        flags.setdefault(key, []).extend(matches if isinstance(matches, list) else [matches])
        text = pattern.sub(" ", text)
    return text


def _cleanup(text: str) -> str:
    text = re.sub(r"\s+", " ", text)
    text = text.strip(" -/:,+")
    return text.strip()


def _detect_peptide(text: str) -> Tuple[str, Dict[str, str]]:
    """Detect peptide-like strings and return category and info."""

    lowered = text.lower()
    if re.search(r"poly(?:-|\()[a-z:,]+", lowered):
        return "peptide", {"type": "polymer"}
    if re.search(r"\b(?:peptide|oligopeptide|polypeptide)\b", lowered):
        return "peptide", {"type": "aa_terms"}

    tokens = re.split(r"[-:\s]+", text)
    if len(tokens) >= 2:
        if all(t.upper() in AA1 for t in tokens):
            return "peptide", {"type": "sequence_like"}
        if all(t[:1].upper() + t[1:].lower() in AA3 for t in tokens if t):
            return "peptide", {"type": "sequence_like"}
    return "small_molecule", {}


def normalize_name(name: str) -> Dict[str, object]:
    """Normalize a chemical name and extract flags.

    Parameters
    ----------
    name:
        Raw input chemical name.

    Returns
    -------
    dict
        Dictionary with normalized fields.
    """

    flags: Dict[str, List[str]] = {}
    text = _unicode_normalize(name)
    text = _fix_spacing(text)
    base_clean = text  # for fallback

    text = _remove_concentrations(text, flags)
    for key in ["isotope", "fluorophore", "biotin", "salt", "hydrate"]:
        text = _detect_and_remove(text, key, flags)
    text = _remove_parenthetical(text, flags)
    text = _detect_and_remove(text, "noise", flags)

    text = _cleanup(text)
    category, peptide_info = _detect_peptide(text)

    if not text:
        text = base_clean

    normalized_name = text
    search_name = _cleanup(text).lower()

    result = {
        "normalized_name": normalized_name,
        "search_name": search_name,
        "category": category,
        "peptide_info": peptide_info,
        "flags": flags,
        "flag_isotope": bool(flags.get("isotope")),
        "flag_fluorophore": bool(flags.get("fluorophore")),
        "flag_biotin": bool(flags.get("biotin")),
        "flag_salt": bool(flags.get("salt")),
        "flag_hydrate": bool(flags.get("hydrate")),
    }
    return result
