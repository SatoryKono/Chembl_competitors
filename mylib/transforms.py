"""Text normalization utilities for chemical names."""

from __future__ import annotations

import logging
import re
import unicodedata
from typing import Dict, List, Tuple

logger = logging.getLogger(__name__)

# Tokens representing salts and mineral acids to strip early in processing
SALT_TOKENS = [
    "hydrochloride",
    "phosphate",
    "mesylate",
    "citrate",
    "tartrate",
    "acetate",
    "sulfate",
    "nitrate",
    "maleate",
    "fumarate",
    "oxalate",
    "sodium",
    "potassium",
    "calcium",
    "lithium",
    "HCl",
    "HBr",
    "HNO3",
    "H2SO4",
]

# Tokens representing hydrate forms and related descriptors
HYDRATE_TOKENS = [
    "monohydrate",
    "dihydrate",
    "trihydrate",
    "tetrahydrate",
    "pentahydrate",
    "hydrate",
    "anhydrous",
]

# Regex patterns for various flags
PATTERNS: Dict[str, re.Pattern[str]] = {
    "fluorophore": re.compile(
        r"""
        \b(
            FITC|
            Alexa(?:\s|-)?Fluor[\s-]?\d+|
            HiLyte(?:\s|-)?(?:Fluor)?[\s-]?\d+|
            DyLight[\s-]?\d+|
            CF[\s-]?\d+|
            Janelia(?:\s|-)?Fluor[\s-]?\d+|
            BODIPY(?:[-/][A-Za-z0-9/]+|\s(?:[A-Za-z]{1,3}|[0-9/]+)){0,2}|
            Cy\d+|
            Rhodamine|
            TRITC|
            DAPI|
            Texas\sRed|
            PE|
            PerCP|
            APC
        )\b
        """,
        re.IGNORECASE | re.VERBOSE,
    ),
    "isotope": re.compile(
        r"(?<!\w)(?:"
        r"\[(?:3H|14C|13C|15N|2H|125I|18F)\]"  # bracketed isotopes
        r"|(?:3H|14C|13C|15N|2H|125I|18F|D|T)"    # bare prefixes and single letters
        r"|d\d+"                                  # deuteration like d5
        r"|U-?13C"                                 # uniform 13C labeling
        r"|tritiated|deuterated"                   # descriptive words
        r")(?!\w)",
        re.IGNORECASE,
    ),
    "biotin": re.compile(r"\bbiotin(?:ylated)?\b", re.IGNORECASE),
    "salt": re.compile(
        r"\b(" + "|".join(map(re.escape, SALT_TOKENS)) + r")\b",
        re.IGNORECASE,
    ),
    "hydrate": re.compile(
        r"\b(" + "|".join(map(re.escape, HYDRATE_TOKENS)) + r")\b",
        re.IGNORECASE,
    ),
}

# Non-structural descriptor tokens to remove in two passes
NOISE_WORDS = [
    "solution",
    "soln",
    "aqueous",
    r"aq\.",
    "stock",
    "buffer",
    "USP",
    "EP",
    "ACS",
    "reagent",
    "analytical",
    "grade",
    "crystal(?:line)?",
    "powder",
    "PBS",
]
NOISE_WITH_ID = r"(?:lot|cat(?:alog)?|code|ref)[\s:_-]*\w+"
PURITY_PATTERN = r"≥?\s*\d{1,2}\s*%?\s*purity"
NOISE_REGEX = re.compile(
    rf"{NOISE_WITH_ID}|\b(?:{'|'.join(NOISE_WORDS)})\b|{PURITY_PATTERN}",
    re.IGNORECASE,
)

STOPWORDS = {"in", "of", "and"}

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


def _flatten_flags(flags: Dict[str, List[str]]) -> str:
    """Flatten selected flag tokens into a pipe-delimited string."""

    order = [
        "fluorophore",
        "isotope",
        "biotin",
        "salt",
        "hydrate",
        "noise",
        "parenthetical",
    ]
    parts: List[str] = []
    for key in order:
        for token in flags.get(key, []):
            parts.append(f"{key}:{token}")
    return "|".join(parts)


def _unicode_normalize(text: str) -> str:
    """Apply Unicode normalization and basic whitespace fixes."""

    text = unicodedata.normalize("NFKC", text)
    text = text.replace("µ", "u")
    # Remove NBSP and zero-width spaces
    text = re.sub(r"[\u00A0\u200B\u200C\u200D\uFEFF]", "", text)
    text = re.sub(r"\s+", " ", text).strip()
    return text


def _fix_spacing(text: str) -> str:
    """Normalize spacing around punctuation and decimals.

    In addition to compacting spaces around ``-``, ``/``, ``:``, and ``+``,
    this function also removes extraneous spaces surrounding commas and
    periods. Decimal numbers such as ``1 . 5`` are collapsed to ``1.5``.
    """

    # Compact common connector characters
    text = re.sub(r"\s*([-/:+,])\s*", r"\1", text)
    # Remove spaces around periods, including decimal numbers
    text = re.sub(r"(?<=\d)\s*\.\s*(?=\d)", ".", text)
    text = re.sub(r"\s*\.\s*", ".", text)
    return text





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


def _detect_and_remove(text: str, key: str, flags: Dict[str, List[str]]) -> str:
    pattern = PATTERNS[key]
    matches = pattern.findall(text)
    if matches:
        flags.setdefault(key, []).extend(matches if isinstance(matches, list) else [matches])
        text = pattern.sub(" ", text)
    return text


def _remove_noise_descriptors(text: str, flags: Dict[str, List[str]]) -> str:
    """Remove non-structural descriptors inside and outside brackets."""

    def _log_noise(match: re.Match[str]) -> str:
        token = match.group(0)
        flags.setdefault("noise", []).append(token)
        return ""

    def _strip_bracket(match: re.Match[str]) -> str:
        inner = match.group(0)[1:-1]
        cleaned = NOISE_REGEX.sub(_log_noise, inner)
        cleaned = cleaned.strip(" ,;:-")
        if cleaned.lower() in STOPWORDS:
            cleaned = ""
        if cleaned:
            return cleaned
        flags.setdefault("parenthetical", []).append(match.group(0))
        return ""

    text = re.sub(r"\([^()]*\)|\[[^\[\]]*\]", _strip_bracket, text)
    text = NOISE_REGEX.sub(_log_noise, text)
    return text


def _cleanup(text: str) -> str:
    """Final whitespace and punctuation cleanup."""

    # Remove errant spaces around connectors that may appear after token removal
    text = _fix_spacing(text)
    # Collapse multiple whitespace characters to single spaces
    text = re.sub(r"\s+", " ", text)
    # Re-run spacing fix in case the previous collapse introduced new gaps
    text = _fix_spacing(text)
    # Consolidate repeated connectors that may result from removals
    text = re.sub(r"([-/:+]){2,}", r"\1", text)
    # Drop leading/trailing punctuation and whitespace
    text = text.strip(" -/:,+")
    return text.strip()


def _detect_peptide(text: str) -> Tuple[str, Dict[str, str]]:
    """Detect peptide-like strings and return category and info."""

    lowered = text.lower()


    # polymer-style notation: poly-Glu:Tyr, poly (Glu, Tyr), poly Glu Tyr
    poly_match = re.search(
        r"\bpoly\b(?:\s*\(\s*|\s+|-)([A-Za-z]{1,3}(?:[,:\s-]+[A-Za-z]{1,3})*)\)?",
        text,
    )
    if poly_match:
        comp_tokens = [t for t in re.split(r"[,:\s-]+", poly_match.group(1)) if t]
        comp_clean = [t.lower() for t in comp_tokens]
        if comp_clean and all(
            t.upper() in AA1 or t[:1].upper() + t[1:].lower() in AA3
            for t in comp_clean
        ):
            return "peptide", {
                "type": "polymer",
                "composition": ":".join(comp_clean),
            }

    if re.search(r"\b(?:peptide|oligopeptide|polypeptide)\b", lowered):
        return "peptide", {"type": "aa_terms"}

    tokens = [t for t in re.split(r"[-:,\s]+", text) if t]
    protect = {"H", "AC", "BOC", "OH", "NH2"}
    tokens_clean = [t for t in tokens if t.upper() not in protect]
    if len(tokens_clean) >= 2:
        if all(t.upper() in AA1 for t in tokens_clean):
            return "peptide", {"type": "sequence_like"}
        if all(
            t[:1].upper() + t[1:].lower() in AA3 for t in tokens_clean if t
        ):

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
        Dictionary with normalized fields. By default ``search_name`` equals
        ``normalized_name``; if they differ an explanatory string is stored in
        ``search_override_reason``.
    """

    flags: Dict[str, List[str]] = {}
    text = _unicode_normalize(name)
    text = _fix_spacing(text)
    base_clean = text  # for fallback

    # Strip fluorophore labels before any other processing to avoid misclassifying
    # numeric components as concentrations or noise.
    text = _detect_and_remove(text, "fluorophore", flags)
    text = _remove_concentrations(text, flags)
    # Remove remaining flagged tokens
    for key in ["salt", "isotope", "biotin", "hydrate"]:
        text = _detect_and_remove(text, key, flags)
    text = _remove_noise_descriptors(text, flags)

    text = _cleanup(text)
    category, peptide_info = _detect_peptide(text)

    status = ""
    flag_empty_after_clean = False
    if not text:
        # Fall back to the base-clean string and ensure it is fully cleaned
        text = _cleanup(base_clean)
        if not text:


            # As a last resort, minimally normalize the original text
            text = _cleanup(_unicode_normalize(name))



        status = "empty_after_clean"
        flag_empty_after_clean = True
        logger.warning("Name empty after cleaning; using fallback: %r", name)




    # Ensure spacing is compact after any late fallbacks
    normalized_name = _fix_spacing(text).lower()


    search_name = normalized_name
    removed_tokens_flat = _flatten_flags(flags)

    result = {
        "normalized_name": normalized_name,
        "search_name": search_name,
        "search_override_reason": "",
        "category": category,
        "peptide_info": peptide_info,
        "flags": flags,
        "removed_tokens_flat": removed_tokens_flat,
        "status": status,
        "flag_isotope": bool(flags.get("isotope")),
        "flag_fluorophore": bool(flags.get("fluorophore")),
        "flag_biotin": bool(flags.get("biotin")),
        "flag_salt": bool(flags.get("salt")),
        "flag_hydrate": bool(flags.get("hydrate")),
        "flag_empty_after_clean": flag_empty_after_clean,
    }
    return result
