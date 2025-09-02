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
            FAM|
            Alexa(?:\s|-)?Fluor[\s-]?\d+|
            HiLyte(?:\s|-)?(?:Fluor)?[\s-]?\d+|
            DyLight[\s-]?\d+|
            CF[\s-]?\d+|
            Janelia(?:\s|-)?Fluor[\s-]?\d+|
            BODIPY(?:[-/][A-Za-z0-9/]+|\s(?:[A-Za-z]{1,3}|[0-9/]+)){0,2}|
            Cy\d+|
            Rhodamine|
            AMC|AFC|ACC|pNP|
            BHQ\d*|
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
        r"""
        (
            \[\s*\d{1,3}\s*[A-Z][a-z]?\s*\]                           # bracketed forms
            |(?<![A-Za-z0-9])(?:3H|2H|D|T|13C|14C|15N|18F|32P|86Rb|125I)(?![A-Za-z0-9]) # bare prefixes
            |\bd\d+\b                                                    # d-number deuteration
            |\b(?:deuterated|tritiated|U-?13C)\b                          # words
        )
        """,
        re.IGNORECASE | re.VERBOSE,
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

# ---------------------------------------------------------------------------
# Cyclic nucleotide detection
# ---------------------------------------------------------------------------

# Patterns covering common representations of cyclic nucleotides. These
# capture variants such as ``3',5'-cAMP``, ``c-di-gmp`` and ``cGAMP`` while
# allowing for stray spaces or mixed delimiters.
CNUC_PATTERNS = [
    re.compile(r"(?i)\b(?:3'?[, ]?\s*5'?)\s*-\s*c\s*[acgut]\s*(?:mp|dp|tp)\b"),
    re.compile(r"(?i)\bc\s*[acgut]\s*(?:mp|dp|tp)\b"),
    re.compile(r"(?i)\bc-?di-?\s*(?:amp|gmp|cmp|ump|gamp)\b"),
    re.compile(r"(?i)\b(?:(?:2'[, ]?3'|3'[, ]?5')\s*-\s*)?c\s*gamp\b"),
]


def is_cyclic_nucleotide(text: str) -> bool:
    """Return ``True`` if ``text`` appears to denote a cyclic nucleotide."""

    return any(p.search(text) for p in CNUC_PATTERNS)


def canonicalize_cyclic_nucleotide(text: str) -> str:
    """Canonicalize cyclic nucleotide names.

    Parameters
    ----------
    text:
        Raw text that already passed through spacing normalization.

    Returns
    -------
    str
        Canonical form such as ``3',5'-cAMP`` or ``c-di-GMP``.
    """

    # 3',5'-cNMP with optional spaces/commas
    text = re.sub(
        r"(?i)3'[, ]?\s*5'-\s*c\s*([acgut])\s*(mp|dp|tp)",
        lambda m: f"3',5'-c{m.group(1).upper()}{m.group(2).upper()}",
        text,
    )
    # 5'-cNMP prefixes
    text = re.sub(
        r"(?i)5'-\s*c\s*([acgut])\s*(mp|dp|tp)",
        lambda m: f"5'-c{m.group(1).upper()}{m.group(2).upper()}",
        text,
    )
    # Simple cNMP forms
    text = re.sub(
        r"(?i)\bc\s*([acgut])\s*(mp|dp|tp)\b",
        lambda m: f"c{m.group(1).upper()}{m.group(2).upper()}",
        text,
    )
    # c-di-XXX forms
    text = re.sub(
        r"(?i)\bc-?di-?\s*(amp|gmp|cmp|ump|gamp)\b",
        lambda m: f"c-di-{m.group(1).upper()}",
        text,
    )
    # cGAMP with optional prime orientation prefixes
    text = re.sub(
        r"(?i)(?:2'[, ]?3'|3'[, ]?5')?-?c\s*gamp",
        "cGAMP",
        text,
    )
    return _fix_spacing(text)

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

# Keywords and patterns for oligonucleotide detection
OLIGO_KEYWORDS = [
    "oligo",
    "oligonucleotide",
    "primer",
    "probe",
    "aptamer",
    "sirna",
    "shrna",
    "mirna",
    "antisense",
    "aso",
    "morpholino",
    "lna",
    "grna",
    "sgrna",
    "crrna",
    "tracrrna",
    "ribo",
    "g-block",
    "gene fragment",
    "dna",
    "rna",
    "guide",
    "target",
    "protospacer",
    "pam",
]

NUCLEO_PATTERN = re.compile(r"\b[ACGTURYKMSWBDHVN]{8,}\b", re.IGNORECASE)
ROLE_PATTERN = re.compile(
    r"(sense|antisense|guide|tracrrna|crrna)[:\s]+([-ACGTURYKMSWBDHVN\s]+)",
    re.IGNORECASE,
)
SLASH_MOD_PATTERN = re.compile(r"/([^/]+)/")


def _valid_nuc_sequence(seq: str) -> bool:
    """Check if ``seq`` appears to be a nucleotide sequence."""

    seq = seq.upper()
    if re.search(r"[^ACGTURYKMSWBDHVN]", seq):
        return False
    if len(seq) < 8:
        return False
    bases = sum(1 for c in seq if c in "ACGTU")
    degens = sum(1 for c in seq if c in "RYSWKMBDHVN")
    if bases / len(seq) < 0.6:
        return False
    if degens / len(seq) > 0.4:
        return False
    if re.search(r"[EFILMPQZ]", seq):
        return False
    return True


def _has_oligo_signal(text: str) -> bool:
    """Return True if text contains oligonucleotide hints."""

    lower = text.lower()
    if re.search(r"\b(" + "|".join(map(re.escape, OLIGO_KEYWORDS)) + r")\b", lower):
        return True
    if SLASH_MOD_PATTERN.search(text) or re.search(r"5['\u2032]|3['\u2032]", text):
        return True
    for seq in NUCLEO_PATTERN.findall(text):
        if _valid_nuc_sequence(seq):
            return True
    return False


def _flatten_oligo_flags(flags: Dict[str, object]) -> str:
    """Flatten oligo-related flags into a pipe-delimited string."""

    parts: List[str] = []
    if flags.get("oligo_type"):
        parts.append(f"oligo_type:{flags['oligo_type']}")
    for token in flags.get("oligo_mods", []) or []:
        parts.append(f"oligo_mod:{token}")
    for role in flags.get("oligo_roles", []) or []:
        parts.append(f"oligo_role:{role}")
    if "oligo_len_total" in flags:
        parts.append(f"oligo_len_total:{flags['oligo_len_total']}")
    return "|".join(parts)


def _flatten_flags(flags: Dict[str, List[str]]) -> str:
    """Flatten selected flag tokens into a pipe-delimited string."""

    order = [
        "fluorophore",
        "isotope",
        "biotin",
        "salt",
        "hydrate",
        "chromophore",
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
    """Canonicalize spacing around common punctuation.

    The transformations are applied sequentially to mirror the specification
    in the user instructions. Only the listed punctuation characters are
    affected; other in-line punctuation remains untouched.
    """

    # 0) Normalize various dash characters to a simple hyphen
    text = re.sub(r"[\u2013\u2014\u2212]", "-", text)

    # 1) Remove spaces before closing brackets and after opening brackets
    text = re.sub(r"\s+([)\]\}])", r"\1", text)
    text = re.sub(r"([(\[\{])\s+", r"\1", text)

    # 2) Tighten connector punctuation
    text = re.sub(r"\s*-\s*", "-", text)
    text = re.sub(r"\s*/\s*", "/", text)
    text = re.sub(r"\s*:\s*", ":", text)
    text = re.sub(r"\s*\+\s*", "+", text)
    text = re.sub(r"\s*;\s*", "; ", text)

    # 3) Primes/apostrophes cling to adjacent tokens
    text = re.sub(r"\s*(['\u2032])\s*", r"\1", text)

    # 4) Commas: numeric enumerations keep tight format, lists get a space
    text = re.sub(r"(?<=\d)\s*,\s*(?=\d)", ",", text)
    text = re.sub(r"(?<!\d)\s*,\s*(?!\d)", ", ", text)
    text = re.sub(r",\s+([)\]\}])", r",\1", text)

    # 5) Remove space before trailing hyphens
    text = re.sub(r"\s+-\b", "-", text)

    # 6) Collapse repeated spaces and trim
    text = re.sub(r"\s{2,}", " ", text).strip()

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


def _canonicalize_and_strip_isotopes(text: str, flags: Dict[str, List[str]]) -> str:
    """Canonicalize isotope markers and remove them without artifacts.

    The routine first normalizes common isotopic forms—particularly for
    I-125—and then strips both bracketed and bare tokens while logging all
    matches under ``flags['isotope']``.
    """

    # Normalize dash variants to hyphen for consistent matching
    text = re.sub(r"[\u2013\u2014\u2212]", "-", text)

    # Canonicalize various I-125 notations
    text = re.sub(r"(?i)\[\s*i\s*-\s*?125\s*\]", "[125I]", text)
    text = re.sub(r"(?i)\[\s*125\s*-\s*?i\s*\]", "[125I]", text)
    text = re.sub(r"(?i)\[\s*i\s*125\s*\]", "[125I]", text)
    text = re.sub(r"(?i)\[\s*125\s*i\s*\]", "[125I]", text)
    text = re.sub(r"(?i)\bi\s*-\s*?125\b", "125I", text)
    text = re.sub(r"(?i)\b125\s*-\s*?i\b", "125I", text)
    # cautious replacement of "1251" with "125I" in bracketed or iodo context
    text = re.sub(r"(?i)(?<=\[)\s*125\s*1\s*(?=\])", "125I", text)
    text = re.sub(r"(?i)(\biodo[\w-]{0,20})\[\s*125\s*1\s*\]", r"\1[125I]", text)
    text = re.sub(r"(?i)\[\s*125\s*1\s*\](?=[-\s]*iodo)", "[125I]", text)

    # Canonicalize other bracketed isotopes like [3H], [14C], [32P], [86Rb]
    for iso in ["3H", "14C", "32P", "86Rb"]:
        a, el = iso[:-1], iso[-1]
        text = re.sub(rf"(?i)\[\s*{a}\s*{el}\s*\]", f"[{iso}]", text)

    before = text
    matches = PATTERNS["isotope"].findall(before)
    # Detect prefixes like "[125Ityr0]" which aren't caught by the pattern
    if re.search(r"(?i)\[\s*125\s*I", before) and "[125I]" not in matches:
        matches.append("[125I]")

    if matches:
        norm_matches: List[str] = []
        for m in matches:
            token = m
            lower = token.lower()
            if re.match(r"d\d+", lower) or lower in {"deuterated", "tritiated"}:
                token = lower
            elif lower.startswith("u13c") or lower.startswith("u-13c"):
                token = "U-13C"
            else:
                token = token.upper()
            norm_matches.append(token)
        flags.setdefault("isotope", []).extend(norm_matches)

    # Remove bracketed isotopes while preserving necessary hyphens
    iso_brkt = r"\[\s*(?:125I|3H|14C|32P|86Rb)\s*\]"
    text = re.sub(rf"(?i)-\s*{iso_brkt}\s*-", "-", text)
    text = re.sub(rf"(?i)-\s*{iso_brkt}", "-", text)
    text = re.sub(rf"(?i){iso_brkt}\s*-", "-", text)
    text = re.sub(rf"(?i){iso_brkt}", "", text)
    # Remove prefix 125I inside mixed brackets like [125Ityr0]
    text = re.sub(r"(?i)(?<=\[)\s*125\s*I\s*(?=[A-Za-z0-9])", "", text)

    # Remove remaining isotopic tokens
    text = PATTERNS["isotope"].sub(" ", text)

    # Clean up bracket artifacts
    text = re.sub(r"\[\s*\]", " ", text)
    text = re.sub(r"\[\[", "[", text)
    text = re.sub(r"\]\]", "]", text)

    return text


def _detect_and_remove(text: str, key: str, flags: Dict[str, List[str]]) -> str:
    pattern = PATTERNS[key]
    matches = pattern.findall(text)
    if matches:
        if not isinstance(matches, list):
            matches = [matches]
        norm_matches: List[str] = []
        for m in matches:
            token = m
            if key == "isotope":
                lower = token.lower()
                if re.match(r"d\d+", lower) or lower in {"deuterated", "tritiated"}:
                    token = lower
                else:
                    token = token.upper()
            elif key == "fluorophore" and token.isalpha():
                token = token.upper()
            norm_matches.append(token)
        flags.setdefault(key, []).extend(norm_matches)
        text = pattern.sub(" ", text)
    return text


def _detect_chromophore(text: str, flags: Dict[str, List[str]]) -> str:
    """Detect pNA chromophore tokens without mistaking uppercase PNA."""

    def repl(match: re.Match[str]) -> str:
        token = match.group(0)
        if token.isupper():
            return token
        flags.setdefault("chromophore", []).append(token)
        return " "

    return re.sub(r"\bpna\b", repl, text, flags=re.IGNORECASE)


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


def parse_oligo_segments(text: str, flags: Dict[str, List[str]]) -> Tuple[str, Dict[str, object], Dict[str, object]]:
    """Extract oligonucleotide sequences and modifications.

    Parameters
    ----------
    text:
        Raw text after Unicode and spacing normalization.
    flags:
        Global flags dictionary to populate with fluorophores/biotin markers
        found within modification tags.

    Returns
    -------
    cleaned_text, info, extra_flags
        ``cleaned_text`` has modification tokens removed so subsequent stages
        can operate on residual words. ``info`` captures parsed sequences and
        modification metadata. ``extra_flags`` carries oligo-specific flag
        entries to merge into the main ``flags`` mapping.
    """

    mods = {"five_prime": [], "three_prime": [], "internal": [], "backbone": "UNKNOWN"}
    sequences: List[Dict[str, object]] = []
    roles: List[str] = []
    orig_lower = text.lower()

    # Vendor-style modification tags /token/
    for token in SLASH_MOD_PATTERN.findall(text):
        lower = token.lower()
        if lower.startswith("5"):
            mod = token[1:]
            mods["five_prime"].append(mod)
            if "bio" in mod.lower():
                flags.setdefault("biotin", []).append(mod)
        elif lower.startswith("3"):
            mod = token[1:]
            mods["three_prime"].append(mod)
            if "bio" in mod.lower():
                flags.setdefault("biotin", []).append(mod)
        else:
            mods["internal"].append(token)
        if PATTERNS["fluorophore"].search(token):
            flags.setdefault("fluorophore", []).append(token)
        text = text.replace(f"/{token}/", " ")

    # Five-prime dash-style modifications e.g., 5'-FAM-
    def _five_dash(match: re.Match[str]) -> str:
        token = match.group(1)
        mods["five_prime"].append(token)
        if PATTERNS["fluorophore"].search(token):
            flags.setdefault("fluorophore", []).append(token)
        if token.lower().startswith("bio"):
            flags.setdefault("biotin", []).append(token)
        return ""

    text = re.sub(r"5['\u2032]?-(\w{1,6})-", _five_dash, text)

    # Three-prime dash-style modifications e.g., -BHQ1-3'
    def _three_dash(match: re.Match[str]) -> str:
        token = match.group(1)
        if token:
            mods["three_prime"].append(token)
            if PATTERNS["fluorophore"].search(token):
                flags.setdefault("fluorophore", []).append(token)
            if token.lower().startswith("bio"):
                flags.setdefault("biotin", []).append(token)
        return ""

    text = re.sub(r"-(\w{1,6})-3['\u2032]?", _three_dash, text)
    text = re.sub(r"5['\u2032]?|3['\u2032]?", "", text)

    # Backbone PS indicated by '*' or PS tokens
    if "*" in text or re.search(r"\bps\b", text, re.IGNORECASE):
        mods["backbone"] = "PS"
    else:
        mods["backbone"] = "PO"

    # Role-annotated sequences
    for match in ROLE_PATTERN.finditer(text):
        role = match.group(1).lower()
        seq_raw = match.group(2)
        cleaned = re.sub(r"[^ACGTURYKMSWBDHVN]", "", seq_raw)
        if _valid_nuc_sequence(cleaned):
            seq = cleaned.upper()
            sequences.append({"role": role, "seq": seq, "length": len(seq)})
            roles.append(role)
        text = text.replace(match.group(0), " ")

    # Remaining sequences without explicit roles
    for seq_match in re.findall(r"[ACGTURYKMSWBDHVN*-]{8,}", text, re.IGNORECASE):
        cleaned = re.sub(r"[^ACGTURYKMSWBDHVN]", "", seq_match)
        if _valid_nuc_sequence(cleaned):
            seq = cleaned.upper()
            role = "sense" if not roles else f"seq{len(roles)+1}"
            sequences.append({"role": role, "seq": seq, "length": len(seq)})
            roles.append(role)
        text = text.replace(seq_match, " ")

    total_len = sum(s["length"] for s in sequences)
    seq_concat = "".join(s["seq"] for s in sequences)
    has_u = "U" in seq_concat
    has_t = "T" in seq_concat
    base_type = "RNA" if has_u and not has_t else "DNA" if has_t and not has_u else "UNKNOWN"

    oligo_type = base_type
    if re.search(r"sirna", orig_lower) or ("sense" in roles and "antisense" in roles):
        oligo_type = "siRNA"
    elif re.search(r"grna|sgrna|crrna|tracrrna", orig_lower):
        oligo_type = "CRISPR"
    elif re.search(r"aso|antisense", orig_lower) or mods["backbone"] == "PS":
        oligo_type = "ASO"

    subtype = "NONE"
    if re.search(r"aptamer", orig_lower):
        subtype = "aptamer"
    elif re.search(r"primer", orig_lower):
        subtype = "primer"
    elif re.search(r"probe", orig_lower):
        subtype = "probe"

    info = {
        "type": oligo_type,
        "subtype": subtype,
        "sequences": sequences,
        "mods": mods,
    }

    extra_flags: Dict[str, object] = {
        "oligo": 1,
        "oligo_type": oligo_type,
        "oligo_mods": mods["five_prime"] + mods["three_prime"] + mods["internal"],
        "oligo_roles": roles,
        "oligo_len_total": total_len,
    }

    return text, info, extra_flags


def _cleanup(text: str) -> str:
    """Final whitespace and punctuation cleanup."""

    # Remove errant spaces around connectors that may appear after token removal
    text = _fix_spacing(text)
    # Fix decimals such as ``1 . 5`` -> ``1.5``
    text = re.sub(r"(?<=\d)\s*\.\s*(?=\d)", ".", text)
    text = re.sub(r"\s*\.\s*", ".", text)
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

    if re.search(r"\b(?:peptide|oligopeptide|polypeptide|substrate|histone)\b", lowered) or re.search(
        r"from\s+p\d+", lowered
    ):
        return "peptide", {"type": "aa_terms"}

    tokens = [t for t in re.split(r"[-:,\s]+", text) if t]
    protect = {
        "H",
        "AC",
        "BOC",
        "OH",
        "NH2",
        "PNA",
        "AMC",
        "AFC",
        "ACC",
        "PNP",
        "MEO",
        "SUC",
        "PYROGLU",
    }
    tokens_clean = [t for t in tokens if t.upper() not in protect and t.isalpha()]
    for tok in tokens_clean:
        up = tok.upper()
        if tok.lower() in OLIGO_KEYWORDS:
            continue
        if len(up) >= 6 and all(c in AA1 for c in up) and not _valid_nuc_sequence(up):
            return "peptide", {"type": "sequence_like"}
    if len(tokens_clean) >= 2:
        if all(t.upper() in AA1 for t in tokens_clean):
            return "peptide", {"type": "sequence_like"}
        if all(t[:1].upper() + t[1:].lower() in AA3 for t in tokens_clean if t):
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
    base_clean = text

    # Early removal of chromophore-like tags
    text = _detect_chromophore(text, flags)

    # Remove concentrations and other flagged tokens
    text = _remove_concentrations(text, flags)
    text = _canonicalize_and_strip_isotopes(text, flags)
    text = re.sub(r"(?i)(\d+)(HCl|HBr|HNO3|H2SO4)", r"\1 \2", text)
    for key in ["salt", "biotin", "hydrate"]:
        text = _detect_and_remove(text, key, flags)
    text = _remove_noise_descriptors(text, flags)
    text = _fix_spacing(text)

    small_molecule_info: Dict[str, object] = {}
    peptide_info: Dict[str, str] = {}
    oligo_info: Dict[str, object] = {}

    tmp_no_fluor = _detect_and_remove(text, "fluorophore", {})
    if is_cyclic_nucleotide(tmp_no_fluor):
        text = _detect_and_remove(text, "fluorophore", flags)
        text = canonicalize_cyclic_nucleotide(tmp_no_fluor)
        small_molecule_info["subtype"] = "cyclic_nucleotide"
        category = "small_molecule"
    else:
        category, peptide_info = _detect_peptide(text)
        if category != "peptide" and "pna" not in base_clean.lower() and _has_oligo_signal(text):
            text, oligo_info, extra = parse_oligo_segments(text, flags)
            flags.update(extra)
            category = "oligonucleotide"
        text = _detect_and_remove(text, "fluorophore", flags)

    text = _cleanup(text)

    status = ""
    flag_empty_after_clean = False
    if category != "oligonucleotide" and not text:
        text = _cleanup(base_clean)
        if not text:
            text = _cleanup(_unicode_normalize(name))
        status = "empty_after_clean"
        flag_empty_after_clean = True
        logger.warning("Name empty after cleaning; using fallback: %r", name)

    if category == "oligonucleotide":
        seqs = oligo_info.get("sequences", [])
        length = seqs[0]["length"] if seqs else 0
        typ = oligo_info.get("type", "unknown").lower()
        if typ == "sirna":
            normalized_name = f"sirna {length}mer sense/antisense"
        elif typ == "crispr":
            normalized_name = f"crispr grna {length}mer"
        else:
            normalized_name = f"{typ} {length}mer"
        normalized_name = normalized_name.strip()
        search_name = normalized_name
    else:
        if small_molecule_info.get("subtype") == "cyclic_nucleotide":
            normalized_name = _fix_spacing(text)
        else:
            normalized_name = _fix_spacing(text).lower()
        search_name = normalized_name

    removed_tokens_flat = _flatten_flags(flags)
    oligo_tokens_flat = _flatten_oligo_flags(flags)

    result = {
        "normalized_name": normalized_name,
        "search_name": search_name,
        "search_override_reason": "",
        "category": category,
        "peptide_info": peptide_info,
        "oligo_info": oligo_info,
        "small_molecule_info": small_molecule_info,
        "flags": flags,
        "removed_tokens_flat": removed_tokens_flat,
        "oligo_tokens_flat": oligo_tokens_flat,
        "status": status,
        "flag_isotope": bool(flags.get("isotope")),
        "flag_fluorophore": bool(flags.get("fluorophore")),
        "flag_biotin": bool(flags.get("biotin")),
        "flag_salt": bool(flags.get("salt")),
        "flag_hydrate": bool(flags.get("hydrate")),
        "flag_chromophore": bool(flags.get("chromophore")),
        "flag_oligo": bool(flags.get("oligo")),
        "flag_empty_after_clean": flag_empty_after_clean,
    }
    return result
