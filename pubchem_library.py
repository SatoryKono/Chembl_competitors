from __future__ import annotations

"""Helpers for assembling compound name dictionaries.

This module provides :func:`build_compound_name_dictionary` which mirrors the
behaviour of an internal Power Query script used for normalising chemical
compound names.  Only the parts required for building a synonym dictionary are
reimplemented here in pure pandas.
"""

import logging
from typing import Iterable, Optional, Sequence

import pandas as pd

logger = logging.getLogger(__name__)

__all__ = ["build_compound_name_dictionary"]

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

DEFAULT_EXCLUSIONS: Sequence[str] = [
    "(+)-pentazocine",
    "(r)-qnb",
    "af-dx 384",
    "crotonyl-coa",
    "dd-coa",
    "dodecenoyl-coa",
    "ethylketazocine",
    "gtpgammas",
    "l-alpha-phosphatidylinositol",
    "l-dihydroorotate",
    "l-glyceraldehyde",
    "leu-mca",
    "r-pia",
    "rr-src",
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _normalise_str(s: object) -> str:
    """Return ``s`` as a trimmed string; ``None`` becomes ``""``."""

    if pd.isna(s):
        return ""
    return str(s).strip()


def _split_synonyms(value: str) -> list[str]:
    """Split the pipe separated synonym string."""

    if not value:
        return []
    return [part.strip() for part in value.split("|") if part.strip()]


def _parse_inchikey(key: str) -> tuple[Optional[str], Optional[str], Optional[str]]:
    """Split an InChIKey into its three sections.

    Parameters
    ----------
    key:
        Raw InChIKey string.
    """

    if not key:
        return None, None, None
    parts = key.strip().upper().split("-")
    if len(parts) != 3:
        return None, None, None
    return parts[0], parts[1], parts[2]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def build_compound_name_dictionary(
    df: pd.DataFrame,
    *,
    name_col: str = "search_name",
    pubchem_cid_col: str = "pubchem_cid",
    canonical_smiles_col: str = "canonical_smiles",
    inchi_col: str = "inchi",
    inchikey_col: str = "inchi_key",
    iupac_name_col: str = "iupac_name",
    synonyms_col: str = "synonyms",
    exclusions: Optional[Iterable[str]] = None,
) -> pd.DataFrame:
    """Build a synonym dictionary from annotated PubChem data.

    Parameters
    ----------
    df:
        Source table containing at least the columns required by the function.
    name_col, pubchem_cid_col, canonical_smiles_col, inchi_col, inchikey_col,
    iupac_name_col, synonyms_col:
        Column names matching the respective fields in ``df``.
    exclusions:
        Optional iterable of names that should be ignored when building the
        isotope mapping pool.  Comparison is case-insensitive.  When ``None``
        a built in default list is used.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns ``["synonyms", "preferent_name", "pubchem_cid",
        "inchi_key", "merge_index", "reference_synonyms"]``.

    Raises
    ------
    ValueError
        If any of the required columns are missing from ``df``.
    """

    required = {
        name_col,
        pubchem_cid_col,
        canonical_smiles_col,
        inchi_col,
        inchikey_col,
        iupac_name_col,
        synonyms_col,
    }
    missing = required.difference(df.columns)
    if missing:
        msg = f"Missing required columns: {', '.join(sorted(missing))}"
        logger.error(msg)
        raise ValueError(msg)

    df = df.copy()

    # ------------------------------------------------------------------
    # Sanitize input types
    # ------------------------------------------------------------------
    for col in [
        name_col,
        canonical_smiles_col,
        inchi_col,
        inchikey_col,
        iupac_name_col,
        synonyms_col,
    ]:
        df[col] = df[col].map(_normalise_str)
    df[pubchem_cid_col] = df[pubchem_cid_col].map(_normalise_str)

    # ------------------------------------------------------------------
    # Normalisation and synonym filtering
    # ------------------------------------------------------------------
    df[synonyms_col] = df[synonyms_col].str.lower()
    sentinels = {"compound name is too short", "multiply", "unknown"}
    df = df[~df[synonyms_col].isin(sentinels)].copy()
    df["_synonyms_list"] = df[synonyms_col].apply(_split_synonyms)

    # ------------------------------------------------------------------
    # Aggregate all synonyms
    # ------------------------------------------------------------------
    def build_all_syn(row: pd.Series) -> list[str]:
        items: list[str] = []
        items.append(_normalise_str(row[name_col]).lower())
        items.extend(row["_synonyms_list"])
        items.append(_normalise_str(row[iupac_name_col]).lower())
        items = [i for i in items if i]
        return list(dict.fromkeys(items))

    df["all_synonyms"] = df.apply(build_all_syn, axis=1)

    # ------------------------------------------------------------------
    # InChIKey parsing
    # ------------------------------------------------------------------
    key_parts = df[inchikey_col].str.upper().apply(_parse_inchikey)
    df["inchi_key1"], df["inchi_key2"], df["inchi_key3"] = zip(*key_parts)

    # ------------------------------------------------------------------
    # Base columns and deterministic index
    # ------------------------------------------------------------------
    df = df[
        [
            name_col,
            pubchem_cid_col,
            inchikey_col,
            "all_synonyms",
            "inchi_key1",
            "inchi_key2",
            "inchi_key3",
        ]
    ].rename(
        columns={
            name_col: "preferent_name",
            pubchem_cid_col: "pubchem_cid",
            inchikey_col: "inchi_key",
        }
    )
    df = df.reset_index(drop=True)
    df["Index"] = df.index

    # ------------------------------------------------------------------
    # Count per InChIKey1 for entries with InChIKey3 == 'N'
    # ------------------------------------------------------------------
    dedup = df.drop_duplicates(subset=["pubchem_cid"])
    mask = (
        dedup["inchi_key3"].eq("N")
        & dedup["inchi_key1"].notna()
        & dedup["inchi_key1"].ne("")
    )
    counts = (
        dedup[mask].groupby("inchi_key1").size().rename("Count").reset_index()
    )
    df = df.merge(counts, on="inchi_key1", how="left")
    df["Count"] = df["Count"].fillna(0).astype(int)
    df = df.sort_values(["Count", "inchi_key1"], ascending=[False, False])

    # ------------------------------------------------------------------
    # Candidate selection
    # ------------------------------------------------------------------
    df["preferent_name_lc"] = df["preferent_name"].str.lower()
    exc = {
        s.lower() for s in (exclusions if exclusions is not None else DEFAULT_EXCLUSIONS)
    }
    iso_candidates = df[df["Count"].isin([2, 3])].copy()
    tmp_pool = df[
        (df["Count"] != 1)
        & df["inchi_key3"].eq("N")
        & ~df["preferent_name_lc"].isin(exc)
    ][[
        "inchi_key1",
        "preferent_name",
        "pubchem_cid",
        "inchi_key",
        "Index",
    ]]

    iso_joined = iso_candidates.merge(
        tmp_pool,
        on="inchi_key1",
        how="left",
        suffixes=("_iso", ""),
    )
    iso_joined = iso_joined.rename(
        columns={
            "preferent_name_iso": "iso_preferent_name",
            "pubchem_cid_iso": "iso_pubchem_cid",
            "inchi_key_iso": "iso_inchi_key",
            "Index": "merge_index",
            "Index_iso": "Index",
        }
    )

    # ------------------------------------------------------------------
    # Singleton handling
    # ------------------------------------------------------------------
    singles = df[df["Count"] == 1].copy()
    singles["iso_preferent_name"] = singles["preferent_name"]
    singles["iso_pubchem_cid"] = singles["pubchem_cid"]
    singles["iso_inchi_key"] = singles["inchi_key"]
    singles["merge_index"] = singles["Index"]

    # ------------------------------------------------------------------
    # Union
    # ------------------------------------------------------------------
    combined = pd.concat([iso_joined, singles], ignore_index=True, sort=False)

    combined["synonyms"] = combined["all_synonyms"].apply(lambda v: "|".join(v) if v else "")
    combined = combined.drop(columns=["all_synonyms"])
    combined["merge_index"] = combined["merge_index"].fillna(combined["Index"])

    # Cleanup temporary columns
    combined = combined.drop(
        columns=[
            "preferent_name_lc",
            "Count",
            "inchi_key1",
            "inchi_key2",
            "inchi_key3",
            "Index",
            "iso_preferent_name",
            "iso_pubchem_cid",
            "iso_inchi_key",
        ]
    )

    # ------------------------------------------------------------------
    # Grouping and expansion
    # ------------------------------------------------------------------
    def agg_syns(series: pd.Series) -> list[str]:
        items: list[str] = []
        for value in series.dropna():
            items.extend(_split_synonyms(value))
        return list(dict.fromkeys(items))

    grouped = (
        combined.groupby(
            ["preferent_name", "pubchem_cid", "inchi_key", "merge_index"],
            dropna=False,
        )["synonyms"].apply(agg_syns).reset_index()
    )
    grouped["reference_synonyms"] = grouped["synonyms"].apply(lambda lst: "|".join(lst))
    final_df = grouped.explode("synonyms").reset_index(drop=True)
    final_df = final_df[
        [
            "synonyms",
            "preferent_name",
            "pubchem_cid",
            "inchi_key",
            "merge_index",
            "reference_synonyms",
        ]
    ]
    # Remove any accidental duplicate synonym rows that might have been
    # produced during aggregation.  Only the first occurrence is kept to match
    # the behaviour of the original Power Query script.
    final_df = final_df.drop_duplicates(subset=["synonyms"]).reset_index(
        drop=True
    )

    return final_df
