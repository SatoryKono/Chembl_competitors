"""Utilities for querying PubChem by compound name.

This module exposes helpers for looking up chemical information via the
`PUG REST`_ API. An exact name match is attempted first; if no result is
returned the function falls back to a broader search to avoid false
negatives for valid synonyms.

The higher level :func:`fetch_pubchem_record` function retrieves both
basic compound properties and known synonyms.  It is designed to degrade
gracefully when the remote service is unavailable or returns partial
information.

.. _PUG REST: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest
"""

from __future__ import annotations

import logging
from typing import Dict, Optional
from urllib.parse import quote

import requests
import pandas as pd

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PUBCHEM_NAME_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/cids/TXT"
)
PUBCHEM_PROPERTY_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/property/"
    "CanonicalSMILES,InChI,InChIKey,MolecularFormula,MolecularWeight,IUPACName/JSON"
)
PUBCHEM_SYNONYM_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/synonyms/JSON"
)


# ---------------------------------------------------------------------------
# Core functionality
# ---------------------------------------------------------------------------

def fetch_pubchem_cid(name: str, *, session: Optional[requests.Session] = None) -> str:
    """Return PubChem CID for ``name`` using an exact match.

    The function implements the following decision logic:

    * If ``name`` has fewer than five characters the lookup is skipped and
      the string ``"compound name is too short"`` is returned.
    * If the PUG REST service yields no CID ``"unknown"`` is returned.
    * If more than one CID is returned ``"multiply"`` is returned.
    * Otherwise the single CID string is returned.

    Parameters
    ----------
    name:
        Compound name to query.
    session:
        Optional :class:`requests.Session` for connection pooling. If not
        provided a temporary session is created for this call.

    Returns
    -------
    str
        CID string or one of the sentinel messages described above.

    Notes
    -----
    Network errors are logged and return ``"unknown"`` instead of raising an
    exception.  This mirrors the behaviour of the original Power Query
    implementation this project is derived from.
    """

    if len(name) < 5:
        logger.debug("Skipping lookup for short name: %s", name)
        return "compound name is too short"

    sess = session or requests.Session()
    url = PUBCHEM_NAME_URL.format(quote(name))
    logger.debug("Requesting PubChem CID for %s", name)
    try:
        # First attempt an exact name match.
        response = sess.get(url, params={"name_type": "exact"}, timeout=10)
        if response.status_code in {400, 404}:
            # Fallback: retry without the exact-match constraint which may
            # yield results for valid synonyms not recognised as exact names.
            logger.debug(
                "Exact lookup failed for %s, retrying without name_type", name
            )
            response = sess.get(url, timeout=10)
    except requests.RequestException:
        logger.warning("Failed to query PubChem for %s", name)
        return "unknown"

    # PubChem returns HTTP 404 or 400 when no compound matches the query.
    if response.status_code in {400, 404}:
        logger.info("No PubChem entry found for %s", name)
        return "unknown"

    response.raise_for_status()
    lines = [line.strip() for line in response.text.splitlines() if line.strip()]
    if not lines:
        return "unknown"
    if len(lines) > 1:
        return "multiply"
    return lines[0]


def annotate_pubchem_cids(
    df: pd.DataFrame, *, name_column: str = "search_name", session: Optional[requests.Session] = None
) -> pd.DataFrame:
    """Annotate a DataFrame with PubChem CIDs.

    Parameters
    ----------
    df:
        Input DataFrame containing a column with compound names.
    name_column:
        Column name holding the compounds to search in PubChem.
    session:
        Optional :class:`requests.Session` for connection pooling.

    Returns
    -------
    pandas.DataFrame
        Copy of ``df`` with an additional ``pubchem_cid`` column.

    Raises
    ------
    ValueError
        If ``name_column`` is missing from ``df``.
    requests.RequestException
        Propagated from :func:`fetch_pubchem_cid` on network failures.
    """

    if name_column not in df.columns:
        msg = f"Missing required column '{name_column}'"
        logger.error(msg)
        raise ValueError(msg)

    sess = session or requests.Session()
    logger.debug("Annotating %d records with PubChem CIDs", len(df))
    df = df.copy()
    df["pubchem_cid"] = df[name_column].apply(
        lambda n: fetch_pubchem_cid(str(n), session=sess)
    )
    return df


def fetch_pubchem_record(
    name: str, *, session: Optional[requests.Session] = None
) -> Dict[str, str]:
    """Retrieve PubChem metadata for ``name``.

    The function first resolves ``name`` to a CID via :func:`fetch_pubchem_cid`.
    It then queries additional endpoints for chemical properties and
    synonyms. Network issues or 400/404 responses are treated as missing
    data and yield empty strings in the resulting mapping.

    Parameters
    ----------
    name:
        Compound name to look up.
    session:
        Optional :class:`requests.Session` for connection pooling.

    Returns
    -------
    dict
        Mapping with keys ``pubchem_cid``, ``canonical_smiles``, ``inchi``,
        ``inchi_key``, ``molecular_formula``, ``molecular_weight``,
        ``iupac_name`` and ``synonyms``.
    """

    cid = fetch_pubchem_cid(name, session=session)
    record: Dict[str, str] = {
        "pubchem_cid": cid,
        "canonical_smiles": "",
        "inchi": "",
        "inchi_key": "",
        "molecular_formula": "",
        "molecular_weight": "",
        "iupac_name": "",
        "synonyms": "",
    }
    if cid in {"compound name is too short", "unknown", "multiply"}:
        return record

    sess = session or requests.Session()

    try:
        prop_resp = sess.get(PUBCHEM_PROPERTY_URL.format(cid), timeout=10)
        if prop_resp.status_code not in {400, 404}:
            prop_resp.raise_for_status()
            props = (
                prop_resp.json()
                .get("PropertyTable", {})
                .get("Properties", [{}])[0]
            )
            record.update(
                {
                    "canonical_smiles": props.get("CanonicalSMILES", ""),
                    "inchi": props.get("InChI", ""),
                    "inchi_key": props.get("InChIKey", ""),
                    "molecular_formula": props.get("MolecularFormula", ""),
                    "molecular_weight": str(props.get("MolecularWeight", "")),
                    "iupac_name": props.get("IUPACName", ""),
                }
            )
    except requests.RequestException:
        logger.warning("Property lookup failed for CID %s", cid)

    try:
        syn_resp = sess.get(PUBCHEM_SYNONYM_URL.format(cid), timeout=10)
        if syn_resp.status_code not in {400, 404}:
            syn_resp.raise_for_status()
            info = (
                syn_resp.json()
                .get("InformationList", {})
                .get("Information", [])
            )
            if info:
                record["synonyms"] = "|".join(info[0].get("Synonym", []))
    except requests.RequestException:
        logger.warning("Synonym lookup failed for CID %s", cid)

    return record


def annotate_pubchem_info(
    df: pd.DataFrame,
    *,
    name_column: str = "search_name",
    session: Optional[requests.Session] = None,
) -> pd.DataFrame:
    """Annotate ``df`` with PubChem metadata.

    Parameters
    ----------
    df:
        Input DataFrame containing a column with compound names.
    name_column:
        Column name holding the compounds to search in PubChem.
    session:
        Optional :class:`requests.Session` for connection pooling.

    Returns
    -------
    pandas.DataFrame
        ``df`` with additional columns as produced by
        :func:`fetch_pubchem_record`.

    Raises
    ------
    ValueError
        If ``name_column`` is missing from ``df``.
    """

    if name_column not in df.columns:
        msg = f"Missing required column '{name_column}'"
        logger.error(msg)
        raise ValueError(msg)

    sess = session or requests.Session()
    logger.debug("Annotating %d records with PubChem metadata", len(df))
    records = df[name_column].apply(
        lambda n: fetch_pubchem_record(str(n), session=sess)
    )
    record_df = pd.DataFrame(list(records))
    record_df = record_df.reset_index(drop=True)
    return pd.concat([df.reset_index(drop=True), record_df], axis=1)
