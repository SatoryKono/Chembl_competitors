"""Utilities for querying PubChem by compound name.

This module exposes functions to look up PubChem Compound IDs (CIDs)
using the PUG REST API. An exact name match is attempted first; if no
result is returned the function falls back to a broader search to avoid
false negatives for valid synonyms.
"""

from __future__ import annotations

import logging
from typing import Dict, Optional
from urllib.parse import quote

import requests
import pandas as pd
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PUBCHEM_NAME_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/cids/TXT"
# JSON endpoints for richer metadata
PUBCHEM_PROPERTY_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/property/"
    "CanonicalSMILES,InChI,InChIKey,MolecularFormula,MolecularWeight,IUPACName/JSON"
)
PUBCHEM_SYNONYM_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/synonyms/JSON"
)


# ---------------------------------------------------------------------------
# Session helper
# ---------------------------------------------------------------------------


def _get_session(session: Optional[requests.Session] = None) -> requests.Session:
    """Return a requests session with retry logic.

    Parameters
    ----------
    session:
        Existing session to reuse. If ``None`` a new session configured with
        retries for common transient network errors is created.

    Returns
    -------
    requests.Session
        Session instance ready for use.
    """

    if session is not None:
        return session

    retries = Retry(
        total=3,
        backoff_factor=0.3,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
    )
    adapter = HTTPAdapter(max_retries=retries)
    sess = requests.Session()
    sess.mount("http://", adapter)
    sess.mount("https://", adapter)
    return sess


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

    On network errors the function logs the exception and returns ``"unknown"``.
    """

    if len(name) < 5:
        logger.debug("Skipping lookup for short name: %s", name)
        return "compound name is too short"

    sess = _get_session(session)
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
        logger.exception("Failed to query PubChem for %s", name)
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


def fetch_pubchem_record(
    name: str, *, session: Optional[requests.Session] = None
) -> Dict[str, str]:
    """Return metadata for ``name`` from PubChem.

    The lookup first resolves ``name`` to a CID via :func:`fetch_pubchem_cid` and
    subsequently queries the compound's properties and synonyms. The following
    keys are returned:

    ``pubchem_cid``, ``canonical_smiles``, ``inchi``, ``inchi_key``,
    ``molecular_formula``, ``molecular_weight``, ``iupac_name``, ``synonyms``.

    If the CID lookup yields one of the sentinel messages (``"unknown"``,
    ``"multiply"``, ``"compound name is too short"``) that value is propagated
    to all fields.

    Parameters
    ----------
    name:
        Compound name to query.
    session:
        Optional :class:`requests.Session` for connection pooling.

    Returns
    -------
    dict
        Mapping of field name to value as described above. Network failures
        during property or synonym retrieval are logged and yield empty strings.
    """

    cid = fetch_pubchem_cid(name, session=session)
    if not cid.isdigit():
        # Propagate sentinel value across all fields
        return {
            "pubchem_cid": cid,
            "canonical_smiles": cid,
            "inchi": cid,
            "inchi_key": cid,
            "molecular_formula": cid,
            "molecular_weight": cid,
            "iupac_name": cid,
            "synonyms": cid,
        }

    sess = _get_session(session)

    # ------------------------------------------------------------------
    # Fetch compound properties
    # ------------------------------------------------------------------
    prop_data: dict[str, str] = {}
    try:
        prop_resp = sess.get(PUBCHEM_PROPERTY_URL.format(cid), timeout=10)
        if prop_resp.status_code not in {400, 404}:
            prop_resp.raise_for_status()
            try:
                prop_json = prop_resp.json()
            except ValueError:
                prop_json = {}
            props = (
                prop_json.get("PropertyTable", {})
                .get("Properties", [{}])[0]
            )
            prop_data = {
                k: str(props.get(k, ""))
                for k in [
                    "CanonicalSMILES",
                    "InChI",
                    "InChIKey",
                    "MolecularFormula",
                    "MolecularWeight",
                    "IUPACName",
                ]
            }
        else:
            logger.info("No properties found for CID %s", cid)
    except requests.RequestException:
        logger.exception("Failed to fetch properties for CID %s", cid)
        prop_data = {}
    except Exception:
        logger.exception("Unexpected error fetching properties for CID %s", cid)
        prop_data = {}

    # ------------------------------------------------------------------
    # Fetch synonyms
    # ------------------------------------------------------------------
    synonyms = ""
    try:
        syn_resp = sess.get(PUBCHEM_SYNONYM_URL.format(cid), timeout=10)
        if syn_resp.status_code not in {400, 404}:
            syn_resp.raise_for_status()
            try:
                syn_json = syn_resp.json()
            except ValueError:
                syn_json = {}
            syn_list = (
                syn_json.get("InformationList", {})
                .get("Information", [{}])[0]
                .get("Synonym", [])
            )
            synonyms = "|".join(syn_list)
        else:
            logger.info("No synonyms found for CID %s", cid)
    except requests.RequestException:
        logger.exception("Failed to fetch synonyms for CID %s", cid)
        synonyms = ""
    except Exception:
        logger.exception("Unexpected error fetching synonyms for CID %s", cid)
        synonyms = ""

    return {
        "pubchem_cid": cid,
        "canonical_smiles": prop_data.get("CanonicalSMILES", ""),
        "inchi": prop_data.get("InChI", ""),
        "inchi_key": prop_data.get("InChIKey", ""),
        "molecular_formula": prop_data.get("MolecularFormula", ""),
        "molecular_weight": prop_data.get("MolecularWeight", ""),
        "iupac_name": prop_data.get("IUPACName", ""),
        "synonyms": synonyms,
    }


def annotate_pubchem_info(
    df: pd.DataFrame,
    *,
    name_column: str = "search_name",
    session: Optional[requests.Session] = None,
) -> pd.DataFrame:
    """Annotate a DataFrame with PubChem metadata.

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
        Copy of ``df`` with additional PubChem columns described in
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

    sess = _get_session(session)
    logger.debug("Annotating %d records with PubChem metadata", len(df))
    df = df.copy()
    info = df[name_column].apply(lambda n: fetch_pubchem_record(str(n), session=sess))
    info_df = pd.DataFrame(list(info))
    return pd.concat([df, info_df], axis=1)


# Backwards compatibility ---------------------------------------------------
annotate_pubchem_cids = annotate_pubchem_info
