"""PubChem API client utilities.

This module provides functions to interact with the PubChem REST API.
The implementation is a Python translation of a PowerQuery script.
"""
from __future__ import annotations

from dataclasses import dataclass
import logging
import time
from typing import Any, Dict, List, Optional
from urllib.parse import quote

import requests
from requests import Session
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


import pandas as pd

logger = logging.getLogger(__name__)

# A single shared session with retry/backoff for all HTTP calls.  PubChem
# enforces fairly strict rate limits; the retry configuration helps to recover
# from transient failures such as HTTP 5xx responses.
_retry = Retry(
    total=3,
    backoff_factor=1.0,
    status_forcelist=[429, 500, 502, 503, 504],
    allowed_methods=["GET"],
)
_session: Session = requests.Session()
_session.mount("http://", HTTPAdapter(max_retries=_retry))
_session.mount("https://", HTTPAdapter(max_retries=_retry))


def url_encode(text: str) -> str:
    """URL-encode *text* for safe usage in HTTP requests.

    Parameters
    ----------
    text: str
        The string to encode.

    Returns
    -------
    str
        URL-encoded string.
    """
    return quote(text, safe="")

def _cids_from_identifier_list(data: Dict[str, Any]) -> List[str]:
    """Extract CIDs from a JSON ``IdentifierList`` structure."""

    return [str(cid) for cid in data.get("IdentifierList", {}).get("CID", [])]


def get_cid_from_smiles(smiles: str) -> Optional[str]:
    """Retrieve PubChem CID(s) for a SMILES string.

    Parameters
    ----------
    smiles: str
        SMILES representation of a compound.

    Returns
    -------
    str or None
        Pipe-separated list of CIDs or ``None`` if the structure is
        unknown to PubChem.
    """

    safe_smiles = url_encode(smiles)
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/"
        f"{safe_smiles}/cids/JSON"
    )
    response = make_request(url)
    if not response:
        return None
    cids = _cids_from_identifier_list(response)
    unique_cids = sorted(set(cids))
    return "|".join(unique_cids) if unique_cids else None


def get_cid_from_inchi(inchi: str) -> Optional[str]:
    """Retrieve PubChem CID(s) for an InChI string.

    Parameters
    ----------
    inchi:
        InChI representation of a compound.

    Returns
    -------
    str or None
        Pipe-separated list of CIDs or ``None`` if the structure is
        unknown to PubChem.
    """

    safe_inchi = url_encode(inchi)
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/"
        f"{safe_inchi}/cids/JSON"
    )
    response = make_request(url)
    if not response:
        return None
    cids = _cids_from_identifier_list(response)
    unique_cids = sorted(set(cids))
    return "|".join(unique_cids) if unique_cids else None


def get_cid_from_inchikey(inchikey: str) -> Optional[str]:
    """Retrieve PubChem CID(s) for an InChIKey."""

    safe_inchikey = url_encode(inchikey)
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
        f"{safe_inchikey}/cids/JSON"
    )
    response = make_request(url)
    if not response:
        return None
    cids = _cids_from_identifier_list(response)
    unique_cids = sorted(set(cids))
    return "|".join(unique_cids) if unique_cids else None


def make_request(url: str, delay: float = 3.0) -> Optional[Dict[str, Any]]:
    """Make an HTTP GET request and return parsed JSON.

    The function sleeps for ``delay`` seconds before issuing the request in
    order to respect PubChem rate limits.  A shared session configured with
    retries is used to automatically retry transient failures.

    Parameters
    ----------
    url:
        Endpoint URL to query.
    delay:
        Time in seconds to wait before making the request. Defaults to three
        seconds.

    Returns
    -------
    dict or None
        Parsed JSON response or ``None`` when the request fails, the server
        returns a non-success status code, or the payload cannot be decoded.
    """

    time.sleep(delay)
    try:
        response = _session.get(url, timeout=10)
        if response.status_code == 404:
            logger.warning("Request returned 404 for url %s", url)
            return None
        if response.status_code == 400:
            logger.warning("Request returned 400 for url %s", url)
            return None
        response.raise_for_status()
        try:
            return response.json()
        except ValueError:
            logger.warning("Non-JSON response for url %s", url)
            return None
    except requests.RequestException as exc:  # pragma: no cover - network
        logger.error("HTTP request failed for url %s: %s", url, exc)
        return None


def validate_cid(cid: str) -> Optional[str]:
    """Validate PubChem CID.

    Parameters
    ----------
    cid: str
        Candidate CID.

    Returns
    -------
    str or None
        ``cid`` if valid, otherwise ``None`` when the identifier is empty or
        represents an invalid placeholder (``"0"`` or ``"-1"``).
    """
    if cid in {"", "0", "-1"}:
        return None
    return cid


def _extract_cids(bindings: List[Dict[str, Any]]) -> List[str]:
    """Extract CIDs from API bindings."""
    cids: List[str] = []
    for item in bindings:
        cid_field = item.get("cid")
        if isinstance(cid_field, dict):
            cid_value = cid_field.get("value", "")
        else:
            cid_value = str(cid_field)
        cid_value = cid_value.replace(
            "http://rdf.ncbi.nlm.nih.gov/pubchem/compound/CID", ""
        )
        if cid_value:
            cids.append(cid_value)
    return cids


def get_cid(compound_name: str) -> Optional[str]:
    """Retrieve PubChem CID(s) for *compound_name* (exact match).

    Parameters
    ----------
    compound_name: str
        Compound name to query.

    Returns
    -------
    str or None
        Pipe-separated list of CIDs or ``None`` if not found.
    """
    safe_name = url_encode(compound_name)
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/rdf/query?graph=synonym&"
        f"return=cid&format=json&name={safe_name}"
    )
    response = make_request(url)
    if not response:
        return None
    bindings = response.get("results", {}).get("bindings", [])
    cids = _extract_cids(bindings)
    unique_cids = sorted(set(cids))
    return "|".join(unique_cids) if unique_cids else None


def get_all_cid(compound_name: str) -> Optional[str]:
    """Retrieve PubChem CID(s) for *compound_name* (partial match)."""
    safe_name = url_encode(compound_name)
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/rdf/query?graph=synonym&"
        f"return=cid&format=json&name={safe_name}&contain=true"
    )
    response = make_request(url)
    if not response:
        return None
    bindings = response.get("results", {}).get("bindings", [])
    cids = _extract_cids(bindings)
    unique_cids = sorted(set(cids))
    return "|".join(unique_cids) if unique_cids else None


def get_standard_name(cid: str) -> Optional[str]:
    """Retrieve the standard compound name for a given CID."""
    validated = validate_cid(cid)
    if not validated:
        return None
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
        f"{validated}/description/JSON"
    )
    response = make_request(url)
    if not response:
        return None
    info = response.get("InformationList", {}).get("Information", [])
    if not info:
        return None
    return info[0].get("Title")


@dataclass
class Properties:
    """Chemical properties for a PubChem compound."""

    IUPACName: str
    MolecularFormula: str
    iSMILES: str
    cSMILES: str
    InChI: str
    InChIKey: str


def get_properties(cid: str) -> Properties:
    """Retrieve chemical properties for a compound by CID."""
    validated = validate_cid(cid)
    if not validated:
        return Properties(
            "Not Found", "Not Found", "Not Found", "Not Found", "Not Found", "Not Found"
        )
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
        f"{validated}/property/MolecularFormula,IUPACName,IsomericSMILES,"
        "CanonicalSMILES,InChI,InChIKey/JSON"
    )
    response = make_request(url)
    if not response:
        return Properties(
            "Not Found", "Not Found", "Not Found", "Not Found", "Not Found", "Not Found"
        )
    props = response.get("PropertyTable", {}).get("Properties", [])
    if not props:
        return Properties(
            "Not Found", "Not Found", "Not Found", "Not Found", "Not Found", "Not Found"
        )
    item = props[0]
    return Properties(
        item.get("IUPACName", "Not Found"),
        item.get("MolecularFormula", "Not Found"),
        item.get("IsomericSMILES", "Not Found"),
        item.get("CanonicalSMILES", "Not Found"),
        item.get("InChI", "Not Found"),
        item.get("InChIKey", "Not Found"),
    )


def process_compound(compound_name: str) -> Dict[str, str]:
    """Process *compound_name* into a structured record.

    Parameters
    ----------
    compound_name: str
        Name of the compound to look up.

    Returns
    -------
    dict
        Dictionary containing compound details.
    """
    cid = get_cid(compound_name)
    standard = get_standard_name(cid) if cid else None
    props = get_properties(cid) if cid else Properties(
        "Not Found", "Not Found", "Not Found", "Not Found", "Not Found", "Not Found"
    )
    return {
        "Name": compound_name,
        "CID": cid or "Not Found",
        "Standard Name": standard or "Not Found",
        "IUPACName": props.IUPACName,
        "MolecularFormula": props.MolecularFormula,
        "iSMILES": props.iSMILES,
        "cSMILES": props.cSMILES,
        "InChI": props.InChI,
        "InChIKey": props.InChIKey,
    }

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

__all__ = [
    "url_encode",
    "make_request",
    "validate_cid",
    "get_cid",
    "get_all_cid",
    "get_standard_name",
    "get_properties",
    "process_compound",
    "fetch_pubchem_cid", 
    "annotate_pubchem_info",
    "Properties",
]
