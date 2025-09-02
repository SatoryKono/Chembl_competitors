"""Utilities for querying PubChem by compound name.

This module exposes functions to look up PubChem Compound IDs (CIDs)
using the PUG REST API. An exact name match is performed according to
project requirements.
"""

from __future__ import annotations

import logging
from typing import Optional
from urllib.parse import quote_plus

import requests
import pandas as pd

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PUBCHEM_NAME_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/cids/TXT"
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

    Raises
    ------
    requests.RequestException
        On network or HTTP errors other than a 404 response.
    """

    if len(name) < 5:
        logger.debug("Skipping lookup for short name: %s", name)
        return "compound name is too short"

    sess = session or requests.Session()
    url = PUBCHEM_NAME_URL.format(quote_plus(name))
    params = {"name_type": "exact"}
    logger.debug("Requesting PubChem CID for %s", name)
    try:
        response = sess.get(url, params=params, timeout=10)
    except requests.RequestException:
        logger.exception("Failed to query PubChem for %s", name)
        raise

    if response.status_code == 404:
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
    df["pubchem_cid"] = df[name_column].apply(lambda n: fetch_pubchem_cid(str(n), session=sess))
    return df
