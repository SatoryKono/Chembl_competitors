from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parent.parent))

from mylib.transforms import normalize_name


def test_stereochemistry_not_flagged_as_isotope():
    name = "(4s)-4,5-dihydro-2-(6-hydroxybenzothiazolyl)-4-thiazolecarboxylic acid"
    norm, flags = normalize_name(name)
    assert not flags["isotope"], f"Unexpected isotope flag: {flags['isotope']}"
    assert "labeled" not in norm.lower()


def test_isotope_detected():
    name = "(33p)gammaatp"
    norm, flags = normalize_name(name)
    assert "33P" in flags["isotope"]
    assert norm == "gammaatp"


def test_fluorophore_and_noise_detection():
    name = "fitc-labeled (+)-jq1"
    norm, flags = normalize_name(name)
    assert "fitc" in flags["fluorophore"]
    assert "labeled" in flags["noise"]
    assert "jq1" in norm


def test_isotope_hyphen_cleanup():
    """Isotope tags should not leave stray hyphens in the name."""

    name = "(+)-[3h]-3-ppp"
    norm, flags = normalize_name(name)
    assert norm == "(+)-3-ppp"
    assert "3H" in flags["isotope"]

    name2 = "[3h]-(+)-pentazocine"
    norm2, flags2 = normalize_name(name2)
    assert norm2 == "(+)-pentazocine"
    assert "3H" in flags2["isotope"]
