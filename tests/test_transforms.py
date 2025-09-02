"""Unit tests for normalization utilities."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

import pytest

from mylib.transforms import PATTERNS, normalize_name, _fix_spacing


def test_isotope_flag() -> None:
    res = normalize_name("[3H] 8 - oh dpat")
    assert res["flag_isotope"] is True
    assert res["search_name"] == "8-oh dpat"


def test_biotin_and_peptide() -> None:
    res = normalize_name("biotinylated peptide")
    assert res["flag_biotin"] is True
    assert res["category"] == "peptide"
    assert res["search_name"] == "peptide"


def test_sequence_detection() -> None:
    res = normalize_name("Ala-Gly-Ser")
    assert res["category"] == "peptide"
    assert res["peptide_info"]["type"] == "sequence_like"


def test_protective_group_sequence() -> None:
    """Sequences with protective termini are still classified as peptides."""

    res = normalize_name("H-Ala-Gly-OH")
    assert res["category"] == "peptide"
    assert res["peptide_info"]["type"] == "sequence_like"


def test_noise_and_concentration_removal() -> None:
    res = normalize_name("Sample solution 10 mM")
    assert "solution" in res["flags"].get("noise", [])
    assert res["search_name"] == "sample"


def test_parenthetical_noise_two_pass() -> None:
    res = normalize_name("histamine (solution, 10 mM in PBS)")
    noise_flags = res["flags"].get("noise", [])
    assert res["search_name"] == "histamine"
    assert "solution" in noise_flags
    assert "PBS" in noise_flags
    assert any("10 mM" in t for t in noise_flags)


def test_empty_after_clean_status_flag() -> None:
    """Names that empty out after cleaning fallback and set status."""

    res = normalize_name("solution 10 mM")
    assert res["status"] == "empty_after_clean"
    assert res["flag_empty_after_clean"] is True
    assert res["normalized_name"] == "solution 10 mm"
    assert res["search_name"] == "solution 10 mm"


@pytest.mark.parametrize("connector", ["-", "/", "+", ":"])
def test_spacing_compaction_after_flag_removal(connector: str) -> None:
    """Removal of flagged tokens should not leave spaces around connectors."""

    res = normalize_name(f"a FITC {connector} b")
    assert res["search_name"] == f"a{connector}b"


def test_repeated_connector_collapses() -> None:
    res = normalize_name("a - biotin - b")
    assert res["search_name"] == "a-b"


def test_spacing_for_comma_and_decimal() -> None:
    """Spaces around commas and decimals are compacted."""

    res = normalize_name("N , N-dimethyl 1 . 5")
    assert res["search_name"] == "n, n-dimethyl 1.5"


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("5 ' ; 1,3 -diol 1,2- (diol); 1,2- {", "5'; 1,3-diol 1,2-(diol); 1,2-{"),
        ("poly Glu : Tyr", "poly Glu:Tyr"),
        ("8 - oh dpat", "8-oh dpat"),
        ("Na + Cl", "Na+Cl"),
        ("1 ,2 )", "1,2)"),
        ("[{ A } ]", "[{A}]"),
        ("alpha ; beta ;gamma", "alpha; beta; gamma"),
        ("1 ; 2 ; 3", "1; 2; 3"),
        ("1 / 2 / 3", "1/2/3"),
        ("A - B - C", "A-B-C"),
        ("A : B : C", "A:B:C"),
        ("A + B + C", "A+B+C"),
        ("1, 2,3", "1,2,3"),
        ("( A , B )", "(A, B)"),
        ("5 'prime", "5'prime"),
        ("word ′ prime", "word′prime"),
        ("chloro -", "chloro-"),
        ("1,2 -diol", "1,2-diol"),
        ("A ;", "A;"),
        ("(1 ,2 ;3 )", "(1,2; 3)"),
        ("{ 1 , 2 }", "{1,2}"),
    ],
)
def test_canonical_spacing_cases(raw: str, expected: str) -> None:
    """Spacing normalization conforms to canonical punctuation rules."""

    assert _fix_spacing(raw) == expected


def test_salt_tokens_removed_and_logged() -> None:
    res = normalize_name("histamine hydrochloride")
    assert res["search_name"] == "histamine"
    assert res["flags"].get("salt") == ["hydrochloride"]


def test_mineral_acid_tokens() -> None:
    res = normalize_name("dextromethorphan HBr")
    assert res["search_name"] == "dextromethorphan"
    assert res["flags"].get("salt") == ["HBr"]


def test_salt_and_hydrate_combination() -> None:
    res = normalize_name("metformin hydrochloride monohydrate")
    assert res["search_name"] == "metformin"
    assert res["flags"].get("salt") == ["hydrochloride"]
    assert res["flags"].get("hydrate") == ["monohydrate"]


@pytest.mark.parametrize(
    "token",
    [
        "dihydrate",
        "trihydrate",
        "tetrahydrate",
        "pentahydrate",
        "anhydrous",
    ],
)
def test_various_hydrate_tokens(token: str) -> None:
    res = normalize_name(f"glucose {token}")
    assert res["search_name"] == "glucose"
    assert res["flags"].get("hydrate") == [token]


def test_removed_tokens_flat() -> None:
    text = "Alexa Fluor 488 [3H] histamine hydrochloride"
    res = normalize_name(text)
    assert (
        res["removed_tokens_flat"]
        == "fluorophore:Alexa Fluor 488|isotope:[3H]|salt:hydrochloride"
    )


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("] dslet", "dslet"),
        ("[[ampa", "ampa"),
        ("[3H]]-5-ct", "5-ct"),
        ("[ [ 3h ] - progesterone", "progesterone"),
    ],
)
def test_hanging_brackets_removed(raw: str, expected: str) -> None:
    res = normalize_name(raw)
    assert res["search_name"] == expected


@pytest.mark.parametrize("raw", ["9", "14", "1a", "2B", "3 c"])
def test_short_garbage_flag(raw: str) -> None:
    res = normalize_name(raw)
    assert res["status"] == "empty_after_clean"
    assert res["flag_empty_after_clean"] is True



@pytest.mark.parametrize(
    "raw, flag",
    [
        ("FAM-lys", "FAM"),
        ("lys-EDANS", "EDANS"),
    ],
)
def test_single_amino_acid_is_small_molecule(raw: str, flag: str) -> None:
    res = normalize_name(raw)
    assert res["category"] == "small_molecule"
    assert flag in res["flags"].get("fluorophore", [])


def test_removed_tokens_flat_empty() -> None:
    res = normalize_name("aspirin")
    assert res["removed_tokens_flat"] == ""


def test_oligo_tokens_flat_empty() -> None:
    res = normalize_name("aspirin")
    assert res["oligo_tokens_flat"] == ""


def test_search_name_defaults_to_normalized() -> None:
    """search_name should equal normalized_name when no override reason."""

    res = normalize_name("Histamine")
    assert res["search_name"] == res["normalized_name"]
    assert res["search_override_reason"] == ""


@pytest.mark.parametrize(
    "text, expected, tokens",
    [
        ("[3H]-histamine", "histamine", ["[3H]"]),
        ("d5-amphetamine", "amphetamine", ["d5"]),
        ("U-13C glucose", "glucose", ["U-13C"]),
        ("d5 [3H] amphetamine", "amphetamine", ["d5", "[3H]"]),
        ("14C caffeine", "caffeine", ["14C"]),
        ("[14C]caffeine", "caffeine", ["[14C]"]),
        ("125I-insulin", "insulin", ["125I"]),
        ("[125 I] insulin", "insulin", ["[125I]"]),
        ("18F-FDG", "fdg", ["18F"]),
        ("[18F]fluorodeoxyglucose", "fluorodeoxyglucose", ["[18F]"]),
        ("2H water", "water", ["2H"]),
        ("D-amphetamine", "amphetamine", ["D"]),
        ("T-thymidine", "thymidine", ["T"]),
        ("deuterated ethanol", "ethanol", ["deuterated"]),
        ("tritiated thymidine", "thymidine", ["tritiated"]),
        ("d3-deuterated phenol", "phenol", ["d3", "deuterated"]),
        ("[3H][14C] compound", "compound", ["[3H]", "[14C]"]),
        ("d5-125I-amphetamine", "amphetamine", ["d5", "125I"]),
        ("U13C-15N-lysine", "lysine", ["U-13C", "15N"]),
        ("d5 U-13C [3H] sample", "sample", ["d5", "U-13C", "[3H]"]),
        ("[i125]-tyrosine", "tyrosine", ["[125I]"]),
        ("[125-i]tyrosine", "tyrosine", ["[125I]"]),
        ("i-125-tyrosine", "tyrosine", ["125I"]),
        ("iodobenzene[1251]", "iodobenzene", ["[125I]"]),
    ],
)
def test_isotope_variants(text: str, expected: str, tokens: list[str]) -> None:
    """Isotope markers are stripped and logged across representations."""

    res = normalize_name(text)
    assert res["search_name"] == expected
    assert res["flags"].get("isotope") == tokens
    # Ensure no isotopic labels remain after normalization
    assert PATTERNS["isotope"].findall(res["search_name"]) == []


@pytest.mark.parametrize(
    "text, expected, token",
    [
        ("poly-Glu:Tyr Alexa Fluor 488", "poly-glu:tyr", "Alexa Fluor 488"),
        ("HiLyte Fluor 555 peptide", "peptide", "HiLyte Fluor 555"),
        ("DyLight-650 antibody", "antibody", "DyLight-650"),
        ("peptide CF568", "peptide", "CF568"),
        ("Janelia Fluor 549 ligand", "ligand", "Janelia Fluor 549"),
        ("BODIPY-581/591 conjugate", "conjugate", "BODIPY-581/591"),
    ],
)
def test_expanded_fluorophore_tokens(text: str, expected: str, token: str) -> None:
    """Expanded fluorophore patterns are removed prior to other processing."""

    res = normalize_name(text)
    assert res["search_name"] == expected
    assert res["flags"].get("fluorophore") == [token]

@pytest.mark.parametrize(
    "text, composition",
    [
        ("poly-Glu:Tyr", "glu:tyr"),
        ("poly (Glu, Tyr)", "glu:tyr"),
        ("poly Glu Tyr", "glu:tyr"),
    ],
)
def test_poly_peptide_detection(text: str, composition: str) -> None:
    res = normalize_name(text)
    assert res["category"] == "peptide"
    assert res["peptide_info"] == {"type": "polymer", "composition": composition}


def test_polymer_non_peptide() -> None:
    res = normalize_name("polymer support resin")
    assert res["category"] == "small_molecule"


@pytest.mark.parametrize(
    "text, subtype",
    [
        ("10-acetyl-3,7-dihydroxyphenoxazin", "dye"),
        ("acetyl CoA", "cofactor"),
        ("33P-gammaATP", "nucleotide"),
        ("2-mesATP", "nucleotide"),
        ("2-mesADP", "nucleotide"),
        ("acetyl choline", "choline"),
        ("4-mu-glcNAc", "fluorogenic_glycoside"),
        (
            "4-methylumbelliferyl N-acetyl-beta-d-glucosaminide",
            "fluorogenic_glycoside",
        ),
    ],
)
def test_small_molecule_guards(text: str, subtype: str) -> None:
    res = normalize_name(text)
    assert res["category"] == "small_molecule"
    assert res["small_molecule_info"].get("subtype") == subtype


def test_isotope_small_molecule() -> None:
    text = "[3H]-pyrrolidine-2-carboxylic acid biphenyl-2-ylamide"
    res = normalize_name(text)
    assert res["category"] == "small_molecule"
    assert res["flags"].get("isotope") == ["[3H]"]


def test_nucleotide_not_peptide() -> None:
    res = normalize_name("ATP")
    assert res["category"] == "small_molecule"
    assert res["small_molecule_info"].get("subtype") == "nucleotide"


@pytest.mark.parametrize(
    "text", [
        "gly-pro-pna",
        "pyroglu-pro-arg-pna",
        "meo-suc-ala-ala-pro-val-pna",
    ],
)
def test_pna_chromophore_peptides(text: str) -> None:
    res = normalize_name(text)
    assert res["category"] == "peptide"
    assert res["flags"].get("chromophore") == ["pna"]


def test_amc_peptide() -> None:
    res = normalize_name("rhkackac-AMC")
    assert res["category"] == "peptide"
    assert res["flags"].get("fluorophore") == ["AMC"]


def test_plain_peptide_sequence() -> None:
    res = normalize_name("rhkackac")
    assert res["category"] == "peptide"


@pytest.mark.parametrize("seq", ["rrrdddsddd", "grsrsrsrsrsr"])
def test_poly_rich_sequences_are_peptides(seq: str) -> None:
    res = normalize_name(seq)
    assert res["category"] == "peptide"


def test_peptide_with_salt_and_pna() -> None:
    res = normalize_name("h-lys-ala-pna.2HCl")
    assert res["category"] == "peptide"
    assert res["flags"].get("salt") == ["HCl"]
    assert res["flags"].get("chromophore") == ["pna"]


def test_histone_peptide_keyword() -> None:
    res = normalize_name("peptide histone h4 fragment")
    assert res["category"] == "peptide"


def test_from_p_number_peptide() -> None:
    res = normalize_name("rhkackac from p53")
    assert res["category"] == "peptide"


def test_pna_not_oligo() -> None:
    res = normalize_name("PNA ACGTACGT")
    assert res["category"] != "oligonucleotide"


def test_simple_rna_with_mods() -> None:
    res = normalize_name("5'-FAM-ACGUACGUACGU-3'")
    assert res["category"] == "oligonucleotide"
    assert res["oligo_info"]["type"] == "RNA"
    assert res["flags"].get("fluorophore") == ["FAM"]
    assert res["oligo_info"]["mods"]["five_prime"] == ["FAM"]
    assert res["normalized_name"] == "rna 12mer"


def test_sirna_two_strands() -> None:
    text = "siRNA sense ACGUACGUACGUACGUACGU; antisense ACUACGUACGUACGUACGUA"
    res = normalize_name(text)
    assert res["category"] == "oligonucleotide"
    assert res["oligo_info"]["type"] == "siRNA"
    roles = {s["role"] for s in res["oligo_info"]["sequences"]}
    assert {"sense", "antisense"} <= roles
    assert res["normalized_name"].startswith("sirna")


def test_crispr_two_parts() -> None:
    text = "crRNA: GTTTTAGAGCTA; tracrRNA: AAAACCCGGGTT"
    res = normalize_name(text)
    assert res["category"] == "oligonucleotide"
    roles = {s["role"] for s in res["oligo_info"]["sequences"]}
    assert {"crrna", "tracrrna"} <= roles
    assert res["oligo_info"]["type"] == "CRISPR"


def test_aso_ps_backbone() -> None:
    text = "T*G*C*A*T*G*C*A* antisense oligonucleotide PS"
    res = normalize_name(text)
    assert res["category"] == "oligonucleotide"
    assert res["oligo_info"]["type"] == "ASO"
    assert res["oligo_info"]["mods"]["backbone"] == "PS"


def test_vendor_tags_phos_bio() -> None:
    text = "/5Phos/ACGTNNNNACGT/3Bio/"
    res = normalize_name(text)
    assert res["category"] == "oligonucleotide"
    mods = res["oligo_info"]["mods"]
    assert mods["five_prime"] == ["Phos"]
    assert mods["three_prime"] == ["Bio"]
    assert res["flags"].get("biotin") == ["Bio"]
    assert res["normalized_name"] == "dna 12mer"


def test_dna_probe_with_sequence() -> None:
    text = "oligo DNA probe: ACGTNNNNACGT"
    res = normalize_name(text)
    assert res["category"] == "oligonucleotide"
    assert res["oligo_info"]["type"] == "DNA"


def test_cyclic_nucleotide_is_small_molecule() -> None:
    res = normalize_name("3 ' , 5 ' - [ 3h ] camp")
    assert res["category"] == "small_molecule"
    assert res["small_molecule_info"] == {"subtype": "cyclic_nucleotide"}
    assert res["normalized_name"] == "3',5'-cAMP"
    assert res["flags"].get("isotope") == ["[3H]"]


def test_cyclic_nucleotide_fluorophore() -> None:
    res = normalize_name("fam-3 ' , 5 '-camp")
    assert res["category"] == "small_molecule"
    assert res["normalized_name"] == "3',5'-cAMP"
    assert res["flags"].get("fluorophore") == ["FAM"]


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("3 ' , 5 ' - [ 3h ] cgmp", "3',5'-cGMP"),
        ("5 '-cgmp", "5'-cGMP"),
    ],
)
def test_cyclic_nucleotide_variants(raw: str, expected: str) -> None:
    res = normalize_name(raw)
    assert res["category"] == "small_molecule"
    assert res["normalized_name"] == expected


def test_sgrna_context() -> None:
    text = "sgRNA 20nt: GACUACGUACGUACGUACGU"
    res = normalize_name(text)
    assert res["category"] == "oligonucleotide"
    assert res["oligo_info"]["type"] == "CRISPR"
