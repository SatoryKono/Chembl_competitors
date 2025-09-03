# Chemical Name Normalizer

This project provides utilities to preprocess and normalize chemical names in bulk.

## PubChem Metadata Lookup

The repository also ships a helper script to annotate compound names with
PubChem metadata. For each value in a column named ``search_name`` the tool
attempts an exact name query against the PubChem PUG REST service. If no
record is found the lookup falls back to a broader synonym search. When a
single CID is resolved the following fields are retrieved:

* ``pubchem_cid``
* ``canonical_smiles``
* ``inchi``
* ``inchi_key``
* ``molecular_formula``
* ``molecular_weight``
* ``iupac_name``
* ``synonyms`` – pipe-separated list

If multiple or no matches are found, or if a name has fewer than five
characters, the corresponding sentinel value (``"multiply"``, ``"unknown`` or
``"compound name is too short"``) is returned for all PubChem columns.

Requests are performed using a session configured with automatic retries so
that transient network failures (e.g., connection resets or HTTP 5xx errors)
do not abort the entire lookup.

### Installation

```bash
pip install -r requirements.txt
```

Optional quality tools:

```bash
pip install black ruff mypy pytest
```

### Usage

```bash
python main.py --input examples1.csv --output out.csv
```

Additional arguments:

- ``--sep`` – CSV delimiter (default `,`).
- ``--encoding`` – file encoding (default `utf-8`).
- ``--log-level`` – logging level.

### Example

Input ``examples1.csv``:

```csv
search_name
aspirin
water
NaCl
```

Run:

```bash
python main.py --input examples1.csv --output out.csv
```

Output ``out.csv``:

```csv
search_name,pubchem_cid,canonical_smiles,inchi,inchi_key,molecular_formula,molecular_weight,iupac_name,synonyms
aspirin,2244,SMILES,InChI,KEY,C9H8O4,180.16,Name,aspirin|acetylsalicylic acid
water,962,,,H2O,18.02,,water
NaCl,compound name is too short,compound name is too short,compound name is too short,compound name is too short,compound name is too short,compound name is too short,compound name is too short,compound name is too short
```

## Installation


During normalization the library canonicalizes spaces around punctuation so
that artifacts such as `5 ' ; 1,3 -diol` become `5'; 1,3-diol`. Brackets lose
padding, connectors (`-`, `/`, `:`, `+`) have no surrounding spaces, semicolons
and commas carry a single space on the right, and primes cling to neighboring
tokens.

Orphaned or empty brackets are removed once all annotations have been stripped.


```bash
pip install -r requirements.txt
```

Optional quality tools:

```bash
pip install black ruff mypy pytest
```

## Usage

```bash
python main.py --input examples1.csv --output out.csv
```

### Arguments

- `--input`: path to input CSV containing an `input_name` column.
- `--output`: path for the output CSV.
- `--sep`: CSV delimiter (default `,`).
- `--encoding`: file encoding (default `utf-8`).
- `--log-level`: logging level.

If the input file contains unescaped commas within chemical names, the loader
falls back to a line-by-line parser. Ensure the first line is the header
`input_name` and each subsequent line contains a single name.

## Development

Recommended commands:

```bash
black .
ruff check .
mypy .
pytest
```

## Testing

Unit tests live under `tests/` and can be run with `pytest`.

## Example

Input:
```
[3H] 8 - oh dpat
biotinylated peptide
```

Output includes columns:
- `normalized_name`
- `search_name`
- `search_override_reason`
- `category`
- `peptide_info`
- `oligo_info`
- `flags`
- `removed_tokens_flat`
- `oligo_tokens_flat`
- `status`
- `flag_isotope`
- `flag_fluorophore`
- `flag_biotin`
- `flag_salt`
- `flag_hydrate`
- `flag_oligo`
- `flag_empty_after_clean`

Isotopic labels such as `[3H]`, `[125I]`, `14C`, `18F`, bare prefixes like `D`
or `T`, `d5`-style deuteration, `deuterated`/`tritiated` descriptors, and
`U-13C` are removed from the normalized name and logged under
`flags.isotope`.

Salts and mineral acids such as hydrochloride, HCl, HBr, HNO3 or H2SO4 are
removed from the normalized name and logged under `flags.salt`. Flattened
tokens are also provided in `removed_tokens_flat` using a
`<flag>:<token>|<flag>:<token>` format.

Non-structural descriptors like `solution`, `soln`, `stock`, `buffer`,
`USP/EP/ACS`, `reagent`, `analytical grade`, `crystalline`, `powder`, or
purity annotations (e.g., `≥95% purity`) are stripped in two passes—first within
parentheses/brackets and then globally—with the removed terms collected under
`flags.noise`.

If aggressive cleaning removes all content from a name, the pipeline falls back
to a minimally cleaned version of the original text. In such cases the output
includes `status = empty_after_clean` and sets the boolean indicator

`flag_empty_after_clean` to `True` so these rows can be reviewed manually. The
same status is raised for short "garbage" tokens (single characters, two-digit
numbers, or a digit followed by `a`/`b`/`c`).


`search_name` always matches `normalized_name` unless a documented override
occurs. The reason for any override is recorded in `search_override_reason`.

Fluorophore labels such as **Alexa Fluor**, **HiLyte Fluor**, **DyLight**,
**CF** dye series, **Janelia Fluor**, or **BODIPY** families are stripped
early in the pipeline and logged under `flags.fluorophore` so that base
chemical names remain intact.

Chromophore tags like **pNA** are removed and recorded under
`flags.chromophore`, preventing peptide substrates from being mistaken
for oligonucleotides.


Terminal fluorophore tags are retained for peptides consisting of a single
amino-acid residue (e.g., `FAM-lys`), though the matched fluorophores are still
logged in `flags.fluorophore`. Internal tags within sequences, such as
`lys(AMC)`, are preserved.


Oligonucleotides are recognised via sequence patterns or keywords such as
`oligo`, `primer`, `siRNA`, `gRNA`, and CRISPR-specific terms. The parser
extracts 5′/3′ end modifiers, internal modifications, backbone type, and
roles like sense/antisense or guide/tracr. Matches set `category =
oligonucleotide`, populate `oligo_info`, and record details in
`flags.oligo_*` with a flattened summary in `oligo_tokens_flat`.

Cyclic nucleotides including `cAMP`, `cGMP`, `c-di-GMP`, and `cGAMP` are
detected by dedicated patterns, normalised to canonical forms (e.g.
`3',5'-cAMP`), and classified as small molecules with
`small_molecule_info.subtype = cyclic_nucleotide`. Any accompanying isotope
or fluorophore tags are stripped and logged before classification, ensuring
they are not mistaken for oligonucleotides.


Additional "guard" classes prevent common metabolites and dyes from being
mistaken for peptides. Standard nucleotides (e.g., `ATP`, `ADP`, `GTP`),
cofactors like `CoA`, cholines, fluorogenic 4-MU glycosides, and dyes such as
phenoxazines or resorufin are recognised early and classified under
`category = small_molecule` with corresponding `small_molecule_info.subtype`
values (`nucleotide`, `cofactor`, `choline`, `fluorogenic_glycoside`, `dye`).

Peptides are detected via several heuristics: polymer-style prefixes like
`poly-Glu:Tyr`, explicit terms such as `peptide` or `polypeptide`, and
sequences of one- or three-letter amino-acid codes (optionally bearing
protective groups like `H-`, `Ac-`, `Boc-`, `-OH`, or `-NH2`, as well as
chromophore/fluorophore suffixes such as `-pNA` or `-AMC`. When detected,
Peptides are detected via several heuristics: polymer-style prefixes like
`poly-Glu:Tyr`, explicit terms such as `peptide` or `polypeptide`, and
sequences of one- or three-letter amino-acid codes (optionally bearing
protective groups like `H-`, `Ac-`, `Boc-`, `-OH`, or `-NH2`). When detected,
the output sets `category = peptide` and populates `peptide_info` with the
peptide type and, for polymer forms, the normalized composition. Generic
materials like "polymer support resin" are not misclassified as peptides
because amino-acid signatures are required.

## License

MIT
