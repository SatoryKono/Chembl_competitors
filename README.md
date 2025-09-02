# Chemical Name Normalizer

This project provides utilities to preprocess and normalize chemical names in bulk.

## Installation


During normalization the library canonicalizes spaces around punctuation so
that artifacts such as `5 ' ; 1,3 -diol` become `5'; 1,3-diol`. Brackets lose
padding, connectors (`-`, `/`, `:`, `+`) have no surrounding spaces, semicolons
and commas carry a single space on the right, and primes cling to neighboring
tokens.

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
`flag_empty_after_clean` to `True` so these rows can be reviewed manually.

`search_name` always matches `normalized_name` unless a documented override
occurs. The reason for any override is recorded in `search_override_reason`.

Fluorophore labels such as **Alexa Fluor**, **HiLyte Fluor**, **DyLight**,
**CF** dye series, **Janelia Fluor**, or **BODIPY** families are stripped
early in the pipeline and logged under `flags.fluorophore` so that base
chemical names remain intact.

Chromophore tags like **pNA** are removed and recorded under
`flags.chromophore`, preventing peptide substrates from being mistaken
for oligonucleotides.

Oligonucleotides are recognised via sequence patterns or keywords such as
`oligo`, `primer`, `siRNA`, `gRNA`, and CRISPR-specific terms. The parser
extracts 5′/3′ end modifiers, internal modifications, backbone type, and
roles like sense/antisense or guide/tracr. Matches set `category =
oligonucleotide`, populate `oligo_info`, and record details in
`flags.oligo_*` with a flattened summary in `oligo_tokens_flat`.

Peptides are detected via several heuristics: polymer-style prefixes like
`poly-Glu:Tyr`, explicit terms such as `peptide` or `polypeptide`, and
sequences of one- or three-letter amino-acid codes (optionally bearing
protective groups like `H-`, `Ac-`, `Boc-`, `-OH`, or `-NH2`, as well as
chromophore/fluorophore suffixes such as `-pNA` or `-AMC`. When detected,
the output sets `category = peptide` and populates `peptide_info` with the
peptide type and, for polymer forms, the normalized composition. Generic
materials like "polymer support resin" are not misclassified as peptides
because amino-acid signatures are required.

## License

MIT
