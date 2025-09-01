# Chemical Name Normalizer

This project provides utilities to preprocess and normalize chemical names in bulk.

## Installation

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
- `category`
- `peptide_info`
- `flags`
- `removed_tokens_flat`
- `flag_isotope`
- `flag_fluorophore`
- `flag_biotin`
- `flag_salt`
- `flag_hydrate`

Isotopic labels such as `[3H]`, `14C`, `d5`, or `U-13C` are removed from the
normalized name and logged under `flags.isotope`.

Salts and mineral acids such as hydrochloride, HCl, HBr, HNO3 or H2SO4 are
removed from the normalized name and logged under `flags.salt`. Flattened
tokens are also provided in `removed_tokens_flat` using a
`<flag>:<token>|<flag>:<token>` format.

## License

MIT
