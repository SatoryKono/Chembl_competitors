# Chembl Competitors

Utilities for loading and normalising chemical compound names.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash

python main.py --input examples.csv --output normalized_competitor.csv
```

If ``--output`` is omitted the normalised table is printed to ``stdout``.



## Development

Formatting and static checks:

```bash
black .
ruff .
mypy .
```

## Testing

```bash
pytest -q
```
