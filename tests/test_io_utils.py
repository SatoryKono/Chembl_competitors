"""Tests for I/O utilities."""

from __future__ import annotations

from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parents[1]))

from mylib.io_utils import read_input_csv


def test_read_input_csv_unescaped_commas(tmp_path: Path) -> None:
    """Unescaped commas should be handled via line-based fallback."""

    data = "input_name\nN,N-dimethyltryptamine\naspirin\n"
    csv_path = tmp_path / "broken.csv"
    csv_path.write_text(data)

    df = read_input_csv(csv_path)
    assert df["input_name"].tolist() == ["N,N-dimethyltryptamine", "aspirin"]

