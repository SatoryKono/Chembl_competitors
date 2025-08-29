from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parent.parent))

from mylib.io_utils import smart_read_csv


def test_smart_read_csv_detects_cp1251_and_semicolon():
    path = Path(__file__).parent / "data" / "sample_cp1251_semicolon.csv"
    df, enc, sep = smart_read_csv(path)
    assert sep == ";"
    assert enc.lower().startswith("cp1251")
    # NA token should be converted to empty string
    assert df.loc[1, "input_name"] == ""
