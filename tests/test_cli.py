from pathlib import Path
import subprocess
import sys

sys.path.append(str(Path(__file__).resolve().parent.parent))

from mylib.io_utils import smart_read_csv


def test_cli_writes_output(tmp_path: Path) -> None:
    input_path = Path(__file__).parent / "data" / "sample_cp1251_semicolon.csv"
    output_path = tmp_path / "out.csv"
    subprocess.run(
        [
            sys.executable,
            "main.py",
            "--input",
            str(input_path),
            "--output",
            str(output_path),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert output_path.exists()
    df, _, _ = smart_read_csv(output_path)
    assert "normalized_name" in df.columns
    assert df.loc[1, "input_name"] == ""
