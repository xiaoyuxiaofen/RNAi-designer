from __future__ import annotations

import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

from .fasta import FastaRecord


@dataclass(frozen=True)
class AccessibilityConfig:
    method: str = "heuristic"
    executable: str = "RNAplfold"
    window: int = 80
    span: int = 40
    unpaired_length: int = 21
    require_external: bool = False


def compute_accessibility(record: FastaRecord, config: AccessibilityConfig) -> dict[int, float]:
    if config.method == "none":
        return {}
    if config.method == "heuristic":
        return heuristic_accessibility(record.sequence, config.unpaired_length)
    if config.method == "rnaplfold":
        if not rnaplfold_available(config.executable):
            python_values = run_viennarna_python(record, config)
            if python_values is not None:
                return python_values
            if config.require_external:
                raise RuntimeError(
                    f"Neither RNAplfold executable nor ViennaRNA Python bindings were found: {config.executable}"
                )
            return heuristic_accessibility(record.sequence, config.unpaired_length)
        return run_rnaplfold(record, config)
    raise ValueError(f"Unknown accessibility method: {config.method}")


def rnaplfold_available(executable: str = "RNAplfold") -> bool:
    return shutil.which(executable) is not None


def heuristic_accessibility(sequence: str, unpaired_length: int) -> dict[int, float]:
    """Fallback accessibility estimate used only when RNAplfold is unavailable."""

    values: dict[int, float] = {}
    for start in range(0, len(sequence) - unpaired_length + 1):
        window = sequence[start : start + unpaired_length]
        au_fraction = sum(1 for base in window if base in {"A", "T", "U"}) / unpaired_length
        gc_fraction = sum(1 for base in window if base in {"G", "C"}) / unpaired_length
        values[start + 1] = max(0.0, min(1.0, 0.20 + 0.75 * au_fraction - 0.15 * gc_fraction))
    return values


def run_rnaplfold(record: FastaRecord, config: AccessibilityConfig) -> dict[int, float]:
    with TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        command = [
            config.executable,
            "-W",
            str(config.window),
            "-L",
            str(config.span),
            "-u",
            str(config.unpaired_length),
        ]
        fasta_text = f">{record.id}\n{record.sequence.replace('T', 'U')}\n"
        subprocess.run(command, input=fasta_text, text=True, cwd=tmp_path, check=True, capture_output=True)
        lunp_files = sorted(tmp_path.glob("*_lunp"))
        if not lunp_files:
            raise RuntimeError("RNAplfold did not create a *_lunp output file")
        return parse_lunp(lunp_files[0].read_text(encoding="utf-8"), config.unpaired_length)


def run_viennarna_python(record: FastaRecord, config: AccessibilityConfig) -> dict[int, float] | None:
    _add_local_python_packages()
    try:
        import RNA  # type: ignore[import-not-found]
    except ImportError:
        return None

    matrix = RNA.pfl_fold_up(
        record.sequence.replace("T", "U"),
        config.unpaired_length,
        config.window,
        config.span,
    )
    values: dict[int, float] = {}
    for position in range(1, len(record.sequence) - config.unpaired_length + 2):
        try:
            values[position] = max(0.0, min(1.0, float(matrix[position][config.unpaired_length])))
        except (IndexError, TypeError):
            values[position] = 0.0
    return values


def parse_lunp(text: str, unpaired_length: int) -> dict[int, float]:
    values: dict[int, float] = {}
    selected_column = unpaired_length
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        columns = line.split()
        if len(columns) <= selected_column:
            continue
        try:
            position = int(columns[0])
            value = float(columns[selected_column])
        except ValueError:
            continue
        values[position] = max(0.0, min(1.0, value))
    return values


def _add_local_python_packages() -> None:
    project_root = Path(__file__).resolve().parents[2]
    local_packages = project_root / "tools" / "python-packages"
    if local_packages.exists():
        local_path = str(local_packages)
        if local_path not in sys.path:
            sys.path.insert(0, local_path)
