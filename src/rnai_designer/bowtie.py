from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Sequence

from .fasta import FastaRecord, write_fasta


@dataclass(frozen=True)
class BowtieAlignment:
    query_id: str
    transcript_id: str
    position: int
    strand: str
    mismatches: int


def bowtie_available(executable: str = "bowtie") -> bool:
    return resolve_bowtie_executable(executable) is not None


def resolve_bowtie_executable(executable: str = "bowtie") -> str | None:
    path = shutil.which(executable)
    if path:
        return path
    explicit = Path(executable)
    if explicit.exists():
        return str(explicit)
    project_root = Path(__file__).resolve().parents[2]
    candidates = [
        project_root / "tools" / "bowtie" / "bowtie-1.2" / "bowtie-align-s.exe",
        project_root / "tools" / "bowtie" / "bowtie-1.2" / "bowtie.exe",
    ]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate)
    return None


def resolve_bowtie_build_executable(executable: str = "bowtie-build") -> str | None:
    path = shutil.which(executable)
    if path:
        return path
    explicit = Path(executable)
    if explicit.exists():
        return str(explicit)
    project_root = Path(__file__).resolve().parents[2]
    candidates = [
        project_root / "tools" / "bowtie" / "bowtie-1.2" / "bowtie-build-s.exe",
        project_root / "tools" / "bowtie" / "bowtie-1.2" / "bowtie-build.exe",
    ]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate)
    return None


def run_bowtie_v_mode(
    queries: dict[str, str],
    bowtie_index: str,
    mismatches: int,
    executable: str = "bowtie",
) -> list[BowtieAlignment]:
    """Run Bowtie v-mode for siRNA queries against a prebuilt index."""

    resolved_executable = resolve_bowtie_executable(executable)
    if not resolved_executable:
        raise RuntimeError(f"Bowtie executable not found: {executable}")
    with TemporaryDirectory() as tmpdir:
        query_path = Path(tmpdir) / "sirnas.fa"
        write_fasta(
            [FastaRecord(id=query_id, description=query_id, sequence=sequence) for query_id, sequence in queries.items()],
            query_path,
        )
        command = [
            resolved_executable,
            "-f",
            "-a",
            "-v",
            str(mismatches),
            "--best",
            "--strata",
            bowtie_index,
            str(query_path),
        ]
        completed = subprocess.run(command, check=True, capture_output=True, text=True)
    return parse_bowtie_output(completed.stdout.splitlines())


def parse_bowtie_output(lines: Sequence[str]) -> list[BowtieAlignment]:
    alignments: list[BowtieAlignment] = []
    for line in lines:
        if not line.strip():
            continue
        columns = line.rstrip("\n").split("\t")
        if len(columns) < 8:
            continue
        mismatch_text = columns[7]
        mismatches = 0 if mismatch_text == "" else len(mismatch_text.split(","))
        alignments.append(
            BowtieAlignment(
                query_id=columns[0],
                strand=columns[1],
                transcript_id=columns[2],
                position=int(columns[3]) + 1,
                mismatches=mismatches,
            )
        )
    return alignments
