from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

from .fasta import FastaRecord, write_fasta


@dataclass(frozen=True)
class Bowtie2Alignment:
    query_id: str
    transcript_id: str
    position: int
    strand: str
    mismatches: int


def resolve_bowtie2_executable(executable: str = "bowtie2") -> str | None:
    return _resolve_tool(
        executable,
        [
            Path("tools") / "bowtie2" / "bowtie2-2.5.0-mingw-x86_64" / "bowtie2-align-s.exe",
            Path("tools") / "bowtie2" / "bowtie2-2.5.0-mingw-x86_64" / "bowtie2.exe",
        ],
    )


def resolve_bowtie2_build_executable(executable: str = "bowtie2-build") -> str | None:
    return _resolve_tool(
        executable,
        [
            Path("tools") / "bowtie2" / "bowtie2-2.5.0-mingw-x86_64" / "bowtie2-build-s.exe",
            Path("tools") / "bowtie2" / "bowtie2-2.5.0-mingw-x86_64" / "bowtie2-build.exe",
        ],
    )


def build_bowtie2_index(
    transcriptome_fasta: str,
    index_prefix: str,
    executable: str = "bowtie2-build",
) -> None:
    resolved = resolve_bowtie2_build_executable(executable)
    if not resolved:
        raise RuntimeError(f"Bowtie2 build executable not found: {executable}")
    Path(index_prefix).parent.mkdir(parents=True, exist_ok=True)
    subprocess.run([resolved, transcriptome_fasta, index_prefix], check=True, capture_output=True, text=True)


def run_bowtie2(
    queries: dict[str, str],
    bowtie2_index: str,
    executable: str = "bowtie2",
    seed_size: int = 12,
) -> list[Bowtie2Alignment]:
    resolved = resolve_bowtie2_executable(executable)
    if not resolved:
        raise RuntimeError(f"Bowtie2 executable not found: {executable}")
    with TemporaryDirectory() as tmpdir:
        query_path = Path(tmpdir) / "sirnas.fa"
        write_fasta(
            [FastaRecord(id=query_id, description=query_id, sequence=sequence) for query_id, sequence in queries.items()],
            query_path,
        )
        command = [
            resolved,
            "-x",
            bowtie2_index,
            "-f",
            "-U",
            str(query_path),
            "-a",
            "--end-to-end",
            "-N",
            "1",
            "-L",
            str(seed_size),
            "--no-unal",
        ]
        completed = subprocess.run(command, check=True, capture_output=True, text=True)
    return parse_bowtie2_sam(completed.stdout.splitlines())


def parse_bowtie2_sam(lines: list[str]) -> list[Bowtie2Alignment]:
    alignments: list[Bowtie2Alignment] = []
    for line in lines:
        if not line or line.startswith("@"):
            continue
        columns = line.split("\t")
        if len(columns) < 11:
            continue
        flag = int(columns[1])
        if flag & 4:
            continue
        nm = 0
        for column in columns[11:]:
            if column.startswith("NM:i:"):
                nm = int(column.removeprefix("NM:i:"))
                break
        alignments.append(
            Bowtie2Alignment(
                query_id=columns[0],
                transcript_id=columns[2],
                position=int(columns[3]),
                strand="-" if flag & 16 else "+",
                mismatches=nm,
            )
        )
    return alignments


def _resolve_tool(executable: str, local_candidates: list[Path]) -> str | None:
    path = shutil.which(executable)
    if path:
        return path
    explicit = Path(executable)
    if explicit.exists():
        return str(explicit)
    project_root = Path(__file__).resolve().parents[2]
    for candidate in local_candidates:
        full_path = project_root / candidate
        if full_path.exists():
            return str(full_path)
    return None
