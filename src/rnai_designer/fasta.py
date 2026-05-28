from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class FastaRecord:
    """A single FASTA record."""

    id: str
    description: str
    sequence: str


def read_fasta(path: str | Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    current_header: str | None = None
    chunks: list[str] = []

    with Path(path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    records.append(_make_record(current_header, chunks))
                current_header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)

    if current_header is not None:
        records.append(_make_record(current_header, chunks))

    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def write_fasta(records: Iterable[FastaRecord], path: str | Path, line_width: int = 80) -> None:
    with Path(path).open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            header = record.id
            if record.description and record.description != record.id:
                header = f"{record.id} {record.description}"
            handle.write(f">{header}\n")
            for start in range(0, len(record.sequence), line_width):
                handle.write(f"{record.sequence[start:start + line_width]}\n")


def _make_record(header: str, chunks: list[str]) -> FastaRecord:
    record_id = header.split()[0]
    sequence = "".join(chunks).upper().replace("U", "T")
    allowed = {"A", "C", "G", "T", "N"}
    invalid = sorted(set(sequence) - allowed)
    if invalid:
        raise ValueError(f"Record {record_id} contains unsupported bases: {', '.join(invalid)}")
    return FastaRecord(id=record_id, description=header, sequence=sequence)
