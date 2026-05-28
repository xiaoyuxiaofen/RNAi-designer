from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from .fasta import FastaRecord


@dataclass(frozen=True)
class OffTargetHit:
    query: str
    transcript_id: str
    position: int
    mismatches: int
    mismatch_positions: tuple[int, ...]
    seed_mismatches: int
    risk_score: float


def scan_offtargets(
    queries: Sequence[str],
    transcriptome: Sequence[FastaRecord],
    target_ids: set[str],
    max_mismatches: int,
    seed_size: int = 12,
) -> list[OffTargetHit]:
    """Scan non-target transcripts for siRNA-like hits with bounded mismatches."""

    if not queries:
        return []
    query_length = len(queries[0])
    if any(len(query) != query_length for query in queries):
        raise ValueError("All off-target queries must have the same length")
    if max_mismatches < 0:
        raise ValueError("max_mismatches must be non-negative")
    if seed_size < 1 or seed_size > query_length:
        raise ValueError("seed_size must be between 1 and query length")

    query_set = set(queries)
    hits: list[OffTargetHit] = []
    if max_mismatches == 0:
        for record in transcriptome:
            if record.id in target_ids:
                continue
            seen_in_record: set[tuple[str, int]] = set()
            for start in range(0, len(record.sequence) - query_length + 1):
                window = record.sequence[start : start + query_length]
                if window in query_set and (window, start) not in seen_in_record:
                    hits.append(OffTargetHit(window, record.id, start + 1, 0, tuple(), 0, 1.0))
                    seen_in_record.add((window, start))
        return hits

    seed_index: dict[tuple[int, str], list[str]] = {}
    seed_ranges = _seed_ranges(query_length, max_mismatches, seed_size)
    for query in query_set:
        for seed_start, seed_end in seed_ranges:
            seed_index.setdefault((seed_start, query[seed_start:seed_end]), []).append(query)

    for record in transcriptome:
        if record.id in target_ids:
            continue
        seen_in_record: set[tuple[str, int]] = set()
        for start in range(0, len(record.sequence) - query_length + 1):
            window = record.sequence[start : start + query_length]
            possible: list[str] = []
            for seed_start, seed_end in seed_ranges:
                possible.extend(seed_index.get((seed_start, window[seed_start:seed_end]), []))
            for query in possible:
                key = (query, start)
                if key in seen_in_record:
                    continue
                mismatch_positions = mismatch_positions_for(query, window, stop_after=max_mismatches)
                if len(mismatch_positions) <= max_mismatches:
                    seed_mismatches = sum(1 for pos in mismatch_positions if 2 <= pos <= 8)
                    hits.append(
                        OffTargetHit(
                            query=query,
                            transcript_id=record.id,
                            position=start + 1,
                            mismatches=len(mismatch_positions),
                            mismatch_positions=tuple(mismatch_positions),
                            seed_mismatches=seed_mismatches,
                            risk_score=offtarget_risk_score(len(mismatch_positions), seed_mismatches),
                        )
                    )
                    seen_in_record.add(key)
    return hits


def _seed_ranges(query_length: int, max_mismatches: int, requested_seed_size: int) -> list[tuple[int, int]]:
    chunks = max_mismatches + 1
    base = query_length // chunks
    ranges: list[tuple[int, int]] = []
    start = 0
    for chunk in range(chunks):
        end = start + base
        if chunk < query_length % chunks:
            end += 1
        ranges.append((start, end))
        start = end
    if all(end - start >= requested_seed_size for start, end in ranges):
        return [(start, start + requested_seed_size) for start, _ in ranges]
    return ranges


def hamming_distance(left: str, right: str, stop_after: int | None = None) -> int:
    if len(left) != len(right):
        raise ValueError("Hamming distance requires strings of equal length")
    distance = 0
    for a, b in zip(left, right):
        if a != b:
            distance += 1
            if stop_after is not None and distance > stop_after:
                return distance
    return distance


def mismatch_positions_for(left: str, right: str, stop_after: int | None = None) -> list[int]:
    if len(left) != len(right):
        raise ValueError("Mismatch positions require strings of equal length")
    positions: list[int] = []
    for index, (a, b) in enumerate(zip(left, right), start=1):
        if a != b:
            positions.append(index)
            if stop_after is not None and len(positions) > stop_after:
                return positions
    return positions


def offtarget_risk_score(mismatches: int, seed_mismatches: int) -> float:
    score = 1.0 - (0.18 * mismatches) - (0.22 * seed_mismatches)
    return max(0.0, score)
