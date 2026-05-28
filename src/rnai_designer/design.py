from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from .alignment import needleman_wunsch, project_reference_window
from .fasta import FastaRecord
from .offtarget import OffTargetHit, scan_offtargets
from .sirna import SirnaFeatures, gc_fraction, kmers, sirna_features


@dataclass(frozen=True)
class DesignConfig:
    min_len: int = 300
    max_len: int = 450
    step: int = 25
    sirna_size: int = 21
    max_target_mismatches: int = 0
    min_alignment_identity: float = 0.85
    min_shared_fraction: float = 0.70
    min_efficient_sirnas: int = 5
    efficient_sirna_score: float = 1.5
    max_offtarget_transcripts: int = 0
    max_offtarget_mismatches: int = 1
    offtarget_seed_size: int = 12
    min_gc: float = 30.0
    max_gc: float = 60.0
    max_homopolymer: int = 8
    max_candidates: int = 25
    spacer: str = "NNNNNNNNNNNNNNNNNNNN"


@dataclass(frozen=True)
class TargetProjectionStats:
    target_id: str
    alignment_identity: float
    projected_length: int
    shared_sirnas: int
    shared_fraction: float


@dataclass(frozen=True)
class Candidate:
    rank: int
    candidate_id: str
    reference_id: str
    start: int
    end: int
    length: int
    gc_percent: float
    total_sirnas: int
    efficient_sirnas: int
    mean_sirna_score: float
    mean_accessibility: float
    min_accessibility: float
    min_alignment_identity: float
    min_shared_count: int
    min_shared_fraction: float
    target_stats: dict[str, TargetProjectionStats]
    offtarget_transcript_count: int
    offtarget_sirna_count: int
    offtarget_ids: tuple[str, ...]
    sequence: str
    score: float


@dataclass(frozen=True)
class SirnaDetail:
    candidate_id: str
    sirna: str
    offset: int
    features: SirnaFeatures
    accessibility: float
    target_presence: dict[str, bool]
    offtarget_hits: tuple[OffTargetHit, ...]


@dataclass(frozen=True)
class DensityPoint:
    candidate_id: str
    reference_position: int
    total_sirnas_covering: int
    efficient_sirnas_covering: int
    offtarget_sirnas_covering: int


def design_candidates(
    targets: Sequence[FastaRecord],
    transcriptome: Sequence[FastaRecord],
    config: DesignConfig,
    accessibility_by_offset: dict[int, float] | None = None,
) -> list[Candidate]:
    _validate_config(config)
    if len(targets) < 1:
        raise ValueError("--targets must contain at least one target sequence")

    reference = targets[0]
    target_ids = {record.id for record in targets}
    alignments = [
        needleman_wunsch(reference.id, reference.sequence, target.id, target.sequence)
        for target in targets[1:]
    ]

    raw_candidates: list[Candidate] = []
    for length in range(config.min_len, config.max_len + 1):
        for start0 in range(0, len(reference.sequence) - length + 1, config.step):
            sequence = reference.sequence[start0 : start0 + length]
            if not _passes_sequence_filters(sequence, config):
                continue

            candidate_sirnas = sorted(set(kmers(sequence, config.sirna_size)))
            if not candidate_sirnas:
                continue
            feature_by_sirna = {sirna: sirna_features(sirna) for sirna in candidate_sirnas}
            sirna_scores = [features.efficiency_score for features in feature_by_sirna.values()]
            accessibility_values = _window_accessibilities(start0, sequence, config.sirna_size, accessibility_by_offset)
            mean_accessibility = _mean(accessibility_values)
            min_accessibility = min(accessibility_values) if accessibility_values else 1.0
            efficient_sirnas = sum(score >= config.efficient_sirna_score for score in sirna_scores)
            if efficient_sirnas < config.min_efficient_sirnas:
                continue

            target_stats: dict[str, TargetProjectionStats] = {}
            reference_stats = TargetProjectionStats(
                target_id=reference.id,
                alignment_identity=1.0,
                projected_length=length,
                shared_sirnas=len(candidate_sirnas),
                shared_fraction=1.0,
            )
            target_stats[reference.id] = reference_stats

            for alignment in alignments:
                projection = project_reference_window(alignment, start0, start0 + length)
                shared = _count_matching_sirnas(
                    candidate_sirnas,
                    projection.sequence,
                    config.sirna_size,
                    config.max_target_mismatches,
                )
                shared_fraction = shared / len(candidate_sirnas)
                target_stats[projection.target_id] = TargetProjectionStats(
                    target_id=projection.target_id,
                    alignment_identity=projection.aligned_identity,
                    projected_length=projection.ungapped_length,
                    shared_sirnas=shared,
                    shared_fraction=shared_fraction,
                )

            min_alignment_identity = min(stats.alignment_identity for stats in target_stats.values())
            min_shared_count = min(stats.shared_sirnas for stats in target_stats.values())
            min_shared_fraction = min(stats.shared_fraction for stats in target_stats.values())
            if min_alignment_identity < config.min_alignment_identity:
                continue
            if min_shared_fraction < config.min_shared_fraction:
                continue

            offtarget_hits: list[OffTargetHit] = []
            if config.max_offtarget_transcripts >= 0:
                offtarget_hits = scan_offtargets(
                    candidate_sirnas,
                    transcriptome,
                    target_ids=target_ids,
                    max_mismatches=config.max_offtarget_mismatches,
                    seed_size=config.offtarget_seed_size,
                )
            offtarget_ids = tuple(sorted({hit.transcript_id for hit in offtarget_hits}))
            if config.max_offtarget_transcripts >= 0 and len(offtarget_ids) > config.max_offtarget_transcripts:
                continue

            mean_sirna_score = sum(sirna_scores) / len(sirna_scores)
            score = _score_candidate(
                min_alignment_identity=min_alignment_identity,
                min_shared_fraction=min_shared_fraction,
                efficient_sirnas=efficient_sirnas,
                total_sirnas=len(candidate_sirnas),
                mean_sirna_score=mean_sirna_score,
                mean_accessibility=mean_accessibility,
                offtarget_transcript_count=len(offtarget_ids),
                offtarget_sirna_count=len(offtarget_hits),
                sequence=sequence,
            )
            raw_candidates.append(
                Candidate(
                    rank=0,
                    candidate_id="",
                    reference_id=reference.id,
                    start=start0 + 1,
                    end=start0 + length,
                    length=length,
                    gc_percent=_gc_percent(sequence),
                    total_sirnas=len(candidate_sirnas),
                    efficient_sirnas=efficient_sirnas,
                    mean_sirna_score=mean_sirna_score,
                    mean_accessibility=mean_accessibility,
                    min_accessibility=min_accessibility,
                    min_alignment_identity=min_alignment_identity,
                    min_shared_count=min_shared_count,
                    min_shared_fraction=min_shared_fraction,
                    target_stats=target_stats,
                    offtarget_transcript_count=len(offtarget_ids),
                    offtarget_sirna_count=len(offtarget_hits),
                    offtarget_ids=offtarget_ids,
                    sequence=sequence,
                    score=score,
                )
            )

    ranked = sorted(
        raw_candidates,
        key=lambda item: (
            -item.score,
            item.offtarget_transcript_count,
            -item.efficient_sirnas,
            -item.min_alignment_identity,
            -item.length,
            item.start,
        ),
    )
    candidates: list[Candidate] = []
    for rank, candidate in enumerate(ranked[: config.max_candidates], start=1):
        candidate_id = f"{reference.id}_RNAi_{rank:02d}_{candidate.start}_{candidate.end}"
        candidates.append(_with_rank_and_id(candidate, rank, candidate_id))
    return candidates


def collect_sirna_details(
    candidates: Sequence[Candidate],
    targets: Sequence[FastaRecord],
    transcriptome: Sequence[FastaRecord],
    config: DesignConfig,
    accessibility_by_offset: dict[int, float] | None = None,
) -> list[SirnaDetail]:
    target_sequences = {record.id: record.sequence for record in targets}
    target_ids = {record.id for record in targets}
    details: list[SirnaDetail] = []
    for candidate in candidates:
        candidate_kmers = sorted(set(kmers(candidate.sequence, config.sirna_size)))
        hits: list[OffTargetHit] = []
        if config.max_offtarget_transcripts >= 0:
            hits = scan_offtargets(
                candidate_kmers,
                transcriptome,
                target_ids=target_ids,
                max_mismatches=config.max_offtarget_mismatches,
                seed_size=config.offtarget_seed_size,
            )
        hits_by_sirna: dict[str, list[OffTargetHit]] = {}
        for hit in hits:
            hits_by_sirna.setdefault(hit.query, []).append(hit)
        for offset, sirna in _sirnas_with_offsets(candidate.sequence, config.sirna_size):
            reference_offset = candidate.start + offset - 1
            details.append(
                SirnaDetail(
                    candidate_id=candidate.candidate_id,
                    sirna=sirna,
                    offset=offset,
                    features=sirna_features(sirna),
                    accessibility=_accessibility_at(reference_offset, accessibility_by_offset),
                    target_presence={
                        target_id: _sirna_matches_sequence(
                            sirna,
                            sequence,
                            config.max_target_mismatches,
                        )
                        for target_id, sequence in target_sequences.items()
                    },
                    offtarget_hits=tuple(hits_by_sirna.get(sirna, [])),
                )
            )
    return details


def collect_density_points(
    details: Sequence[SirnaDetail],
    sirna_size: int,
    efficient_sirna_score: float,
) -> list[DensityPoint]:
    by_candidate: dict[str, list[SirnaDetail]] = {}
    for detail in details:
        by_candidate.setdefault(detail.candidate_id, []).append(detail)

    points: list[DensityPoint] = []
    for candidate_id, candidate_details in by_candidate.items():
        if not candidate_details:
            continue
        max_position = max(detail.offset + sirna_size - 1 for detail in candidate_details)
        for position in range(1, max_position + 1):
            covering = [
                detail
                for detail in candidate_details
                if detail.offset <= position <= detail.offset + sirna_size - 1
            ]
            points.append(
                DensityPoint(
                    candidate_id=candidate_id,
                    reference_position=position,
                    total_sirnas_covering=len(covering),
                    efficient_sirnas_covering=sum(
                        detail.features.efficiency_score >= efficient_sirna_score for detail in covering
                    ),
                    offtarget_sirnas_covering=sum(bool(detail.offtarget_hits) for detail in covering),
                )
            )
    return points


def write_candidates_tsv(candidates: Sequence[Candidate], path: str | Path, target_ids: Sequence[str]) -> None:
    columns = [
        "rank",
        "candidate_id",
        "reference_id",
        "start",
        "end",
        "length",
        "gc_percent",
        "total_sirnas",
        "efficient_sirnas",
        "mean_sirna_score",
        "mean_accessibility",
        "min_accessibility",
        "min_alignment_identity",
        "min_shared_count",
        "min_shared_fraction",
        *[f"alignment_identity_{target_id}" for target_id in target_ids],
        *[f"projected_length_{target_id}" for target_id in target_ids],
        *[f"shared_sirnas_{target_id}" for target_id in target_ids],
        *[f"shared_fraction_{target_id}" for target_id in target_ids],
        "offtarget_transcript_count",
        "offtarget_sirna_count",
        "offtarget_ids",
        "score",
        "sequence",
    ]
    with Path(path).open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(columns) + "\n")
        for candidate in candidates:
            row = [
                str(candidate.rank),
                candidate.candidate_id,
                candidate.reference_id,
                str(candidate.start),
                str(candidate.end),
                str(candidate.length),
                f"{candidate.gc_percent:.2f}",
                str(candidate.total_sirnas),
                str(candidate.efficient_sirnas),
                f"{candidate.mean_sirna_score:.3f}",
                f"{candidate.mean_accessibility:.3f}",
                f"{candidate.min_accessibility:.3f}",
                f"{candidate.min_alignment_identity:.3f}",
                str(candidate.min_shared_count),
                f"{candidate.min_shared_fraction:.3f}",
                *[_format_target_stat(candidate, target_id, "alignment_identity") for target_id in target_ids],
                *[_format_target_stat(candidate, target_id, "projected_length") for target_id in target_ids],
                *[_format_target_stat(candidate, target_id, "shared_sirnas") for target_id in target_ids],
                *[_format_target_stat(candidate, target_id, "shared_fraction") for target_id in target_ids],
                str(candidate.offtarget_transcript_count),
                str(candidate.offtarget_sirna_count),
                ",".join(candidate.offtarget_ids),
                f"{candidate.score:.4f}",
                candidate.sequence,
            ]
            handle.write("\t".join(row) + "\n")


def write_sirna_details_tsv(details: Sequence[SirnaDetail], path: str | Path, target_ids: Sequence[str]) -> None:
    columns = [
        "candidate_id",
        "offset",
        "sirna",
        "seed_2_8",
        "gc_percent",
        "efficiency_score",
        "accessibility",
        "antisense_5p_stability",
        "passenger_5p_stability",
        "asymmetry",
        "au_5p_7",
        "guide_preferred",
        *[f"present_in_{target_id}" for target_id in target_ids],
        "offtarget_hit_count",
        "offtarget_transcript_count",
        "max_offtarget_risk",
        "offtarget_hits",
    ]
    with Path(path).open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(columns) + "\n")
        for detail in details:
            offtarget_ids = {hit.transcript_id for hit in detail.offtarget_hits}
            max_risk = max((hit.risk_score for hit in detail.offtarget_hits), default=0.0)
            hit_text = ";".join(
                (
                    f"{hit.transcript_id}:{hit.position}:mm={hit.mismatches}:"
                    f"seed_mm={hit.seed_mismatches}:risk={hit.risk_score:.2f}"
                )
                for hit in detail.offtarget_hits
            )
            row = [
                detail.candidate_id,
                str(detail.offset),
                detail.sirna,
                detail.features.seed,
                f"{100.0 * detail.features.gc_fraction:.2f}",
                f"{detail.features.efficiency_score:.3f}",
                f"{detail.accessibility:.3f}",
                f"{detail.features.antisense_5p_stability:.2f}",
                f"{detail.features.passenger_5p_stability:.2f}",
                f"{detail.features.asymmetry:.2f}",
                str(detail.features.au_5p_7),
                str(detail.features.guide_preferred),
                *[str(detail.target_presence.get(target_id, False)) for target_id in target_ids],
                str(len(detail.offtarget_hits)),
                str(len(offtarget_ids)),
                f"{max_risk:.3f}",
                hit_text,
            ]
            handle.write("\t".join(row) + "\n")


def write_density_tsv(points: Sequence[DensityPoint], path: str | Path) -> None:
    columns = [
        "candidate_id",
        "reference_position",
        "total_sirnas_covering",
        "efficient_sirnas_covering",
        "offtarget_sirnas_covering",
    ]
    with Path(path).open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(columns) + "\n")
        for point in points:
            handle.write(
                "\t".join(
                    [
                        point.candidate_id,
                        str(point.reference_position),
                        str(point.total_sirnas_covering),
                        str(point.efficient_sirnas_covering),
                        str(point.offtarget_sirnas_covering),
                    ]
                )
                + "\n"
            )


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return sequence.translate(table)[::-1].upper()


def build_construct_sequence(sequence: str, spacer: str) -> str:
    return f"{sequence.upper()}{spacer.upper()}{reverse_complement(sequence)}"


def _with_rank_and_id(candidate: Candidate, rank: int, candidate_id: str) -> Candidate:
    return Candidate(
        rank=rank,
        candidate_id=candidate_id,
        reference_id=candidate.reference_id,
        start=candidate.start,
        end=candidate.end,
        length=candidate.length,
        gc_percent=candidate.gc_percent,
        total_sirnas=candidate.total_sirnas,
        efficient_sirnas=candidate.efficient_sirnas,
        mean_sirna_score=candidate.mean_sirna_score,
        mean_accessibility=candidate.mean_accessibility,
        min_accessibility=candidate.min_accessibility,
        min_alignment_identity=candidate.min_alignment_identity,
        min_shared_count=candidate.min_shared_count,
        min_shared_fraction=candidate.min_shared_fraction,
        target_stats=candidate.target_stats,
        offtarget_transcript_count=candidate.offtarget_transcript_count,
        offtarget_sirna_count=candidate.offtarget_sirna_count,
        offtarget_ids=candidate.offtarget_ids,
        sequence=candidate.sequence,
        score=candidate.score,
    )


def _sirnas_with_offsets(sequence: str, k: int) -> list[tuple[int, str]]:
    return [
        (start + 1, sequence[start : start + k])
        for start in range(0, len(sequence) - k + 1)
        if "N" not in sequence[start : start + k]
    ]


def _count_matching_sirnas(candidate_sirnas: Sequence[str], target_sequence: str, k: int, max_mismatches: int) -> int:
    target_sirnas = kmers(target_sequence, k)
    return sum(
        1
        for candidate_sirna in candidate_sirnas
        if any(_hamming_at_most(candidate_sirna, target_sirna, max_mismatches) for target_sirna in target_sirnas)
    )


def _sirna_matches_sequence(sirna: str, target_sequence: str, max_mismatches: int) -> bool:
    return any(_hamming_at_most(sirna, target_sirna, max_mismatches) for target_sirna in kmers(target_sequence, len(sirna)))


def _hamming_at_most(left: str, right: str, max_mismatches: int) -> bool:
    mismatches = 0
    for left_base, right_base in zip(left, right):
        if left_base != right_base:
            mismatches += 1
            if mismatches > max_mismatches:
                return False
    return True


def _window_accessibilities(
    reference_start0: int,
    sequence: str,
    sirna_size: int,
    accessibility_by_offset: dict[int, float] | None,
) -> list[float]:
    if not accessibility_by_offset:
        return [1.0 for _ in range(0, len(sequence) - sirna_size + 1)]
    return [
        _accessibility_at(reference_start0 + offset + 1, accessibility_by_offset)
        for offset in range(0, len(sequence) - sirna_size + 1)
    ]


def _accessibility_at(offset: int, accessibility_by_offset: dict[int, float] | None) -> float:
    if not accessibility_by_offset:
        return 1.0
    return accessibility_by_offset.get(offset, 0.0)


def _mean(values: Sequence[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def _format_target_stat(candidate: Candidate, target_id: str, field: str) -> str:
    stats = candidate.target_stats.get(target_id)
    if stats is None:
        return ""
    value = getattr(stats, field)
    if isinstance(value, float):
        return f"{value:.3f}"
    return str(value)


def _passes_sequence_filters(sequence: str, config: DesignConfig) -> bool:
    gc = _gc_percent(sequence)
    if gc < config.min_gc or gc > config.max_gc:
        return False
    return _longest_homopolymer(sequence) <= config.max_homopolymer


def _gc_percent(sequence: str) -> float:
    return 100.0 * gc_fraction(sequence)


def _longest_homopolymer(sequence: str) -> int:
    longest = 0
    current = 0
    previous = ""
    for base in sequence:
        if base == previous:
            current += 1
        else:
            previous = base
            current = 1
        longest = max(longest, current)
    return longest


def _score_candidate(
    min_alignment_identity: float,
    min_shared_fraction: float,
    efficient_sirnas: int,
    total_sirnas: int,
    mean_sirna_score: float,
    mean_accessibility: float,
    offtarget_transcript_count: int,
    offtarget_sirna_count: int,
    sequence: str,
) -> float:
    efficient_fraction = efficient_sirnas / total_sirnas if total_sirnas else 0.0
    gc_distance = abs(_gc_percent(sequence) - 45.0) / 45.0
    return (
        0.40 * min_alignment_identity
        + 0.25 * min_shared_fraction
        + 0.20 * efficient_fraction
        + 0.05 * mean_sirna_score
        + 0.10 * mean_accessibility
        - 0.35 * offtarget_transcript_count
        - 0.015 * offtarget_sirna_count
        - 0.05 * gc_distance
    )


def _validate_config(config: DesignConfig) -> None:
    if config.min_len < config.sirna_size:
        raise ValueError("--min-len must be at least --sirna-size")
    if config.max_len < config.min_len:
        raise ValueError("--max-len must be greater than or equal to --min-len")
    if config.step < 1:
        raise ValueError("--step must be at least 1")
    if not 0 <= config.min_alignment_identity <= 1:
        raise ValueError("--min-alignment-identity must be between 0 and 1")
    if config.max_target_mismatches < 0:
        raise ValueError("--max-target-mismatches must be non-negative")
    if not 0 <= config.min_shared_fraction <= 1:
        raise ValueError("--min-shared-fraction must be between 0 and 1")
    if config.max_offtarget_transcripts < -1:
        raise ValueError("--max-offtarget-transcripts must be -1 or greater")
    if config.max_offtarget_mismatches < 0:
        raise ValueError("--max-offtarget-mismatches must be non-negative")
    if config.offtarget_seed_size < 1 or config.offtarget_seed_size > config.sirna_size:
        raise ValueError("--offtarget-seed-size must be between 1 and --sirna-size")
