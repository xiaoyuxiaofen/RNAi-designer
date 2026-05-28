from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PairwiseAlignment:
    reference_id: str
    target_id: str
    reference_aligned: str
    target_aligned: str
    identity: float


@dataclass(frozen=True)
class WindowProjection:
    target_id: str
    aligned_identity: float
    ungapped_length: int
    sequence: str


def needleman_wunsch(
    reference_id: str,
    reference: str,
    target_id: str,
    target: str,
    match_score: int = 2,
    mismatch_score: int = -1,
    gap_score: int = -2,
) -> PairwiseAlignment:
    """Global pairwise alignment for haplotype transcripts."""

    rows = len(reference) + 1
    cols = len(target) + 1
    scores = [[0] * cols for _ in range(rows)]
    traceback = [[""] * cols for _ in range(rows)]

    for i in range(1, rows):
        scores[i][0] = scores[i - 1][0] + gap_score
        traceback[i][0] = "up"
    for j in range(1, cols):
        scores[0][j] = scores[0][j - 1] + gap_score
        traceback[0][j] = "left"

    for i in range(1, rows):
        ref_base = reference[i - 1]
        for j in range(1, cols):
            target_base = target[j - 1]
            diagonal = scores[i - 1][j - 1] + (match_score if ref_base == target_base else mismatch_score)
            up = scores[i - 1][j] + gap_score
            left = scores[i][j - 1] + gap_score
            best = max(diagonal, up, left)
            scores[i][j] = best
            if best == diagonal:
                traceback[i][j] = "diag"
            elif best == up:
                traceback[i][j] = "up"
            else:
                traceback[i][j] = "left"

    aligned_ref: list[str] = []
    aligned_target: list[str] = []
    i = len(reference)
    j = len(target)
    while i > 0 or j > 0:
        move = traceback[i][j]
        if move == "diag":
            aligned_ref.append(reference[i - 1])
            aligned_target.append(target[j - 1])
            i -= 1
            j -= 1
        elif move == "up":
            aligned_ref.append(reference[i - 1])
            aligned_target.append("-")
            i -= 1
        else:
            aligned_ref.append("-")
            aligned_target.append(target[j - 1])
            j -= 1

    reference_aligned = "".join(reversed(aligned_ref))
    target_aligned = "".join(reversed(aligned_target))
    return PairwiseAlignment(
        reference_id=reference_id,
        target_id=target_id,
        reference_aligned=reference_aligned,
        target_aligned=target_aligned,
        identity=alignment_identity(reference_aligned, target_aligned),
    )


def alignment_identity(reference_aligned: str, target_aligned: str) -> float:
    comparable = 0
    matches = 0
    for ref_base, target_base in zip(reference_aligned, target_aligned):
        if ref_base == "-" or target_base == "-":
            continue
        comparable += 1
        if ref_base == target_base:
            matches += 1
    if comparable == 0:
        return 0.0
    return matches / comparable


def project_reference_window(
    alignment: PairwiseAlignment,
    reference_start: int,
    reference_end: int,
) -> WindowProjection:
    """Project a 0-based half-open reference window onto an aligned target."""

    ref_pos = 0
    in_window = False
    target_chars: list[str] = []
    comparable = 0
    matches = 0
    for ref_base, target_base in zip(alignment.reference_aligned, alignment.target_aligned):
        if ref_base != "-":
            if ref_pos == reference_start:
                in_window = True
            if ref_pos == reference_end:
                in_window = False
            ref_pos += 1

        if not in_window:
            continue
        if target_base != "-":
            target_chars.append(target_base)
        if ref_base != "-" and target_base != "-":
            comparable += 1
            if ref_base == target_base:
                matches += 1

    identity = matches / comparable if comparable else 0.0
    return WindowProjection(
        target_id=alignment.target_id,
        aligned_identity=identity,
        ungapped_length=len(target_chars),
        sequence="".join(target_chars),
    )
