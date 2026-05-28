from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SirnaFeatures:
    sequence: str
    gc_fraction: float
    antisense_5p_stability: float
    passenger_5p_stability: float
    asymmetry: float
    au_5p_7: int
    seed: str
    efficiency_score: float
    guide_preferred: bool


def kmers(sequence: str, k: int) -> list[str]:
    return [
        sequence[start : start + k]
        for start in range(0, len(sequence) - k + 1)
        if "N" not in sequence[start : start + k]
    ]


def sirna_efficiency_score(sequence: str) -> float:
    return sirna_features(sequence).efficiency_score


def sirna_features(sequence: str) -> SirnaFeatures:
    """Heuristic plant RNAi siRNA features inspired by si-Fi-like criteria."""

    if not sequence:
        return SirnaFeatures("", 0.0, 0.0, 0.0, 0.0, 0, "", 0.0, False)
    score = 0.0
    gc = gc_fraction(sequence)
    if 0.30 <= gc <= 0.52:
        score += 1.0
    elif 0.25 <= gc <= 0.60:
        score += 0.5

    # Favor lower 5' antisense stability and useful A/U positions.
    if sequence[0] in {"A", "T"}:
        score += 0.5
    if sequence[-1] in {"G", "C"}:
        score += 0.5
    antisense_5p = terminal_stability(sequence[:4])
    passenger_5p = terminal_stability(reverse_complement(sequence)[:4])
    asymmetry = passenger_5p - antisense_5p
    if asymmetry >= 1.0:
        score += 0.75
    elif asymmetry > 0:
        score += 0.35
    au_5p_7 = sum(1 for base in sequence[:7] if base in {"A", "T"})
    if au_5p_7 >= 4:
        score += 0.5
    for position in (2, 9, 12, 18):
        if position < len(sequence) and sequence[position] in {"A", "T"}:
            score += 0.25
    if "GGGG" in sequence or "CCCC" in sequence:
        score -= 0.5
    return SirnaFeatures(
        sequence=sequence,
        gc_fraction=gc,
        antisense_5p_stability=antisense_5p,
        passenger_5p_stability=passenger_5p,
        asymmetry=asymmetry,
        au_5p_7=au_5p_7,
        seed=sequence[1:8],
        efficiency_score=max(0.0, score),
        guide_preferred=asymmetry > 0,
    )


def gc_fraction(sequence: str) -> float:
    return (sequence.count("G") + sequence.count("C")) / len(sequence) if sequence else 0.0


def terminal_stability(sequence: str) -> float:
    """Approximate duplex end stability; lower values mean weaker A/U-rich ends."""

    return sum(2.0 if base in {"A", "T"} else 3.0 for base in sequence.upper())


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return sequence.translate(table)[::-1].upper()
