from __future__ import annotations

import argparse
from pathlib import Path

from .accessibility import AccessibilityConfig, compute_accessibility
from .bowtie import run_bowtie_v_mode
from .bowtie2 import build_bowtie2_index, run_bowtie2
from .deps import check_dependencies, format_dependency_statuses
from .design import (
    DesignConfig,
    build_construct_sequence,
    collect_density_points,
    collect_sirna_details,
    design_candidates,
    write_candidates_tsv,
    write_density_tsv,
    write_sirna_details_tsv,
)
from .fasta import FastaRecord, read_fasta, write_fasta


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.check_deps:
        print(format_dependency_statuses(check_dependencies()))
        return 0

    missing = [name for name in ("targets", "transcriptome", "out_prefix") if getattr(args, name) is None]
    if missing:
        parser.error(f"the following arguments are required unless --check-deps is used: {', '.join('--' + item.replace('_', '-') for item in missing)}")

    targets = read_fasta(args.targets)
    transcriptome = read_fasta(args.transcriptome)
    config = DesignConfig(
        min_len=args.min_len,
        max_len=args.max_len,
        step=args.step,
        sirna_size=args.sirna_size,
        max_target_mismatches=args.max_target_mismatches,
        min_alignment_identity=args.min_alignment_identity,
        min_shared_fraction=args.min_shared_fraction,
        min_efficient_sirnas=args.min_efficient_sirnas,
        efficient_sirna_score=args.efficient_sirna_score,
        max_offtarget_transcripts=args.max_offtarget_transcripts,
        max_offtarget_mismatches=args.max_offtarget_mismatches,
        offtarget_seed_size=args.offtarget_seed_size,
        min_gc=args.min_gc,
        max_gc=args.max_gc,
        max_homopolymer=args.max_homopolymer,
        max_candidates=args.max_candidates,
        spacer=args.spacer,
    )

    accessibility_config = AccessibilityConfig(
        method=args.accessibility_method,
        executable=args.rnaplfold_executable,
        window=args.rnaplfold_window,
        span=args.rnaplfold_span,
        unpaired_length=args.sirna_size,
        require_external=args.require_rnaplfold,
    )
    accessibility_by_offset = compute_accessibility(targets[0], accessibility_config)
    bowtie2_index = args.bowtie2_index
    if args.build_bowtie2_index:
        bowtie2_index = args.build_bowtie2_index
        build_bowtie2_index(
            args.transcriptome,
            bowtie2_index,
            executable=args.bowtie2_build_executable,
        )

    candidates = design_candidates(targets, transcriptome, config, accessibility_by_offset=accessibility_by_offset)
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    tsv_path = out_prefix.with_name(out_prefix.name + "_candidates.tsv")
    fasta_path = out_prefix.with_name(out_prefix.name + "_candidates.fa")
    construct_path = out_prefix.with_name(out_prefix.name + "_constructs.fa")
    sirna_path = out_prefix.with_name(out_prefix.name + "_sirnas.tsv")
    density_path = out_prefix.with_name(out_prefix.name + "_density.tsv")

    write_candidates_tsv(candidates, tsv_path, [record.id for record in targets])
    sirna_details = collect_sirna_details(
        candidates,
        targets,
        transcriptome,
        config,
        accessibility_by_offset=accessibility_by_offset,
    )
    density_points = collect_density_points(sirna_details, config.sirna_size, config.efficient_sirna_score)
    write_sirna_details_tsv(sirna_details, sirna_path, [record.id for record in targets])
    write_density_tsv(density_points, density_path)
    write_fasta(
        [
            FastaRecord(
                id=candidate.candidate_id,
                description=(
                    f"{candidate.candidate_id} reference={candidate.reference_id} "
                    f"start={candidate.start} end={candidate.end} score={candidate.score:.4f}"
                ),
                sequence=candidate.sequence,
            )
            for candidate in candidates
        ],
        fasta_path,
    )
    write_fasta(
        [
            FastaRecord(
                id=f"{candidate.candidate_id}_sense_spacer_antisense",
                description=(
                    f"{candidate.candidate_id}_sense_spacer_antisense "
                    f"spacer_length={len(config.spacer)}"
                ),
                sequence=build_construct_sequence(candidate.sequence, config.spacer),
            )
            for candidate in candidates
        ],
        construct_path,
    )
    if args.bowtie_index:
        bowtie_path = out_prefix.with_name(out_prefix.name + "_bowtie.tsv")
        query_by_id = {
            f"{detail.candidate_id}_{detail.offset}": detail.sirna
            for detail in sirna_details
        }
        bowtie_hits = run_bowtie_v_mode(
            query_by_id,
            bowtie_index=args.bowtie_index,
            mismatches=config.max_offtarget_mismatches,
            executable=args.bowtie_executable,
        )
        with bowtie_path.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write("query_id\ttranscript_id\tposition\tstrand\tmismatches\n")
            for hit in bowtie_hits:
                handle.write(
                    f"{hit.query_id}\t{hit.transcript_id}\t{hit.position}\t{hit.strand}\t{hit.mismatches}\n"
                )
    if bowtie2_index:
        bowtie2_path = out_prefix.with_name(out_prefix.name + "_bowtie2.tsv")
        query_by_id = {
            f"{detail.candidate_id}_{detail.offset}": detail.sirna
            for detail in sirna_details
        }
        bowtie2_hits = run_bowtie2(
            query_by_id,
            bowtie2_index=bowtie2_index,
            executable=args.bowtie2_executable,
            seed_size=config.offtarget_seed_size,
        )
        with bowtie2_path.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write("query_id\ttranscript_id\tposition\tstrand\tmismatches\n")
            for hit in bowtie2_hits:
                handle.write(
                    f"{hit.query_id}\t{hit.transcript_id}\t{hit.position}\t{hit.strand}\t{hit.mismatches}\n"
                )

    print(f"Wrote {len(candidates)} candidates")
    print(f"TSV: {tsv_path}")
    print(f"FASTA: {fasta_path}")
    print(f"Construct FASTA: {construct_path}")
    print(f"siRNA TSV: {sirna_path}")
    print(f"Density TSV: {density_path}")
    if args.bowtie_index:
        print(f"Bowtie TSV: {bowtie_path}")
    if bowtie2_index:
        print(f"Bowtie2 TSV: {bowtie2_path}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Design haplotype-aware plant RNAi trigger candidates from target transcripts."
    )
    parser.add_argument("--targets", help="FASTA with one or more target transcripts")
    parser.add_argument("--transcriptome", help="Whole-transcriptome FASTA used for off-target screening")
    parser.add_argument("--out-prefix", help="Output path prefix")
    parser.add_argument("--check-deps", action="store_true", help="Print dependency status and exit")
    parser.add_argument("--min-len", type=int, default=300, help="Minimum candidate fragment length")
    parser.add_argument("--max-len", type=int, default=450, help="Maximum candidate fragment length")
    parser.add_argument("--step", type=int, default=25, help="Sliding-window step size")
    parser.add_argument("--sirna-size", type=int, default=21, help="siRNA k-mer size")
    parser.add_argument(
        "--max-target-mismatches",
        type=int,
        default=0,
        help="Mismatches allowed when counting whether a candidate siRNA is shared by target sequences",
    )
    parser.add_argument(
        "--min-alignment-identity",
        type=float,
        default=0.85,
        help="Minimum pairwise alignment identity for every haplotype-projected candidate window",
    )
    parser.add_argument(
        "--min-shared-fraction",
        type=float,
        default=0.70,
        help="Minimum fraction of candidate siRNAs present in every aligned target haplotype window",
    )
    parser.add_argument(
        "--min-efficient-sirnas",
        type=int,
        default=5,
        help="Minimum number of candidate siRNAs passing the heuristic efficiency threshold",
    )
    parser.add_argument(
        "--efficient-sirna-score",
        type=float,
        default=1.5,
        help="Minimum per-siRNA heuristic efficiency score counted as efficient",
    )
    parser.add_argument(
        "--max-offtarget-transcripts",
        type=int,
        default=0,
        help="Maximum number of non-target transcript IDs hit by candidate siRNAs; use -1 to skip internal scanning",
    )
    parser.add_argument(
        "--max-offtarget-mismatches",
        type=int,
        default=1,
        help="Maximum mismatches allowed when scanning candidate siRNAs against non-target transcripts",
    )
    parser.add_argument(
        "--offtarget-seed-size",
        type=int,
        default=12,
        help="Seed size used to accelerate approximate off-target scanning",
    )
    parser.add_argument("--min-gc", type=float, default=30.0, help="Minimum GC percentage")
    parser.add_argument("--max-gc", type=float, default=60.0, help="Maximum GC percentage")
    parser.add_argument("--max-homopolymer", type=int, default=8, help="Maximum allowed homopolymer run")
    parser.add_argument("--max-candidates", type=int, default=25, help="Maximum ranked candidates to write")
    parser.add_argument(
        "--spacer",
        default="NNNNNNNNNNNNNNNNNNNN",
        help="Spacer or intron placeholder sequence for construct FASTA output",
    )
    parser.add_argument(
        "--bowtie-index",
        help="Optional prebuilt Bowtie index basename for external siRNA hit verification",
    )
    parser.add_argument(
        "--bowtie-executable",
        default="bowtie",
        help="Bowtie executable used when --bowtie-index is supplied",
    )
    parser.add_argument("--bowtie2-index", help="Optional prebuilt Bowtie2 index basename for external verification")
    parser.add_argument(
        "--build-bowtie2-index",
        help="Build a Bowtie2 index from --transcriptome at this basename, then use it for verification",
    )
    parser.add_argument("--bowtie2-executable", default="bowtie2", help="Bowtie2 executable or local aligner path")
    parser.add_argument(
        "--bowtie2-build-executable",
        default="bowtie2-build",
        help="Bowtie2-build executable or local build path",
    )
    parser.add_argument(
        "--accessibility-method",
        choices=["none", "heuristic", "rnaplfold"],
        default="heuristic",
        help="Target accessibility backend; use rnaplfold for si-Fi-style ViennaRNA accessibility",
    )
    parser.add_argument(
        "--rnaplfold-executable",
        default="RNAplfold",
        help="RNAplfold executable used with --accessibility-method rnaplfold",
    )
    parser.add_argument("--rnaplfold-window", type=int, default=80, help="RNAplfold -W local window")
    parser.add_argument("--rnaplfold-span", type=int, default=40, help="RNAplfold -L maximum base-pair span")
    parser.add_argument(
        "--require-rnaplfold",
        action="store_true",
        help="Fail if --accessibility-method rnaplfold is selected but RNAplfold is unavailable",
    )
    return parser


if __name__ == "__main__":
    raise SystemExit(main())
