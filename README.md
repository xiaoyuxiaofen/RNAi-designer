# RNAi Designer

[中文文档 / Chinese Guide](docs/中文使用说明.md)

Haplotype-aware RNAi trigger candidate designer for plant ihpRNA / hpRNA workflows.

This tool takes target transcript sequences from two or more haplotypes plus a whole-transcriptome FASTA database, then reports 200-500 bp candidate RNAi trigger regions that:

- sit inside aligned, conserved haplotype regions
- produce many candidate 21 nt siRNAs
- have few or no approximate siRNA-like hits in non-target transcripts
- pass GC, low-complexity, and siRNA efficiency filters
- can be exported as candidate FASTA and a ranked TSV table

It is meant to cover the "one-click candidate fragment extraction" step with a si-Fi-like local workflow: haplotype sequence alignment, siRNA generation, off-target scanning, and ranked long RNAi trigger output.

For strict si-Fi-style runs, install Bowtie and ViennaRNA, build a Bowtie transcriptome index, then use:

```powershell
$env:PYTHONPATH = "src"
python -m rnai_designer.cli `
  --targets targets.fa `
  --transcriptome all_transcripts.fa `
  --out-prefix results/PagLBD4 `
  --bowtie-index indexes/all_transcripts `
  --accessibility-method rnaplfold `
  --require-rnaplfold
```

On Windows, this repository can also use local bundled tools under `tools/`:

- `tools/python-packages/ViennaRNA`: ViennaRNA Python bindings used for RNAplfold-equivalent unpaired probabilities.
- `tools/bowtie2/bowtie2-2.5.0-mingw-x86_64`: Bowtie2 external verification backend.
- `tools/bowtie/bowtie-1.2` and `tools/bowtie-legacy/bowtie-1.2-legacy`: Bowtie 1 downloads are kept for compatibility, but these old Windows binaries may require missing legacy DLLs on modern Windows.

The fully local Windows command verified in this workspace is:

```powershell
$env:PYTHONPATH = "src"
python -m rnai_designer.cli `
  --targets examples/targets.fa `
  --transcriptome examples/transcriptome.fa `
  --out-prefix results/PagLBD4 `
  --accessibility-method rnaplfold `
  --require-rnaplfold `
  --build-bowtie2-index results/bt2/transcriptome
```

## Quick Start

```powershell
$env:PYTHONPATH = "src"
python -m rnai_designer.cli `
  --targets examples/targets.fa `
  --transcriptome examples/transcriptome.fa `
  --out-prefix results/PagLBD4 `
  --min-len 300 `
  --max-len 450 `
  --step 25
```

If `python` is not available in the current Codex desktop shell, use the bundled interpreter path shown by the workspace dependency tool.

## One-Click Setup

Windows:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\bootstrap.ps1
scripts\run.ps1 --check-deps
```

Linux/macOS:

```bash
bash scripts/bootstrap.sh
bash scripts/run.sh --check-deps
```

Build a distributable bundle:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\build_bundle.ps1
```

The Windows bundle can include local tools under `tools/`, including ViennaRNA Python bindings and Bowtie2 MinGW binaries. Platform-specific native binaries are not interchangeable across Windows, Linux, and macOS, so build one bundle per operating system when offline deployment is required.

## Run Tests

```powershell
$env:PYTHONPATH = "src"
python -m unittest discover -s tests -v
```

Outputs:

- `results/PagLBD4_candidates.tsv`: ranked candidate table
- `results/PagLBD4_candidates.fa`: candidate stem sequences
- `results/PagLBD4_constructs.fa`: sense-spacer-antisense sequences using the configured spacer
- `results/PagLBD4_sirnas.tsv`: per-siRNA feature and off-target table
- `results/PagLBD4_density.tsv`: per-position siRNA density table

## Input Files

`--targets` should contain the target haplotype transcripts:

```fasta
>Hap1_LBD4
ATG...
>Hap2_LBD4
ATG...
```

`--transcriptome` should contain the whole transcriptome, including the target transcripts and possible paralogs:

```fasta
>Hap1_LBD4
ATG...
>Hap2_LBD4
ATG...
>Other_LBD
ATG...
```

Target records are excluded from off-target counting by FASTA ID.

## Main Scoring Logic

The first target record is treated as the reference haplotype. Every other target record is globally aligned to it with an internal Needleman-Wunsch aligner. For each reference window, the tool projects the same aligned region onto every target haplotype and keeps only windows that are conserved enough.

For each retained aligned window, the tool creates all unique 21 nt siRNA candidates, applies a heuristic efficiency score, and scans the whole transcriptome for non-target hits. Off-target scanning supports bounded mismatches, so this is stricter than exact k-mer screening and closer to the logic used by si-Fi/Bowtie-style workflows.

The internal off-target scanner reports mismatch positions, seed-region mismatches, and a risk score for every hit. The siRNA feature table includes GC, seed sequence, terminal stability, strand asymmetry, A/U enrichment at the guide 5' end, and whether the guide strand is predicted to be preferred.

Target-site accessibility is available through ViennaRNA `RNAplfold` or local ViennaRNA Python bindings. Use `--accessibility-method rnaplfold` for RNAplfold-equivalent probabilities, or `--require-rnaplfold` to fail rather than falling back to the local heuristic. If Bowtie is installed and you already have a transcriptome index, add `--bowtie-index INDEX_BASENAME` to export an external Bowtie v-mode hit table for the selected candidates. If Bowtie 1 is unavailable on Windows, use `--build-bowtie2-index` or `--bowtie2-index` for local external alignment verification.

Default settings are conservative but adjustable:

```text
window length: 300-450 bp
siRNA size: 21 nt
minimum alignment identity in every haplotype window: 0.85
minimum shared target 21-mer fraction: 0.70
minimum efficient siRNAs: 5
maximum off-target mismatches: 1
maximum off-target transcripts: 0
target accessibility: heuristic by default, RNAplfold optional
GC range: 30-60%
```

The TSV output includes per-haplotype alignment identity, projected length, shared siRNA counts, shared siRNA fraction, efficient siRNA count, approximate off-target transcript count, and final score.

## Useful Options

```text
--min-len / --max-len       Candidate fragment size range
--step                      Sliding window step size
--sirna-size                siRNA k-mer size, default 21
--max-target-mismatches     Allowed mismatches when counting siRNAs shared by target sequences
--min-alignment-identity    Required aligned identity in every haplotype window
--min-shared-fraction       Required fraction of window siRNAs present in each target
--min-efficient-sirnas      Required number of high-scoring siRNAs
--efficient-sirna-score     Per-siRNA threshold counted as efficient
--max-offtarget-transcripts Allowed number of non-target transcript IDs hit; -1 skips internal scanning
--max-offtarget-mismatches  Allowed mismatches in non-target siRNA scans
--offtarget-seed-size       Seed size for approximate off-target scanning
--max-candidates            Number of ranked candidates to write
--spacer                    Spacer/intron placeholder used in construct FASTA
--bowtie-index              Optional Bowtie index basename for external verification
--bowtie-executable         Bowtie executable name or path
--bowtie2-index             Optional Bowtie2 index basename for external verification
--build-bowtie2-index       Build and use a Bowtie2 index from --transcriptome
--bowtie2-executable        Bowtie2 executable name or path
--bowtie2-build-executable  Bowtie2-build executable name or path
--accessibility-method      none, heuristic, or rnaplfold
--rnaplfold-executable      RNAplfold executable name or path
--rnaplfold-window          RNAplfold -W local window
--rnaplfold-span            RNAplfold -L maximum base-pair span
--require-rnaplfold         Fail if RNAplfold is unavailable
```
