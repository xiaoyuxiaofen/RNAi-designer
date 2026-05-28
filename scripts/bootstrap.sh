#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON="${PYTHON:-python3}"

if [ ! -x "$ROOT/.venv/bin/python" ]; then
  "$PYTHON" -m venv "$ROOT/.venv"
fi

"$ROOT/.venv/bin/python" -m pip install --upgrade pip
"$ROOT/.venv/bin/python" -m pip install -e "$ROOT"

TOOLS="$ROOT/tools"
mkdir -p "$TOOLS/python-packages"

if [ "${SKIP_DOWNLOADS:-0}" != "1" ]; then
  "$ROOT/.venv/bin/python" -m pip install --target "$TOOLS/python-packages" ViennaRNA==2.7.2
  cat <<'MSG'

Bowtie/Bowtie2:
  Linux/macOS users can install Bowtie2 through conda, brew, apt, or from SourceForge.
  The app can also use tools/bowtie2 if you place a platform-matching Bowtie2 binary there.
MSG
fi

"$ROOT/.venv/bin/python" -m rnai_designer.cli --check-deps
echo
echo "Ready. Use scripts/run.sh to launch RNAi Designer."
