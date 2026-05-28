#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [ ! -x "$ROOT/.venv/bin/python" ]; then
  "$ROOT/scripts/bootstrap.sh"
fi

"$ROOT/.venv/bin/python" -m rnai_designer.cli "$@"
