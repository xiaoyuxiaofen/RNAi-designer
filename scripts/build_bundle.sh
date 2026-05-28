#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NAME="${1:-rnai-designer-unix}"
DIST="$ROOT/dist"
STAGE="$DIST/$NAME"
ZIP="$DIST/$NAME.zip"

rm -rf "$STAGE" "$ZIP"
mkdir -p "$STAGE"

for item in src tests examples scripts README.md pyproject.toml .gitignore; do
  if [ -e "$ROOT/$item" ]; then
    cp -R "$ROOT/$item" "$STAGE/"
  fi
done

if [ -d "$ROOT/tools" ]; then
  cp -R "$ROOT/tools" "$STAGE/"
fi

(cd "$DIST" && zip -qr "$ZIP" "$NAME")
echo "$ZIP"
