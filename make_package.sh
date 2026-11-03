#!/bin/bash
# 渡す zip を作る (261027)。開発の記録 (README.md・log-format.md・classifier/README.md)・開発用のファイルは入れない。
# 使い方: bash make_package.sh <出力先のフォルダ>
set -euo pipefail
SRC="$(cd "$(dirname "$0")" && pwd)"
OUT="${1:?出力先のフォルダを指定してください}"
# 相対パスだと下の cd で一時フォルダ内に zip が書かれ、後始末で消えるため絶対パスにする (261029)
OUT="$(cd "$OUT" && pwd)"
VER="$(sed -n '1s/.* \([0-9][0-9.]*\) .*/\1/p' "$SRC/VERSION")"
NAME="16S_pipeline_${VER}"
# 作業コピーは終了時 (成功・失敗・中断とも) に必ず消す。以前は消さずに 1 回 2.9 GiB 残っていた (261029)
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
STAGE="$WORK/$NAME"
mkdir -p "$STAGE"
rsync -a \
  --exclude '.DS_Store' --exclude '._*' --exclude '__pycache__' --exclude '*.pyc' \
  --exclude 'README.md' --exclude 'log-format.md' --exclude '.gitignore' --exclude 'make_package.sh' \
  "$SRC/" "$STAGE/"
[ -f "$SRC/classifier/README.md" ] && rm -f "$STAGE/classifier/README.md"
(cd "$(dirname "$STAGE")" && zip -r -q -X "$OUT/$NAME.zip" "$NAME")
echo "$OUT/$NAME.zip"
(cd "$(dirname "$STAGE")" && find "$NAME" -type f | sort)
