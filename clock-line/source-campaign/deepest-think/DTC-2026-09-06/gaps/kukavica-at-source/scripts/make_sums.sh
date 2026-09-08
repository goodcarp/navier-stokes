#!/bin/bash
# Recompute SHA256SUMS over the whole seat folder (excluding the manifests themselves).
set -eu
D="$(cd "$(dirname "$0")/.." && pwd)"
cd "$D"
find . -type f ! -name 'SHA256SUMS' ! -name 'SHA256SUMS.round1' ! -name '.DS_Store' -print0 \
  | sort -z | xargs -0 shasum -a 256 > SHA256SUMS
wc -l < SHA256SUMS
