#!/bin/bash
cd "$(dirname "$0")/.." || exit 1
find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 shasum -a 256 > SHA256SUMS
wc -l < SHA256SUMS
