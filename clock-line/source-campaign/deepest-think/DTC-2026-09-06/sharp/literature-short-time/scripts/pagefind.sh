#!/bin/bash
# find the PDF page index containing a given string
pdf="$1"; shift
n=$(pdfinfo "$pdf" | awk '/^Pages:/{print $2}')
for ((p=1;p<=n;p++)); do
  t=$(pdftotext -layout -f $p -l $p "$pdf" - 2>/dev/null)
  for pat in "$@"; do
    if grep -qF "$pat" <<<"$t"; then echo "p$p :: $pat"; fi
  done
done
