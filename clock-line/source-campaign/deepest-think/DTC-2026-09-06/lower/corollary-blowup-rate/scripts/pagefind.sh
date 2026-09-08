#!/bin/bash
# pagefind.sh <pdf> <string>  -- prints the 1-based PDF page numbers containing <string>
pdf="$1"; pat="$2"
n=$(pdfinfo "$pdf" | awk '/^Pages:/{print $2}')
for ((p=1;p<=n;p++)); do
  if pdftotext -layout -f $p -l $p "$pdf" - 2>/dev/null | grep -qF "$pat"; then echo "$p"; fi
done
