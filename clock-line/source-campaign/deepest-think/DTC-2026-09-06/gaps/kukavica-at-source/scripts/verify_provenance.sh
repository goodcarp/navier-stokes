#!/bin/bash
# Re-derives the provenance facts for pdf/kukavica-JDE194-2003-a28-authorcopy.pdf.
# 1. our bytes reproduce the Wayback CDX SHA-1 digest for the archived object
# 2. an independent second snapshot of the same object is byte-identical
set -eu
cd "$(dirname "$0")/.."
UA='Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/125.0 Safari/537.36'
F=pdf/kukavica-JDE194-2003-a28-authorcopy.pdf
echo "sha256 local : $(shasum -a 256 "$F" | awk '{print $1}')"
python3 - "$F" <<'PY'
import hashlib,base64,sys
d=open(sys.argv[1],'rb').read()
print("sha1-b32 local:",base64.b32encode(hashlib.sha1(d).digest()).decode().rstrip('='),"bytes:",len(d))
PY
echo "CDX rows for the archived object:"
curl -sS -m 60 -A "$UA" "http://web.archive.org/cdx/search/cdx?url=www-bcf.usc.edu/~kukavica/pdf/a28.pdf&output=text&limit=50"
echo "second snapshot fetch + compare:"
curl -sSL -m 120 -A "$UA" -o /tmp/a28_snap2.pdf "https://web.archive.org/web/20170628022520id_/http://www-bcf.usc.edu/%7Ekukavica/pdf/a28.pdf"
shasum -a 256 "$F" /tmp/a28_snap2.pdf
