#!/bin/bash
# Probes for the Elsevier VERSION OF RECORD of Kukavica, JDE 194 (2003) 39-50.
# Re-runnable. Writes logs/probe_vor.tsv (route, url, http, bytes, ctype, saved).
set -u
D="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$D/logs/probe_vor.tsv"
RAW="$D/raw2"
mkdir -p "$RAW"
UA="Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/128.0.0.0 Safari/537.36"
printf "route\turl\thttp\tbytes\tcontent_type\tsaved\n" > "$OUT"
probe () { # route url outfile [extra curl args...]
  local route="$1" url="$2" out="$3"; shift 3
  local code ctype bytes
  code=$(curl -sL -m 45 -A "$UA" \
    -H 'Accept: text/html,application/xhtml+xml,application/pdf,*/*' \
    -H 'Accept-Language: en-US,en;q=0.9' \
    "$@" -o "$RAW/$out" -w '%{http_code}\t%{content_type}' "$url" 2>/dev/null) || code="000\t-"
  bytes=$(wc -c < "$RAW/$out" | tr -d ' ')
  printf "%s\t%s\t%s\t%s\t%s\n" "$route" "$url" "$(printf '%b' "$code" | cut -f1)" "$bytes" "$(printf '%b' "$code" | cut -f2)" >> "$OUT"
  printf "%-28s %-8s %-10s %s\n" "$route" "$(printf '%b' "$code" | cut -f1)" "$bytes" "$out"
}
probe sciencedirect-landing "https://www.sciencedirect.com/science/article/pii/S0022039603001530"            sd_landing.html
probe sciencedirect-pdf     "https://www.sciencedirect.com/science/article/pii/S0022039603001530/pdf"         sd_pdf.bin
probe sciencedirect-pdfft   "https://www.sciencedirect.com/science/article/pii/S0022039603001530/pdfft"       sd_pdfft.bin
probe doi-resolve           "https://doi.org/10.1016/S0022-0396(03)00153-0"                                   doi_resolve.html
probe elsevier-api-doi      "https://api.elsevier.com/content/article/doi/10.1016/S0022-0396(03)00153-0"      els_api.xml
probe openaire              "https://api.openaire.eu/search/publications?doi=10.1016/S0022-0396(03)00153-0"   openaire.xml
probe core-search           "https://api.core.ac.uk/v3/search/works?q=%22bounded%20initial%20data%22%20Kukavica" core.json
probe base-search           "https://www.base-search.net/Search/Results?lookfor=Kukavica+bounded+initial+data" base.html
probe fatcat-doi            "https://api.fatcat.wiki/v0/release/lookup?doi=10.1016/s0022-0396(03)00153-0&expand=files" fatcat.json
probe scholar-archive       "https://scholar.archive.org/search?q=%22bounded+initial+data%22+Kukavica"        iascholar.html
probe researchgate          "https://www.researchgate.net/search/publication?q=Kukavica%20bounded%20initial%20data" rg.html
probe s2-openaccess         "https://api.semanticscholar.org/graph/v1/paper/DOI:10.1016/S0022-0396(03)00153-0?fields=title,openAccessPdf,externalIds" s2.json
probe openalex              "https://api.openalex.org/works/doi:10.1016/s0022-0396(03)00153-0"                openalex.json
probe unpaywall             "https://api.unpaywall.org/v2/10.1016/S0022-0396(03)00153-0?email=contact@example.invalid" unpaywall.json
echo "--- log: $OUT"
