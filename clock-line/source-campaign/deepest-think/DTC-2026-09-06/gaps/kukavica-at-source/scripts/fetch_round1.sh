#!/bin/bash
# Round 1 acquisition probes for Kukavica, JDE 194 (2003) 39-50, DOI 10.1016/S0022-0396(03)00153-0
# Every route logged: URL, HTTP code, bytes, content-type. No Sci-Hub.
set -u
UA='Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/125.0 Safari/537.36'
LOG=log/round1.tsv
echo -e "route\turl\thttp\tbytes\tcontent_type\tsaved_as" > "$LOG"
try () { # name url outfile [extra curl args...]
  n="$1"; u="$2"; o="$3"; shift 3
  r=$(curl -sSL -m 45 -A "$UA" "$@" -o "raw/$o" -w "%{http_code}\t%{size_download}\t%{content_type}" "$u" 2>>log/round1.err || echo -e "ERR\t0\t-")
  echo -e "$n\t$u\t$r\traw/$o" >> "$LOG"
}
try homepage_dornsife "https://dornsife.usc.edu/profile/igor-kukavica/" hp_dornsife.html
try homepage_bcf "http://www-bcf.usc.edu/~kukavica/" hp_bcf.html
try homepage_dreamhost "https://kukavica.dreamhosters.com/" hp_dreamhost.html
try homepage_dreamhost_pubs "https://kukavica.dreamhosters.com/publications.html" hp_dreamhost_pubs.html
try wayback_cdx_dreamhost "http://web.archive.org/cdx/search/cdx?url=kukavica.dreamhosters.com*&output=text&limit=400&collapse=urlkey&filter=statuscode:200" wb_cdx_dreamhost.txt
try wayback_cdx_bcf "http://web.archive.org/cdx/search/cdx?url=www-bcf.usc.edu/~kukavica*&output=text&limit=400&collapse=urlkey" wb_cdx_bcf.txt
try wayback_cdx_math_usc "http://web.archive.org/cdx/search/cdx?url=math.usc.edu/~kukavica*&output=text&limit=400&collapse=urlkey" wb_cdx_mathusc.txt
try fatcat_lookup "https://api.fatcat.wiki/v0/release/lookup?doi=10.1016/s0022-0396(03)00153-0&expand=files" fatcat.json
try openaire "https://api.openaire.eu/search/publications?doi=10.1016/S0022-0396(03)00153-0" openaire.xml
try s2 "https://api.semanticscholar.org/graph/v1/paper/DOI:10.1016/S0022-0396(03)00153-0?fields=title,externalIds,openAccessPdf,isOpenAccess,abstract" s2.json
try core_search "https://core.ac.uk/search?q=%22local%20uniqueness%20of%20weak%20solutions%20of%20the%20Navier-Stokes%20system%20with%20bounded%20initial%20data%22" core_search.html
try base_search "https://www.base-search.net/Search/Results?lookfor=Kukavica+local+uniqueness+weak+solutions+Navier-Stokes+bounded+initial+data" base_search.html
try sciencedirect_landing "https://www.sciencedirect.com/science/article/pii/S0022039603001530" sd_landing.html
try sciencedirect_pdf "https://www.sciencedirect.com/science/article/pii/S0022039603001530/pdfft?isDTMRedir=true&download=true" sd_pdf.pdf
try elsevier_tdm_plain "https://api.elsevier.com/content/article/PII:S0022039603001530?httpAccept=text/plain" elsevier_plain.txt
