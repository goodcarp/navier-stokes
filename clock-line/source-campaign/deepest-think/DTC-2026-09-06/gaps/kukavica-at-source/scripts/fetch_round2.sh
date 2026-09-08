#!/bin/bash
# Round 2: the remaining routes named in the brief, logged for completeness after the
# author-homepage (Wayback) route already succeeded.  Sci-Hub is NOT used.
set -u
UA='Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/125.0 Safari/537.36'
LOG=log/round2.tsv
echo -e "route\turl\thttp\tbytes\tcontent_type\tsaved_as" > "$LOG"
try () { n="$1"; u="$2"; o="$3"; shift 3
  r=$(curl -sSL -m 30 -A "$UA" "$@" -o "raw/$o" -w "%{http_code}\t%{size_download}\t%{content_type}" "$u" 2>>log/round2.err || echo -e "ERR\t0\t-")
  echo -e "$n\t$u\t$r\traw/$o" >> "$LOG"; }
try arxiv_author_search "http://export.arxiv.org/api/query?search_query=au:%22Kukavica%22+AND+abs:%22bounded+initial+data%22&max_results=30" arxiv_kukavica.xml
try arxiv_listing_2003 "http://export.arxiv.org/api/query?search_query=au:%22Kukavica_I%22&start=0&max_results=100&sortBy=submittedDate&sortOrder=ascending" arxiv_kukavica_all.xml
try researchgate "https://www.researchgate.net/search/publication?q=On%20local%20uniqueness%20of%20weak%20solutions%20of%20the%20Navier-Stokes%20system%20with%20bounded%20initial%20data" rg.html
try ia_scholar "https://scholar.archive.org/search?q=%22On+local+uniqueness+of+weak+solutions+of+the+Navier-Stokes+system+with+bounded+initial+data%22" ia_scholar.html
try openalex_oa "https://api.openalex.org/works/doi:10.1016/S0022-0396(03)00153-0" openalex.json
try usc_repository "https://digitallibrary.usc.edu/search?q=Kukavica%20Navier-Stokes" usc_repo.html
try wayback_index_live "https://web.archive.org/web/20170627135430/http://www-bcf.usc.edu/~kukavica/pdf/index.html" wb_index_rendered.html
