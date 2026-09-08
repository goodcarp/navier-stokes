#!/bin/bash
# Downloads the open-access citing works of BFG (OpenAlex cites:W2608160632)
cd "$(dirname "$0")/.."
UA="Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/126 Safari/537.36"
dl(){ curl -sL -A "$UA" -o "pdf/$1" "$2"; f=$(file -b "pdf/$1"); echo "$1 :: $f"; case "$f" in PDF*) pdftotext -layout "pdf/$1" "txt/${1%.pdf}.txt";; esac; }
dl cite-grujic-xu-asymptotic-criticality-1911.00974.pdf https://arxiv.org/pdf/1911.00974
dl cite-grujic-xu-hyperdissipation-2012.05692.pdf https://arxiv.org/pdf/2012.05692
dl cite-bradshaw-grujic-oscillations-1801.09040.pdf https://arxiv.org/pdf/1801.09040
dl cite-local-pressure-expansion-2001.11526.pdf https://arxiv.org/pdf/2001.11526
dl cite-axisym-review-2101.04905.pdf https://arxiv.org/pdf/2101.04905
# (blocked, saved as html/) dl cite-sparseness-regularity-iop.pdf https://iopscience.iop.org/article/10.1088/1361-6544/ac62de/pdf
dl cite-separation-iop.pdf https://iopscience.iop.org/article/10.1088/1361-6544/ad68b9/pdf
dl cite-sqg-2025-springer.pdf https://link.springer.com/content/pdf/10.1007/s00021-025-00947-x.pdf
dl cite-asymptotic-criticality-2024-springer.pdf https://link.springer.com/content/pdf/10.1007/s00021-024-00888-x.pdf
