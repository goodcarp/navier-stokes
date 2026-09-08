"""Download only source text needed for a static theorem-scope inspection.

No downloaded source is imported or executed. URLs are pinned to one commit.
"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen
import hashlib
import json

COMMIT = "d0124689230b58b4f86e7b90ac59de06404b3b6b"
FILES = [
    "Challenge.lean", "Solution.lean", "README.md", "formalization.yaml",
    "comparator.json", "scripts/PrintAxioms.lean", "EulerBlowup/Main.lean",
    "EulerBlowup/Theorem01Prime.lean", "EulerBlowup/Theorem01.lean",
]
OUT = Path("work/sources/proof-scope")


def fetch(rel):
    url = f"https://raw.githubusercontent.com/tristanbuckmaster/fluid_lean/{COMMIT}/euler-blowup/{rel}"
    with urlopen(url, timeout=30) as response:
        data = response.read()
    dest = OUT / rel
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(data)
    return {"path": rel, "url": url, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


if __name__ == "__main__":
    with ThreadPoolExecutor(max_workers=4) as pool:
        records = list(pool.map(fetch, FILES))
    (OUT / "manifest.json").write_text(json.dumps({"commit": COMMIT, "files": records}, indent=2) + "\n")
    print(f"Downloaded {len(records)} source files, {sum(r['bytes'] for r in records)} bytes. Nothing executed.")
