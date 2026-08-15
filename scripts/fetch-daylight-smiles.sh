#!/bin/sh
# Fetch the Daylight SMILES documentation pages and extract candidate SMILES
# strings for the local dialect instrument.
#
# Provenance: Daylight Chemical Information Systems, SMILES theory manual and
# tutorial examples (daylight.com). The pages are copyrighted documentation
# staged for local validation only and must never be committed or
# redistributed (discussion/194 S6b2). The staging directory lies under
# materials/, which is covered by the repository .gitignore.
set -eu

staging="$(dirname "$0")/../materials/formats/daylight"
base="https://www.daylight.com"

mkdir -p "$staging"
curl -fsS "$base/dayhtml/doc/theory/theory.smiles.html" -o "$staging/theory.smiles.html"
curl -fsS "$base/dayhtml_tutorials/languages/smiles/smiles_examples.html" \
    -o "$staging/smiles_examples.html"

python3 - "$staging" <<'PY'
import html
import re
import sys

staging = sys.argv[1]
candidates = []
token = re.compile(r'^[A-Za-z0-9@+\-\[\]()=#$:/\\%.*]{2,}$')
for page in ("theory.smiles.html", "smiles_examples.html"):
    text = open(f"{staging}/{page}", encoding="utf-8", errors="replace").read()
    # Markup-delimited candidates: tt/b/pre/code blocks and table cells.
    blocks = re.findall(
        r"<(?:tt|b|code|pre|td)[^>]*>(.*?)</(?:tt|b|code|pre|td)>",
        text,
        re.IGNORECASE | re.DOTALL,
    )
    for block in blocks:
        block = re.sub(r"<[^>]+>", " ", block)
        for word in html.unescape(block).split():
            if token.match(word):
                candidates.append(f"{page}\t{word}")
seen = set()
with open(f"{staging}/candidates.tsv", "w", encoding="utf-8") as out:
    for line in candidates:
        if line not in seen:
            seen.add(line)
            out.write(line + "\n")
print(f"{len(seen)} candidates")
PY
