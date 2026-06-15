"""
scrape_funplantgenes.py
-----------------------
Scrapes funplantgenes.henau.edu.cn for all Triticum aestivum gene entries
and builds a JSON mapping of functional name -> TraesCS ID.

Output:
    out/funplantgenes_mapping.json

Usage:
    python scrape_funplantgenes.py
"""

import json
import re
import time
from pathlib import Path

import requests
from bs4 import BeautifulSoup

BASE_URL    = "https://funplantgenes.henau.edu.cn/categories/triticum-aestivum/page/{}/"
FIRST_PAGE  = "https://funplantgenes.henau.edu.cn/categories/triticum-aestivum/"
TOTAL_PAGES = 164
OUTPUT      = Path("out/funplantgenes_mapping.json")
DELAY       = 0.5   # seconds between requests to be polite


def parse_page(html: str) -> list[tuple[str, str]]:
    """
    Extract (gene_name, traescs_id) pairs from a listing page.
    Each entry looks like: 'ALI-1 Awn Length Inhibitor 1 ; TraesCS5A02G542800 ; Triticum aestivum'
    """
    soup    = BeautifulSoup(html, "html.parser")
    results = []

    for link in soup.find_all("a", href=re.compile(r"/genes/triticum_aestivum/")):
        text = link.get_text(" ", strip=True)
        # Extract TraesCS ID — matches TraesCS...02G... pattern
        traescs_match = re.search(r"(TraesCS\w+)", text)
        if not traescs_match:
            continue
        traescs_id = traescs_match.group(1)

        # Gene name is the first token before any whitespace/description
        # e.g. "ALI-1Awn Length..." or "FT-A1FLOWERING LOCUS T"
        # The href slug is the cleanest source for the name
        slug      = link["href"].rstrip("/").split("/")[-1]
        gene_name = slug.upper().replace("-", "").replace("_", "")

        # Also store the raw text name (first word-like token)
        raw_name_match = re.match(r"^([A-Za-z0-9\-\.]+)", text)
        raw_name = raw_name_match.group(1) if raw_name_match else slug

        results.append((raw_name, traescs_id))

    return results


def scrape():
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    mapping: dict[str, str] = {}
    skipped = 0

    for page_num in range(1, TOTAL_PAGES + 1):
        url = FIRST_PAGE if page_num == 1 else BASE_URL.format(page_num)
        try:
            r = requests.get(url, timeout=15)
            r.raise_for_status()
        except requests.RequestException as e:
            print(f"  [error] page {page_num}: {e}")
            continue

        pairs = parse_page(r.text)
        for name, traescs_id in pairs:
            if name not in mapping:
                mapping[name] = traescs_id
            # Store all aliases — if name already exists keep first hit
        
        found = len(pairs)
        skipped += (0 if found else 1)
        print(f"  page {page_num:>3}/{TOTAL_PAGES}  {found:>3} entries  (total so far: {len(mapping)})")

        time.sleep(DELAY)

    with open(OUTPUT, "w") as f:
        json.dump(mapping, f, indent=2, sort_keys=True)

    print(f"\n[done] {len(mapping)} gene mappings → {OUTPUT}")
    print(f"       {skipped} pages had no parseable entries")


if __name__ == "__main__":
    scrape()
