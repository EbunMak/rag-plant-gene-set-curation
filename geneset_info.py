#!/usr/bin/env python3
"""
Fetch gene set info from PlantGSEA for Triticum aestivum.

Usage:
    python geneset_info.py <gene_sets.gmt>
"""

import requests
import json
import os
import time
from typing import List, Dict, Any
from bs4 import BeautifulSoup
import sys
import urllib3

# Suppress SSL warnings since we're disabling verification for this server
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# Mimic a browser to avoid 403 rejections
HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
        "(KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
    )
}

def read_gmt(gmt_file: str) -> Dict[str, List[str]]:
    gene_sets = {}
    with open(gmt_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            gene_set_name = parts[0]
            genes = parts[2:]
            gene_sets[gene_set_name] = genes
    return gene_sets


def fetch_geneset_info(gene_set_name: str, session: requests.Session) -> Dict[str, Any]:
    url = (
        f"https://systemsbiology.cau.edu.cn/PlantGSEAv2/gene_set_detail.php"
        f"?species=Triticum%20aestivum&geneset={gene_set_name}"
    )
    try:
        response = session.get(url, headers=HEADERS, verify=False, timeout=30)
        if response.status_code == 200:
            soup = BeautifulSoup(response.text, "html.parser")
            table = soup.find("table", {"class": "result"})
            if table:
                rows = table.find_all("tr")
                info = {}
                # extract only the rows "Standard Gene Set Name", "Gene Set type"  and "Full Description/Abstract"
                target_rows = ["Standard Gene Set Name", "Gene Set type", "Full Description/Abstract"]
                for row in rows:
                    cols = row.find_all("td")
                    if len(cols) == 2:
                        key = cols[0].text.strip()
                        value = cols[1].text.strip()
                        if key in target_rows:
                            info[key] = value
                return info
            else:
                print(f"  No result table found for: {gene_set_name}")
                return {}
        else:
            print(f"  HTTP {response.status_code} for: {gene_set_name}")
            return {}
    except requests.exceptions.SSLError as e:
        print(f"  SSL error for {gene_set_name}: {e}")
        return {}
    except requests.exceptions.ConnectionError as e:
        print(f"  Connection error for {gene_set_name}: {e}")
        return {}
    except requests.exceptions.Timeout:
        print(f"  Timeout for: {gene_set_name}")
        return {}


def main():
    if len(sys.argv) < 2:
        print("Usage: python geneset_info.py <gene_sets.gmt>")
        sys.exit(1)

    gmt_file = sys.argv[1]
    if not os.path.isfile(gmt_file):
        print(f"Error: GMT file not found: {gmt_file}")
        sys.exit(1)

    print(f"Reading gene sets from: {gmt_file}")
    gene_sets = read_gmt(gmt_file)
    print(f"Found {len(gene_sets)} gene sets")

    all_info = {}

    # Use a session for connection reuse
    with requests.Session() as session:
        for i, gene_set_name in enumerate(list(gene_sets.keys()), start=1):  # Limit to first 10 gene sets for testing
            print(f"[{i}/{len(gene_sets)}] Fetching: {gene_set_name}")
            info = fetch_geneset_info(gene_set_name, session)
            all_info[gene_set_name] = info
            time.sleep(0.5)  # Be polite to the server

    # Save to JSON
    output_path = "geneset_info.json"
    with open(output_path, "w") as f:
        json.dump(all_info, f, indent=4)

    print(f"\nDone. Info for {len(all_info)} gene sets saved to: {output_path}")


if __name__ == "__main__":
    main()