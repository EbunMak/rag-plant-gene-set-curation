import requests
import time
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

LOG_FILE = "abstract_data.txt"

HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
        "(KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
    )
}

def make_session() -> requests.Session:
    """Create a session with automatic retry on connection errors."""
    session = requests.Session()
    session.headers.update(HEADERS)
    retry = Retry(
        total=5,                  # retry up to 5 times
        backoff_factor=2,         # wait 2, 4, 8, 16, 32 seconds between retries
        status_forcelist=[429, 500, 502, 503, 504],  # retry on these HTTP codes
        allowed_methods=["GET"]
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def log_abstract_data(query, total_pages, page_limit, total_abstracts):
    with open(LOG_FILE, "a") as f:
        f.write(f"Query: {query}\n")
        f.write(f"Total available pages: {total_pages}\n")
        f.write(f"Page limit used: {page_limit}\n")
        f.write(f"Total abstracts retrieved: {total_abstracts}\n")
        f.write("-" * 40 + "\n")


def safe_get(session: requests.Session, url: str, params: dict, retries: int = 5) -> requests.Response:
    """
    GET with manual retry on ConnectionResetError (104),
    which isn't caught by urllib3's Retry by default.
    """
    delay = 2
    for attempt in range(1, retries + 1):
        try:
            r = session.get(url, params=params, timeout=30)
            r.raise_for_status()
            return r
        except (requests.exceptions.ConnectionError, requests.exceptions.ChunkedEncodingError) as e:
            if attempt == retries:
                raise
            print(f"  Connection error (attempt {attempt}/{retries}): {e}. Retrying in {delay}s...")
            time.sleep(delay)
            delay *= 2  # exponential backoff
        except requests.exceptions.HTTPError as e:
            if r.status_code == 429:
                print(f"  Rate limited (429). Waiting {delay}s before retry...")
                time.sleep(delay)
                delay *= 2
            else:
                raise


class Pubtator:
    BASE_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api"

    @staticmethod
    def find_entity_ID(entity_details: str, bioconcept: str = None, limit: int = 100):
        session = make_session()
        url = f"{Pubtator.BASE_URL}/entity/autocomplete/"
        params = {
            "query": entity_details,
            "concept": bioconcept,
            "limit": limit
        }
        r = safe_get(session, url, {k: v for k, v in params.items() if v is not None})
        return r.json()

    @staticmethod
    def find_related_entity(entity_id: str, relation_type: str = None, entity_type: str = None):
        session = make_session()
        url = f"{Pubtator.BASE_URL}/relations"
        params = {
            "e1": entity_id,
            "type": relation_type,
            "e2": entity_type
        }
        r = safe_get(session, url, {k: v for k, v in params.items() if v is not None})
        return r.json()

    @staticmethod
    def search_pubtator_ID(query: str = "", relation: str = None, limit: int = 25):
        session = make_session()
        results = []
        page = 1
        num_of_pages = 1
        page_limit = limit
        query = query + " in T aestivum OR Triticum aestivum OR wheat AND genes"
        total_available_pages = 1

        while True:
            url = f"{Pubtator.BASE_URL}/search/"
            params = {
                "text": relation if relation else query,
                "page": page
            }

            r = safe_get(session, url, params)
            data = r.json()

            if page == 1:
                total_available_pages = data.get("total_pages", 1)
                num_of_pages = min(total_available_pages, page_limit)
                print(f"Parsing through {num_of_pages} pages")

            curr_result = data.get("results", [])
            pmids = [item["pmid"] for item in curr_result if "pmid" in item]
            results.extend(pmids)

            page += 1
            if page > num_of_pages:
                break

            time.sleep(1.0)  # increased from 0.3 to 1.0 to avoid connection resets

        log_abstract_data(
            query=query,
            total_pages=total_available_pages,
            page_limit=page_limit,
            total_abstracts=len(results)
        )

        return results

    @staticmethod
    def export_abstract(pmid: str, check_for_genes=True):
        session = make_session()
        url = f"{Pubtator.BASE_URL}/publications/export/biocjson?pmids={pmid}"
        r = safe_get(session, url, params={})
        data = r.json()

        result = {
            "pmid": pmid,
            "title": None,
            "journal": None,
            "abstract": None,
            "genes": []
        }

        pub = data["PubTator3"][0]
        passages = pub.get("passages", [])
        result["journal"] = pub.get("journal", None)

        for p in passages:
            p_type = p.get("infons", {}).get("type")
            if p_type == "title":
                result["title"] = p.get("text")
            elif p_type == "abstract":
                result["abstract"] = p.get("text")

            for ann in p.get("annotations", []):
                infons = ann.get("infons", {})
                if infons.get("type", "").lower() == "gene":
                    gene_entry = {
                        "name": infons.get("name"),
                        "identifier": infons.get("identifier"),
                        "accession": infons.get("accession"),
                        "text": ann.get("text"),
                        "location": ann.get("locations", [{}])[0].get("offset", None)
                    }
                    result["genes"].append(gene_entry)

        if check_for_genes:
            if not result["genes"]:
                return None
        return result
    
    @staticmethod
    def export_full_text(pmid: str, check_for_genes=True):
        """
        Retrieve metadata + full text annotations for a given PMID using PubTator3 full text export.
        
        :param pmid: PubMed ID
        :return: dict with title, journal, sections (full text organized by section), and annotations
        """
        url = f"{Pubtator.BASE_URL}/publications/export/biocjson?pmids={pmid}&full=true"
        r = requests.get(url)
        r.raise_for_status()
        data = r.json()

        # Base result container
        result = {
            "pmid": pmid,
            "title": None,
            "journal": None,
            "pmcid": None,
            "sections": {}
        }

        # Extract core metadata
        pub = data["PubTator3"][0]
        pub_infons = pub.get("infons", {})
        result["journal"] = pub_infons.get("journal")
        result["pmcid"] = pub_infons.get("article-id_pmc")
        
        passages = pub.get("passages", [])
        
        # Extract title
        for p in passages:
            p_infons = p.get("infons", {})
            section_type = p_infons.get("section_type")
            p_type = p_infons.get("type")
            
            if section_type == "TITLE" or p_type == "title":
                result["title"] = p.get("text")
                break

        # Extract and organize only abstract, results and discussion sections (if available)
        full_text_sections = {
            "ABSTRACT": "Abstract",
            "RESULTS": "Results",
            "DISCUSS": "Discussion"
        }
        
        full_text_parts = []
        
        for p in passages:
            p_infons = p.get("infons", {})
            section_type = p_infons.get("section_type")

            if section_type in ["TITLE", "front"]:
                continue
                
            #  main full text sections
            if section_type in full_text_sections:
                section_name = full_text_sections[section_type]
                result["sections"][section_type] = {
                    "title": section_name,
                    "text": p.get("text", "")
                }
                full_text_parts.append(f"\n\n## {section_name}\n{p.get('text', '')}")

            elif section_type and len(p.get("text", "")) > 50:
                # print(f"Found additional section type: {section_type} with text length {len(p.get('text', ''))}. Adding to sections.")
                result["sections"][section_type] = {
                    "title": section_type,
                    "text": p.get("text", "")
                }
                full_text_parts.append(f"\n\n## {section_type}\n{p.get('text', '')}")
        
        # Concatenate full text
        # result["full_text"] = "".join(full_text_parts).strip()
        
        return result