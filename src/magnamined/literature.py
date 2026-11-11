"""
Europe PMC & PubMed REST API Client.
Prioritizes structured JATS XML retrieval and targeted gene snippet windowing.
"""

import re
import xml.etree.ElementTree as ET
from typing import Dict, Any, List, Optional
import requests


class EuropePMCClient:
    """Client for Europe PMC REST API with JATS XML extraction and snippet windowing."""

    BASE_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest"

    def __init__(self, session: Optional[requests.Session] = None):
        self.session = session or requests.Session()
        self.session.headers.update({
            "User-Agent": "Magna-Mined/0.1.0 (Research Literature Curation Engine)"
        })

    def search_literature(
        self,
        query: str,
        max_results: int = 25,
        open_access_only: bool = False
    ) -> List[Dict[str, Any]]:
        """
        Search Europe PMC for articles matching query.
        
        Args:
            query: Search query (e.g. gene identifiers or biological processes)
            max_results: Maximum records to return
            open_access_only: If True, restricts to freely accessible full text
        """
        full_query = query
        if open_access_only:
            full_query = f"({query}) AND (OPEN_ACCESS:y)"

        params = {
            "query": full_query,
            "format": "json",
            "pageSize": min(max_results, 100),
            "resultType": "core"
        }

        url = f"{self.BASE_URL}/search"
        try:
            response = self.session.get(url, params=params, timeout=15)
            response.raise_for_status()
            data = response.json()
            results = data.get("resultList", {}).get("result", [])
        except requests.RequestException as e:
            return []

        articles = []
        for item in results:
            articles.append({
                "pmid": item.get("pmid"),
                "pmcid": item.get("pmcid"),
                "doi": item.get("doi"),
                "title": item.get("title", "").strip().rstrip("."),
                "journal": item.get("journalTitle"),
                "pub_year": int(item.get("pubYear")) if str(item.get("pubYear", "")).isdigit() else None,
                "abstract": item.get("abstractText", ""),
                "is_open_access": item.get("isOpenAccess") == "Y",
                "in_epmc": item.get("inEPMC") == "Y",
                "pdf_url": self._pick_pdf_url(item)
            })
        return articles

    def stream_papers(
        self,
        query: str = '("Magnaporthe oryzae" OR "Pyricularia oryzae") AND (mutant OR knockout OR deletion OR CRISPR OR disruption OR overexpression OR RNAi OR complementation OR appressorium OR conidiation)',
        max_papers: int = 100,
        page_size: int = 50,
        open_access_only: bool = False,
        pub_year: Optional[int] = None,
        year_from: Optional[int] = None,
        year_to: Optional[int] = None,
        date_from: Optional[str] = None,
        date_to: Optional[str] = None,
    ):
        """
        Stream articles from Europe PMC using cursor-based pagination.
        Supports optional exact year (`pub_year`), inclusive year range (`year_from`..`year_to`),
        or ISO date window (`date_from`..`date_to` via `FIRST_PDATE`, for quarterly/6-month syncs).
        Pass `max_papers <= 0` to stream all matching results until exhausted.
        """
        full_query = query
        if date_from or date_to:
            start_d = date_from or "1900-01-01"
            end_d = date_to or "2100-12-31"
            full_query = f"({full_query}) AND (FIRST_PDATE:[{start_d} TO {end_d}])"
        elif pub_year is not None:
            full_query = f"({full_query}) AND (PUB_YEAR:{pub_year})"
        elif year_from is not None or year_to is not None:
            start_yr = year_from if year_from is not None else 1900
            end_yr = year_to if year_to is not None else 2100
            full_query = f"({full_query}) AND (PUB_YEAR:[{start_yr} TO {end_yr}])"
        if open_access_only:
            full_query = f"({full_query}) AND (OPEN_ACCESS:y)"

        cursor_mark = "*"
        yielded_count = 0
        unlimited = max_papers is None or max_papers <= 0

        while unlimited or yielded_count < max_papers:
            batch_limit = page_size if unlimited else min(page_size, max_papers - yielded_count)
            params = {
                "query": full_query,
                "format": "json",
                "pageSize": batch_limit,
                "cursorMark": cursor_mark,
                "resultType": "core"
            }
            url = f"{self.BASE_URL}/search"
            data = None
            for attempt in range(3):
                try:
                    response = self.session.get(url, params=params, timeout=45)
                    response.raise_for_status()
                    data = response.json()
                    break
                except requests.RequestException:
                    if attempt < 2:
                        import time
                        time.sleep(1.0 * (attempt + 1))
                    else:
                        data = None
            if data is None:
                break

            results = data.get("resultList", {}).get("result", [])
            if not results:
                break

            for item in results:
                yield {
                    "pmid": item.get("pmid"),
                    "pmcid": item.get("pmcid"),
                    "doi": item.get("doi"),
                    "title": item.get("title", "").strip().rstrip("."),
                    "journal": item.get("journalTitle"),
                    "pub_year": int(item.get("pubYear")) if str(item.get("pubYear", "")).isdigit() else None,
                    "abstract": item.get("abstractText", ""),
                    "is_open_access": item.get("isOpenAccess") == "Y",
                    "in_epmc": item.get("inEPMC") == "Y",
                    "pdf_url": self._pick_pdf_url(item)
                }
                yielded_count += 1
                if not unlimited and yielded_count >= max_papers:
                    break

            next_cursor = data.get("nextCursorMark")
            if not next_cursor or next_cursor == cursor_mark:
                break
            cursor_mark = next_cursor

    def search_gene(
        self,
        gene_id: str,
        gene_name: Optional[str] = None,
        aliases: Optional[List[str]] = None,
        max_results: int = 20
    ) -> List[Dict[str, Any]]:
        """
        Consolidated query for a specific Magnaporthe gene.
        Handles MGG systematic ID, common aliases, and organism prefixes.
        """
        terms = [f'"{gene_id}"']
        if gene_name:
            terms.append(f'"{gene_name}"')
            if not gene_name.upper().startswith("MO"):
                terms.append(f'"Mo{gene_name}"')
        if aliases:
            for a in aliases:
                quoted = f'"{a}"'
                if quoted not in terms:
                    terms.append(quoted)

        gene_clause = " OR ".join(terms)
        organism_clause = '("Magnaporthe oryzae" OR "Pyricularia oryzae" OR "rice blast")'
        query = f"({gene_clause}) AND {organism_clause}"

        return self.search_literature(query, max_results=max_results)

    @staticmethod
    def _pick_pdf_url(item: Dict[str, Any]) -> Optional[str]:
        """
        Pick an open-access publisher PDF link from a Europe PMC core record, if any.
        Europe PMC's own '?pdf=render' links are skipped: they reject scripted clients (HTTP 403).
        """
        urls = (item.get("fullTextUrlList") or {}).get("fullTextUrl", [])
        for u in urls:
            url = u.get("url") or ""
            if "europepmc.org" in url:
                continue
            if u.get("documentStyle") == "pdf" and u.get("availability") in ("Open access", "Free"):
                return url
        return None

    def find_unpaywall_pdf(self, doi: Optional[str], email: Optional[str] = None) -> Optional[str]:
        """
        Look up an open-access PDF location for a DOI via the Unpaywall API.
        Unpaywall requires a contact email; set MAGNAMINED_CONTACT_EMAIL or pass `email`.
        """
        import os
        email = email or os.getenv("MAGNAMINED_CONTACT_EMAIL")
        if not doi or not email:
            return None
        try:
            resp = self.session.get(
                f"https://api.unpaywall.org/v2/{doi}", params={"email": email}, timeout=15
            )
            if resp.status_code != 200:
                return None
            data = resp.json()
        except (requests.RequestException, ValueError):
            return None
        locations = [data.get("best_oa_location")] + (data.get("oa_locations") or [])
        for loc in locations:
            if loc and loc.get("url_for_pdf"):
                return loc["url_for_pdf"]
        return None

    def fetch_fulltext_pdf(
        self,
        pmcid: Optional[str] = None,
        pdf_url: Optional[str] = None,
        doi: Optional[str] = None,
        max_bytes: int = 30 * 1024 * 1024
    ) -> Optional[Dict[str, str]]:
        """
        Download an open-access PDF and convert it to layout-aware Markdown using
        pymupdf4llm (avoids two-column scrambling of naive PDF text extraction).

        Sources tried in order: explicit publisher pdf_url, then Unpaywall (by DOI).
        Returns the same section dict shape as fetch_fulltext_xml (results/methods
        are populated from headings where found; the full text is always in 'other'),
        or None if no PDF could be retrieved/converted.
        """
        urls: List[str] = []
        if pdf_url:
            urls.append(pdf_url)
        unpaywall_url = self.find_unpaywall_pdf(doi)
        if unpaywall_url and unpaywall_url not in urls:
            urls.append(unpaywall_url)

        pdf_bytes = None
        for url in urls:
            try:
                resp = self.session.get(url, timeout=20)
                if (
                    resp.status_code == 200
                    and resp.content[:5] == b"%PDF-"
                    and len(resp.content) <= max_bytes
                ):
                    pdf_bytes = resp.content
                    break
            except requests.RequestException:
                continue

        if not pdf_bytes:
            return None

        return self.pdf_to_sections(pdf_bytes)

    @classmethod
    def pdf_to_sections(cls, pdf_bytes: bytes) -> Optional[Dict[str, str]]:
        """Convert PDF bytes into the standard section dict (None on failure)."""
        markdown = cls.pdf_bytes_to_markdown(pdf_bytes)
        if not markdown:
            return None
        return cls.split_markdown_sections(markdown)

    @staticmethod
    def pdf_bytes_to_markdown(pdf_bytes: bytes) -> Optional[str]:
        """Convert raw PDF bytes to Markdown. Returns None if pymupdf4llm is unavailable/fails."""
        try:
            import pymupdf
            import pymupdf4llm
        except ImportError:
            return None
        try:
            doc = pymupdf.open(stream=pdf_bytes, filetype="pdf")
            try:
                try:
                    # Skip slow Tesseract OCR on figure images in digital PDFs
                    md = pymupdf4llm.to_markdown(doc, use_ocr=False, show_progress=False)
                except TypeError:
                    md = pymupdf4llm.to_markdown(doc)
            finally:
                doc.close()
            return md if md and md.strip() else None
        except Exception:
            return None

    @staticmethod
    def split_markdown_sections(markdown: str) -> Dict[str, str]:
        """Split converted PDF Markdown into the standard section dict using headings."""
        sections: Dict[str, str] = {
            "title": "", "abstract": "", "methods": "",
            "results": "", "discussion": "", "other": markdown
        }
        current = None
        buckets: Dict[str, List[str]] = {"abstract": [], "methods": [], "results": [], "discussion": []}
        for line in markdown.splitlines():
            stripped = line.strip()
            heading = re.match(r"^(?:#{1,6}\s*)?(?:\*\*\s*)?(?:\d+\.?\s*)?([A-Za-z &]+?)\**\s*$", stripped)
            if heading and len(stripped) < 60:
                name = heading.group(1).strip().lower()
                if name.startswith("abstract"):
                    current = "abstract"
                    continue
                if "method" in name or name.startswith("materials"):
                    current = "methods"
                    continue
                if name.startswith("result"):
                    current = "results"
                    continue
                if name.startswith("discussion") or name.startswith("conclusion"):
                    current = "discussion"
                    continue
                if name.startswith("reference") or name.startswith("acknowledg"):
                    current = None
                    continue
            if current:
                buckets[current].append(line)
        for key, lines in buckets.items():
            sections[key] = "\n".join(lines).strip()
        return sections

    def fetch_fulltext_xml(self, pmcid: str) -> Optional[Dict[str, str]]:
        """
        Fetch and parse clean JATS XML for an open access article.
        Returns dictionary of section_name -> text.
        """
        clean_pmcid = pmcid if pmcid.startswith("PMC") else f"PMC{pmcid}"
        url = f"{self.BASE_URL}/{clean_pmcid}/fullTextXML"

        try:
            response = self.session.get(url, timeout=20)
            if response.status_code != 200 or not response.content:
                return None

            root = ET.fromstring(response.content)
            sections: Dict[str, str] = {
                "title": "",
                "abstract": "",
                "methods": "",
                "results": "",
                "discussion": "",
                "other": ""
            }

            # Extract title
            title_node = root.find(".//article-title")
            if title_node is not None:
                sections["title"] = "".join(title_node.itertext()).strip()

            # Extract abstract
            for abs_node in root.findall(".//abstract"):
                sections["abstract"] += " ".join("".join(abs_node.itertext()).split()) + "\n"

            # Extract body sections
            for sec in root.findall(".//body//sec"):
                sec_type = sec.attrib.get("sec-type", "").lower()
                title_elem = sec.find("title")
                title_text = "".join(title_elem.itertext()).lower() if title_elem is not None else ""
                body_text = " ".join("".join(sec.itertext()).split())

                if "method" in sec_type or "method" in title_text or "material" in title_text:
                    sections["methods"] += body_text + "\n"
                elif "result" in sec_type or "result" in title_text:
                    sections["results"] += body_text + "\n"
                elif "discuss" in sec_type or "discuss" in title_text:
                    sections["discussion"] += body_text + "\n"
                else:
                    sections["other"] += body_text + "\n"

            # Also capture floating tables and back-matter notes/appendices (accession tables & symbols)
            for extra_tag in (".//floats-group//table-wrap", ".//back//notes", ".//back//app-group", ".//back//ack"):
                for node in root.findall(extra_tag):
                    extra_text = " ".join("".join(node.itertext()).split())
                    if extra_text:
                        sections["other"] += extra_text + "\n"

            return sections

        except Exception:
            return None

    @staticmethod
    def extract_targeted_snippets(
        text_content: str,
        gene_id: str,
        aliases: Optional[List[str]] = None,
        window_chars: int = 1500
    ) -> List[str]:
        """
        Extract focused text windows around mentions of target gene or mutant.
        Prevents context-window bloat by discarding irrelevant sections.
        """
        search_terms = [re.escape(gene_id)]
        if aliases:
            for alias in aliases:
                if alias:
                    search_terms.append(re.escape(alias))
                    if not alias.upper().startswith("MO"):
                        search_terms.append(re.escape(f"Mo{alias}"))
                    search_terms.append(re.escape(f"Δ{alias}"))
                    search_terms.append(re.escape(f"delta-{alias}"))

        pattern = re.compile(r"|".join(search_terms), re.IGNORECASE)
        snippets = []

        for match in pattern.finditer(text_content):
            start = max(0, match.start() - window_chars // 2)
            end = min(len(text_content), match.end() + window_chars // 2)
            snippet = text_content[start:end].strip()
            snippets.append(snippet)

        return snippets
