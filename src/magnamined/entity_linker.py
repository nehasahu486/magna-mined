"""
Entity Disambiguation & Alias Resolution Engine for Magnaporthe oryzae.
Maps gene symbols, organism-prefixed names (MoGene), mutant designations (Δgene, delta-gene),
and legacy annotation IDs to canonical MGG_##### identifiers.
"""

import re
from typing import Dict, List, Optional, Set, Tuple
from dataclasses import dataclass, field

from .database import DatabaseManager


# Common biological and experimental stopwords to reject from entity recognition
BIOLOGICAL_STOPWORDS = {
    "DNA", "RNA", "MRNA", "TRNA", "RRNA", "CDNA", "DSRNA", "SIRNA",
    "PCR", "QPCR", "RT-PCR", "GFP", "RFP", "YFP", "CFP", "FLAG", "HA", "MYC", "GST",
    "WT", "WILDTYPE", "WILD-TYPE", "VECTOR", "CONTROL", "MOCK",
    "YPD", "PDA", "CM", "MM", "LB", "OD", "OD600", "RPM", "MIC",
    "H2O", "H2O2", "SDS", "PBS", "EDTA", "PEG", "DAPI", "DTT", "DMSO",
    "BLAST", "NCBI", "PDB", "GO", "KEGG", "PFAM", "COG",
    "ACTIN", "TUBULIN", "GAPDH", "UBIQUITIN", "HIS", "TAG",
    "TOTAL", "HIGH", "LOW", "FOLD", "LEVEL", "RATE", "TIME", "DAYS", "HOURS"
}

# Generic prose words that can appear next to (MGG_#####) without being gene symbols
INLINE_NON_GENE_WORDS = BIOLOGICAL_STOPWORDS | {
    "GENE", "GENES", "PROTEIN", "PROTEINS", "HOMOLOG", "HOMOLOGS", "HOMOLOGUE",
    "ORTHOLOG", "ORTHOLOGS", "PARALOG", "FACTOR", "FACTORS", "KINASE", "ENZYME",
    "DOMAIN", "DOMAINS", "MUTANT", "MUTANTS", "STRAIN", "STRAINS", "ACCESSION",
    "GENBANK", "LOCUS", "NUMBER", "ID", "FIGURE", "FIG", "TABLE", "LIKE",
    "PUTATIVE", "CANDIDATE", "HYPOTHETICAL", "PREDICTED", "UNKNOWN", "ENCODING",
    "CLUSTER", "REGION", "SEQUENCE", "CLONE", "CONTIG", "CHROMOSOME", "ORYZAE"
}

# Known legacy Broad Institute annotation IDs mapped to MGG 8.0 canonical IDs
LEGACY_MGG_MAP: Dict[str, str] = {
    "MGG_00692": "MGG_09565",  # PMK1 (Broad v6 -> MGG 8.0)
    "MGG_01130": "MGG_09898",  # MAC1 (Broad v6 -> MGG 8.0)
    "MGG_02846": "MGG_12958",  # MST12 (Broad v6 -> MGG 8.0)
    "MGG_00659": "MGG_06368",  # CPKA (Broad v6 -> MGG 8.0)
}


@dataclass
class GeneMention:
    """A detected mention of a gene or mutant in text."""
    raw_text: str
    start_char: int
    end_char: int
    canonical_id: str
    canonical_name: Optional[str] = None
    confidence: float = 1.0
    is_mutant: bool = False


class EntityLinker:
    """
    High-throughput gene entity disambiguation and alias resolution engine.
    Pre-indexes all 12,500+ M. oryzae genes, standard names, and morphological variations.
    """

    def __init__(self, db: Optional[DatabaseManager] = None):
        self.db = db or DatabaseManager()
        self.exact_mgg: Dict[str, Dict[str, str]] = {}  # "MGG_01001" -> {"gene": ..., "gene_name": ...}
        self.symbol_to_mggs: Dict[str, List[str]] = {}  # "ENG1" -> ["MGG_01001"]
        self.mgg_to_symbols: Dict[str, List[str]] = {}  # "MGG_01001" -> ["ENG1", "MoENG1", ...]
        self.compiled_regex: Optional[re.Pattern] = None

        self._build_index()

    def _build_index(self) -> None:
        """Build in-memory lookup index from database."""
        with self.db.get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT gene, gene_name FROM genes")
            rows = cursor.fetchall()

        all_symbols: Set[str] = set()

        for row in rows:
            gene = row["gene"].strip().upper()
            gene_name = (row["gene_name"] or "").strip()

            self.exact_mgg[gene] = {
                "gene": row["gene"],
                "gene_name": row["gene_name"]
            }

            aliases_for_this_gene = []

            if gene_name:
                # Handle slash, semicolon, or comma separated aliases (e.g. ZIP8/CREB4/ZIP6)
                raw_tokens = re.split(r"[/,;]", gene_name)
                for t in raw_tokens:
                    cleaned = t.strip()
                    if (
                        cleaned
                        and len(cleaned) >= 2
                        and not cleaned.isdigit()
                        and cleaned.upper() not in BIOLOGICAL_STOPWORDS
                    ):
                        upper_symbol = cleaned.upper()
                        all_symbols.add(cleaned)
                        aliases_for_this_gene.append(cleaned)

                        if upper_symbol not in self.symbol_to_mggs:
                            self.symbol_to_mggs[upper_symbol] = []
                        if gene not in self.symbol_to_mggs[upper_symbol]:
                            self.symbol_to_mggs[upper_symbol].append(gene)

            self.mgg_to_symbols[gene] = aliases_for_this_gene

        # Also load persisted synonyms from `gene_aliases` (e.g. UniProt-synced aliases)
        for gene, db_aliases in self.db.get_all_gene_aliases().items():
            if gene not in self.exact_mgg:
                continue
            bucket = self.mgg_to_symbols.setdefault(gene, [])
            for al in db_aliases:
                cleaned = al.strip()
                if (
                    cleaned
                    and len(cleaned) >= 2
                    and not cleaned.isdigit()
                    and cleaned.upper() not in BIOLOGICAL_STOPWORDS
                ):
                    upper_symbol = cleaned.upper()
                    all_symbols.add(cleaned)
                    if cleaned not in bucket:
                        bucket.append(cleaned)
                    self.symbol_to_mggs.setdefault(upper_symbol, [])
                    if gene not in self.symbol_to_mggs[upper_symbol]:
                        self.symbol_to_mggs[upper_symbol].append(gene)
                    # Also index stripped symbol (without Mo/Po prefix) so normalize_token resolves ΔMoXxx
                    stripped_sym = re.sub(r"^(?:Mo|Po)[-_]?", "", cleaned, flags=re.IGNORECASE)
                    if stripped_sym and len(stripped_sym) >= 3 and stripped_sym.upper() not in BIOLOGICAL_STOPWORDS:
                        all_symbols.add(stripped_sym)
                        su = stripped_sym.upper()
                        self.symbol_to_mggs.setdefault(su, [])
                        if gene not in self.symbol_to_mggs[su]:
                            self.symbol_to_mggs[su].append(gene)
            if not self.exact_mgg[gene].get("gene_name") and bucket:
                self.exact_mgg[gene]["gene_name"] = bucket[0]

        # Also register known legacy mappings
        for old_id, new_id in LEGACY_MGG_MAP.items():
            if new_id in self.exact_mgg:
                target_name = self.exact_mgg[new_id].get("gene_name")
                self.exact_mgg[old_id] = {
                    "gene": new_id,
                    "gene_name": target_name
                }

        # Build combined regex for text scanning
        sorted_symbols = sorted(all_symbols, key=len, reverse=True)
        escaped_symbols = "|".join(re.escape(s) for s in sorted_symbols) if sorted_symbols else ""

        # Matches:
        # 1. MGG IDs: MGG_#####, MGG#####, MGG-##### (with optional T0 suffix or delta prefix)
        # 2. Gene symbols with optional organism prefix (Mo, Po) and mutant prefix/suffix (Δ, delta-, -null)
        mgg_pattern = r"(?:[Δ∆]|delta[-_]?)?MGG[_-]?\d{5}(?:[T_\-]\d+)?"
        if escaped_symbols:
            symbol_pattern = (
                r"(?:[Δ∆]|delta[-_]?)?(?:Mo|Po|MO|PO)?(?:[-_])?"
                r"(?:" + escaped_symbols + r")"
                r"(?:[Δ∆]|[-_]null|[-_]mutant|[-_]KO)?"
            )
            full_pattern = rf"\b(?:{mgg_pattern})\b|(?<!\w)(?:{symbol_pattern})(?!\w)"
        else:
            full_pattern = rf"\b(?:{mgg_pattern})\b"

        self.compiled_regex = re.compile(full_pattern, re.IGNORECASE)

    def normalize_token(self, token: str) -> str:
        """Strip organism prefixes, mutant indicators, and normalize MGG formats."""
        s = token.strip()
        # Remove delta prefix (handles both Greek Δ U+0394 and math ∆ U+2206)
        s = re.sub(r"^(?:[Δ∆]|delta[-_]?)", "", s, flags=re.IGNORECASE)
        # Remove mutant suffix
        s = re.sub(r"(?:[Δ∆]|[-_]null|[-_]mutant|[-_]ko)$", "", s, flags=re.IGNORECASE)

        # Check if it's an MGG ID: MGG_#####, MGG-#####, MGG#####, optionally with transcript suffix
        mgg_match = re.match(r"^MGG[_-]?(\d{5})(?:[T_\-\.]\d+)?$", s, re.IGNORECASE)
        if mgg_match:
            return f"MGG_{mgg_match.group(1)}"

        # Remove Mo / Po prefix (case-insensitive, optional hyphen)
        s = re.sub(r"^(?:Mo|Po|MO|PO)[-_]?", "", s, flags=re.IGNORECASE)
        return s.strip()

    def resolve(self, mention: str, context: Optional[str] = None) -> Optional[str]:
        """
        Resolve a raw mention (e.g. 'Δmopmk1', 'MoENG1', 'MGG_01001T0', 'MGG-01001')
        to canonical MGG_##### ID.
        If mention is ambiguous (e.g. ZIP9), uses context string to disambiguate.
        """
        if not mention:
            return None

        # Check legacy map first
        mention_upper = mention.strip().upper()
        if mention_upper in LEGACY_MGG_MAP:
            return LEGACY_MGG_MAP[mention_upper]

        norm = self.normalize_token(mention).upper()

        # Check if it's already an MGG ID
        if norm in self.exact_mgg:
            return self.exact_mgg[norm]["gene"]

        # Check legacy map on normalized token
        if norm in LEGACY_MGG_MAP:
            return LEGACY_MGG_MAP[norm]

        # Check symbol map
        if norm in self.symbol_to_mggs:
            candidates = self.symbol_to_mggs[norm]
            if len(candidates) == 1:
                return candidates[0]

            # Ambiguity: try to disambiguate using context
            if context:
                context_upper = context.upper()
                for c in candidates:
                    if c in context_upper:
                        return c

            # Default to first candidate if disambiguation inconclusive
            return candidates[0]

        return None

    @staticmethod
    def _is_plausible_gene_symbol(sym: str) -> bool:
        """Check whether a candidate inline token looks like a fungal gene symbol (e.g. MoAcb1, MoPct1, Clp1)."""
        if not sym or len(sym) < 2 or len(sym) > 12 or sym.isdigit():
            return False
        if sym.upper() in INLINE_NON_GENE_WORDS:
            return False
        if re.match(r"^MGG[_-]?\d{5}", sym, re.IGNORECASE):
            return False
        # Must look like a symbol: Mo/Po/Mr prefix, >=2 uppercase letters, or ends with digit(s)
        has_org_prefix = bool(re.match(r"^(?:Mo|Po|Mr)[A-Z]", sym))
        upper_count = sum(1 for c in sym if c.isupper())
        ends_with_digit = bool(re.search(r"\d+$", sym))
        return has_org_prefix or upper_count >= 2 or ends_with_digit

    def extract_inline_aliases(self, text: str) -> Dict[str, List[str]]:
        """
        Discover paper-local 'GeneName (MGG_#####)' or 'MGG_##### (GeneName)' bindings
        for genes that are unnamed in the baseline catalog (e.g. 'MoAcb1 (MGG_06177)').
        Returns {canonical_mgg_id: [alias, stripped_alias, ...]}.
        """
        if not text:
            return {}

        patterns = [
            # Symbol (MGG_12345) or Symbol [MGG_12345] or Symbol (gene ID: MGG_12345)
            re.compile(
                r"\b([A-Za-z][A-Za-z0-9-]{1,11})\s*[\(\[]\s*(?:(?:gene|locus|id|accession)\s*[:#]?\s*)?(MGG[_-]?\d{5})(?:[T_\-]\d+)?\s*[\)\]]",
                re.IGNORECASE,
            ),
            # MGG_12345 (Symbol) or MGG_12345 [Symbol]
            re.compile(
                r"\b(MGG[_-]?\d{5})(?:[T_\-]\d+)?\s*[\(\[]\s*([A-Za-z][A-Za-z0-9-]{1,11})\s*[\)\]]",
                re.IGNORECASE,
            ),
            # Symbol/MGG_12345 or MGG_12345/Symbol
            re.compile(
                r"\b([A-Za-z][A-Za-z0-9-]{1,11})/(MGG[_-]?\d{5})(?:[T_\-]\d+)?\b",
                re.IGNORECASE,
            ),
            # Symbol, encoded by MGG_12345
            re.compile(
                r"\b([A-Za-z][A-Za-z0-9-]{1,11}),?\s+(?:encoded by|designated|named)\s+(?:gene\s+)?(MGG[_-]?\d{5})(?:[T_\-]\d+)?\b",
                re.IGNORECASE,
            ),
        ]

        discovered: Dict[str, List[str]] = {}

        def _register(sym: str, mgg_raw: str) -> None:
            canonical_mgg = self.normalize_token(mgg_raw).upper()
            if canonical_mgg not in self.exact_mgg:
                return
            if not self._is_plausible_gene_symbol(sym):
                return
            bucket = discovered.setdefault(canonical_mgg, [])
            for cand in (sym, self.normalize_token(sym)):
                if cand and len(cand) >= 2 and cand.upper() not in INLINE_NON_GENE_WORDS:
                    if cand not in bucket:
                        bucket.append(cand)

        for m in patterns[0].finditer(text):
            _register(m.group(1), m.group(2))
        for m in patterns[1].finditer(text):
            _register(m.group(2), m.group(1))
        for m in patterns[2].finditer(text):
            _register(m.group(1), m.group(2))
        for m in patterns[3].finditer(text):
            _register(m.group(1), m.group(2))

        return discovered

    def find_mentions(
        self,
        text: str,
        extra_aliases: Optional[Dict[str, List[str]]] = None
    ) -> List[GeneMention]:
        """
        Scan text for all gene and mutant mentions.
        Optionally accepts paper-local `extra_aliases` ({MGG_#####: [symbols...]})
        from `extract_inline_aliases()`.
        Returns deduplicated, span-ordered GeneMention objects with canonical MGG IDs.
        """
        if not text or not self.compiled_regex:
            return []

        mentions: List[GeneMention] = []
        seen_spans: Set[Tuple[int, int]] = set()

        def _record_match(raw_match: str, start: int, end: int, forced_id: Optional[str] = None, forced_name: Optional[str] = None) -> None:
            norm_token = self.normalize_token(raw_match).upper()
            if norm_token in BIOLOGICAL_STOPWORDS or raw_match.upper() in BIOLOGICAL_STOPWORDS:
                return

            canonical_id = forced_id or self.resolve(raw_match, context=text)
            if not canonical_id:
                return

            is_mut = bool(re.search(r"[Δ∆]|delta|null|mutant|ko", raw_match, re.IGNORECASE))
            canonical_info = self.exact_mgg.get(canonical_id, {})
            canonical_name = canonical_info.get("gene_name") or forced_name

            if raw_match.upper().startswith("MGG"):
                conf = 1.0
            elif is_mut or re.match(r"^(?:Mo|Po)", raw_match, re.IGNORECASE):
                conf = 0.95
            else:
                conf = 0.90

            span = (start, end)
            if span not in seen_spans:
                seen_spans.add(span)
                mentions.append(GeneMention(
                    raw_text=raw_match,
                    start_char=start,
                    end_char=end,
                    canonical_id=canonical_id,
                    canonical_name=canonical_name,
                    confidence=conf,
                    is_mutant=is_mut
                ))

        for match in self.compiled_regex.finditer(text):
            _record_match(match.group(), *match.span())

        if extra_aliases:
            local_sym_to_id: Dict[str, Tuple[str, str]] = {}
            for mgg_id, syms in extra_aliases.items():
                for s in syms:
                    local_sym_to_id[self.normalize_token(s).upper()] = (mgg_id, syms[0])
            sorted_local = sorted(
                {s for syms in extra_aliases.values() for s in syms},
                key=len,
                reverse=True,
            )
            if sorted_local:
                esc_local = "|".join(re.escape(s) for s in sorted_local)
                local_pat = re.compile(
                    rf"(?<!\w)(?:[Δ∆]|delta[-_]?)?(?:Mo|Po|MO|PO)?(?:[-_])?(?:{esc_local})(?:[Δ∆]|[-_](?:null|mutant|ko))?(?!\w)",
                    re.IGNORECASE,
                )
                for match in local_pat.finditer(text):
                    raw_m = match.group()
                    norm_m = self.normalize_token(raw_m).upper()
                    if norm_m in local_sym_to_id:
                        cid, cname = local_sym_to_id[norm_m]
                        _record_match(raw_m, *match.span(), forced_id=cid, forced_name=cname)

        # Sort mentions by appearance order
        mentions.sort(key=lambda m: m.start_char)
        return mentions

    def get_aliases_for_gene(self, gene_id: str) -> List[str]:
        """
        Return comprehensive list of search terms and aliases for a gene
        (useful for Europe PMC queries and window snippet extraction).
        """
        canonical_mgg = self.normalize_token(gene_id).upper()
        aliases = [canonical_mgg]

        # Add legacy IDs
        for old_id, new_id in LEGACY_MGG_MAP.items():
            if new_id == canonical_mgg and old_id not in aliases:
                aliases.append(old_id)

        symbols = self.mgg_to_symbols.get(canonical_mgg, [])
        for sym in symbols:
            if sym not in aliases:
                aliases.append(sym)
            mo_sym = f"Mo{sym}"
            if mo_sym not in aliases:
                aliases.append(mo_sym)

        return aliases


def sync_uniprot_aliases(db: DatabaseManager, timeout: int = 45) -> Dict[str, int]:
    """
    Fetch M. oryzae gene names & synonyms mapped to MGG_##### ordered locus names
    from UniProt (taxid:242507 and taxid:318829) and populate `gene_aliases`.
    """
    import requests

    pairs: List[Tuple[str, str]] = []
    seen: Set[Tuple[str, str]] = set()
    rows_fetched = 0

    for taxid in ("242507", "318829"):
        url = "https://rest.uniprot.org/uniprotkb/stream"
        params = {
            "query": f"organism_id:{taxid} AND (gene:*)",
            "format": "tsv",
            "fields": "gene_names,gene_primary,gene_synonym,gene_oln",
        }
        try:
            resp = requests.get(url, params=params, timeout=timeout)
            resp.raise_for_status()
        except Exception:
            continue

        lines = resp.text.splitlines()
        for line in lines[1:]:
            rows_fetched += 1
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            gene_names_all, primary, synonyms, oln_str = parts[0], parts[1], parts[2], parts[3]
            combined_names = f"{gene_names_all} {oln_str}".strip()
            mgg_matches = re.findall(r"MGG[_-]?\d{5}", combined_names, re.IGNORECASE)
            if not mgg_matches:
                continue
            # Also ignore trivial ORF placeholders like MGG_01234 or MGCH7_...
            raw_syms = re.split(r"[\s,;]+", f"{primary} {synonyms} {gene_names_all}".strip())
            valid_syms = []
            for sym in raw_syms:
                s = sym.strip()
                if (
                    len(s) >= 3
                    and not re.match(r"^(?:MGG|MGCH\d*)[_-]?\d*", s, re.IGNORECASE)
                    and not re.match(r"^[A-Z0-9]{1,2}\d{3,}$", s)
                    and EntityLinker._is_plausible_gene_symbol(s)
                ):
                    if s not in valid_syms:
                        valid_syms.append(s)

            for mgg_raw in mgg_matches:
                mgg_norm = "MGG_" + re.search(r"\d{5}", mgg_raw).group(0)
                for s in valid_syms:
                    key = (mgg_norm, s)
                    if key not in seen:
                        seen.add(key)
                        pairs.append(key)

    inserted = db.add_gene_aliases_bulk(pairs, source="uniprot")
    return {
        "uniprot_rows_scanned": rows_fetched,
        "candidate_alias_pairs": len(pairs),
        "aliases_inserted": inserted,
    }
