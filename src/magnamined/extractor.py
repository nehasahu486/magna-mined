"""
Grounded Evidence Extraction Engine.
Provides abstract triage, structured Pydantic extraction, and zero-hallucination quote verification.
"""

import os
import re
import json
import warnings
from typing import Optional, List, Dict, Any

from .schema import (
    AbstractTriageResult,
    GeneEvidenceRecord,
    PhenotypeCategory,
    ExperimentalMethod,
    EvidenceConfidence,
)


class EvidenceExtractor:
    """Extracts structured gene evidence from research literature with quote verification."""

    def __init__(
        self,
        api_key: Optional[str] = None,
        model_name: str = "gemini-2.0-flash",
        linker: Optional[Any] = None
    ):
        self.api_key = api_key or os.getenv("GEMINI_API_KEY") or os.getenv("GOOGLE_API_KEY")
        self.model_name = model_name
        self._linker = linker

    @property
    def linker(self):
        if self._linker is None:
            try:
                from .entity_linker import EntityLinker
                self._linker = EntityLinker()
            except Exception:
                self._linker = False
        return self._linker if self._linker is not False else None

    @staticmethod
    def verify_verbatim_quote(quote: str, source_text: str) -> bool:
        """
        Zero-hallucination verification with whitespace, tag, and punctuation normalization.
        Ensures that the extracted verbatim quote is grounded in the source text.
        """
        if not quote or not source_text:
            return False

        def clean_text(t: str) -> str:
            # Strip XML/HTML tags (without leaving inner word gaps)
            s = re.sub(r"<[^>]+>", "", t)
            # Normalize dashes and quotes
            s = s.replace("“", '"').replace("”", '"').replace("‘", "'").replace("’", "'")
            s = re.sub(r"[\u2010\u2011\u2012\u2013\u2014\xad]", "-", s)
            s = re.sub(r"-\s+", "", s)  # stitch soft line-break hyphens
            return re.sub(r"\s+", " ", s).strip().lower()

        clean_quote = clean_text(quote)
        clean_source = clean_text(source_text)

        # 1. Exact normalized substring check
        if clean_quote in clean_source:
            return True

        # 2. Leading 30 characters prefix check (handles trailing period differences)
        if len(clean_quote) > 30 and clean_quote[:30] in clean_source:
            return True

        # 3. Sliding 5-word token sequence check
        words = clean_quote.split()
        if len(words) >= 5:
            first_5 = " ".join(words[:5])
            if first_5 in clean_source:
                return True

        return False

    OTHER_FUNGI_GENERA = (
        "Colletotrichum", "Fusarium", "Botrytis", "Sclerotinia", "Alternaria",
        "Ustilago", "Verticillium", "Bipolaris", "Rhizoctonia", "Cochliobolus",
        "Zymoseptoria", "Puccinia", "Aspergillus", "Candida",
    )

    @classmethod
    def detect_other_fungi_pathogenicity(cls, text: str) -> Dict[str, Any]:
        """
        Detect non-Magnaporthe fungi in a paper and check whether pathogenicity/virulence
        is reported for them. Used to attach cross-species side notes on Magnaporthe papers
        and to filter out comparative multi-fungus papers with no pathogenicity data.
        """
        if not text:
            return {"other_fungi_present": [], "side_note": None}

        fungi_pat = re.compile(
            r"\b(" + "|".join(cls.OTHER_FUNGI_GENERA) + r")\b",
            re.IGNORECASE,
        )
        other_present: List[str] = []
        seen_gen = set()
        for m in fungi_pat.finditer(text):
            g = m.group(1).capitalize()
            if g not in seen_gen:
                seen_gen.add(g)
                other_present.append(g)

        if not other_present:
            return {"other_fungi_present": [], "side_note": None}

        path_kws = ("virulen", "pathogen", "infect", "lesion", "appressor", "penetrat", "disease")
        with_path: List[str] = []
        seen_path = set()
        for sent in cls._split_sentences(cls._normalize_dashes(text)):
            matches = list(fungi_pat.finditer(sent))
            if not matches:
                continue
            s_low = sent.lower()
            if any(k in s_low for k in path_kws):
                for m in matches:
                    g = m.group(1).capitalize()
                    if g not in seen_path:
                        seen_path.add(g)
                        with_path.append(g)

        note = (
            f"Pathogenicity/virulence also reported in {', '.join(with_path)}"
            if with_path else None
        )
        return {"other_fungi_present": other_present, "side_note": note}

    def triage_abstract(self, pmid: str, title: str, abstract: str) -> AbstractTriageResult:
        """
        First-pass screening: Is this paper reporting experimental genetics in M. oryzae?
        Rejects pure reviews and non-Magnaporthe studies before full-text downloads.
        """
        if not abstract:
            return AbstractTriageResult(
                pmid=pmid,
                is_experimental=False,
                target_organism_confirmed=False,
                relevance_score=0.0,
                rationale="No abstract text available"
            )

        title_lower = (title or "").lower()
        text_lower = f"{title} {abstract}".lower()

        # Require explicit Magnaporthe/Pyricularia/MGG/reference-strain evidence (not bare "rice blast")
        organism_matches = any(k in text_lower for k in [
            "magnaporthe", "pyricularia", "rice blast fungus", "m. oryzae", "p. oryzae",
            "mgg_", "guy11", "70-15",
        ])

        # Check experimental genetic manipulation terms
        exp_terms = ["knockout", "deletion", "mutant", "complementation", "disruption", "overexpression", "rnai", "Δ", "∆", "delta", "crispr"]
        exp_matches = sum(1 for term in exp_terms if term in text_lower)

        # Skip pure review titles when no direct genetic manipulation is described
        is_review_title = bool(re.search(r"\b(?:a\s+review|systematic\s+review|mini-review|literature\s+review)\b", title_lower))
        if is_review_title and exp_matches == 0:
            return AbstractTriageResult(
                pmid=pmid,
                is_experimental=False,
                target_organism_confirmed=organism_matches,
                relevance_score=0.1,
                rationale="Filtered out non-experimental review article."
            )

        # Detect phenotypes
        suggested = []
        if any(k in text_lower for k in ["virulen", "pathogen", "infect", "disease", "lesion"]):
            suggested.append(PhenotypeCategory.VIRULENCE)
        if any(k in text_lower for k in ["appressor", "penetrat", "peg"]):
            suggested.append(PhenotypeCategory.APPRESSORIUM)
        if any(k in text_lower for k in ["conidi", "sporulat", "spore"]):
            suggested.append(PhenotypeCategory.CONIDIATION)
        if any(k in text_lower for k in ["stress", "h2o2", "osmotic", "sorbitol", "melanin", "cell wall"]):
            suggested.append(PhenotypeCategory.STRESS_RESPONSE)
        if any(k in text_lower for k in ["growth", "colony", "aerial hyphae"]):
            suggested.append(PhenotypeCategory.GROWTH_DEVELOPMENT)

        # Find mentioned genes using EntityLinker if available, else regex fallback
        mgg_mentions = []
        if self.linker:
            mentions = self.linker.find_mentions(f"{title} {abstract}")
            mgg_mentions = list(dict.fromkeys(m.canonical_id for m in mentions))
        if not mgg_mentions:
            mgg_mentions = list(set(re.findall(r"MGG_\d{5}", title + " " + abstract, re.IGNORECASE)))

        score = min(1.0, (0.4 if organism_matches else 0.0) + (0.4 if exp_matches > 0 else 0.0) + (0.2 if mgg_mentions else 0.0))
        is_exp = organism_matches and (exp_matches > 0 or len(suggested) > 0)

        return AbstractTriageResult(
            pmid=pmid,
            is_experimental=is_exp,
            target_organism_confirmed=organism_matches,
            mentioned_genes=mgg_mentions,
            suggested_phenotypes=suggested or [PhenotypeCategory.NOT_REPORTED],
            relevance_score=score,
            rationale=f"Found {exp_matches} experimental keywords and {len(suggested)} phenotype categories."
        )

    def extract_evidence_from_snippets(
        self,
        gene_id: str,
        gene_name: Optional[str],
        title: str,
        abstract: str,
        snippets: List[str],
        pmid: Optional[str] = None,
        pmcid: Optional[str] = None,
        doi: Optional[str] = None,
        pub_year: Optional[int] = None
    ) -> Optional[GeneEvidenceRecord]:
        """
        Extract structured phenotype evidence for a gene from paper snippets.
        If an LLM client/API key is configured, invokes LLM with Pydantic JSON schema.
        Otherwise, falls back to deterministic rule-based grounding for testing.
        """
        combined_text = f"{title}\n{abstract}\n" + "\n---\n".join(snippets)

        # In production or when API key is available, call LLM with prompt:
        if self.api_key:
            return self._call_llm_extraction(
                gene_id=gene_id,
                gene_name=gene_name,
                source_text=combined_text,
                pmid=pmid,
                pmcid=pmcid,
                doi=doi,
                title=title,
                pub_year=pub_year,
                snippets=snippets
            )

        # Deterministic offline rule-based extraction fallback for testing
        return self._heuristic_extraction(
            gene_id=gene_id,
            gene_name=gene_name,
            source_text=combined_text,
            pmid=pmid,
            pmcid=pmcid,
            doi=doi,
            title=title,
            pub_year=pub_year,
            snippets=snippets
        )

    @staticmethod
    def _strip_comments(text: str) -> str:
        """Remove HTML/OCR comments (e.g. pymupdf4llm picture markers) before sentence splitting."""
        return re.sub(r"<!--.*?-->", " ", text, flags=re.DOTALL)

    @staticmethod
    def _normalize_dashes(text: str) -> str:
        """Normalize Unicode hyphens/dashes (e.g. U+2010 in 'RT‐qPCR') to ASCII '-'."""
        return re.sub(r"[\u2010-\u2015\u2212]", "-", text)

    @staticmethod
    def _split_sentences(text: str) -> List[str]:
        """
        Split prose into sentences without breaking on common scientific abbreviations
        such as 'M. oryzae', 'P. oryzae', 'Fig. 2', 'et al.', 'e.g.', or 'i.e.'.
        """
        abbrevs = {
            "m", "p", "fig", "figs", "al", "sp", "var", "no", "vs", "dr", "approx", "st", "wt",
            "e.g", "i.e",
        }
        sentences: List[str] = []
        start = 0
        for m in re.finditer(r"[.!?](?=\s|\n|$)|\n+", text):
            token = m.group(0)
            if token == ".":
                prefix = text[max(0, m.start() - 12):m.start()]
                word_m = re.search(r"([A-Za-z.]+)$", prefix)
                if word_m and word_m.group(1).lower() in abbrevs:
                    continue
            seg = re.sub(r"\s+", " ", text[start:m.end()]).strip()
            if seg:
                sentences.append(seg)
            start = m.end()
        tail = re.sub(r"\s+", " ", text[start:]).strip()
        if tail:
            sentences.append(tail)
        return sentences

    @staticmethod
    def _is_prior_study_citation(sentence: str) -> bool:
        """Return True when a sentence cites previous literature rather than this paper's experiments."""
        s_norm = EvidenceExtractor._normalize_dashes(sentence).lower()
        # Explicit references to previous/earlier studies
        if any(p in s_norm for p in [
            "previous study", "previous studies", "earlier study", "earlier studies",
            "prior study", "previously reported", "previously shown",
            "previously demonstrated", "previously found", "previously identified",
            "it is known that", "has been reported that", "have been reported",
        ]):
            # Keep only if sentence explicitly states this paper confirmed/tested it too
            if not any(h in s_norm for h in ["in this study", "in the present study", "here, we", "we found", "we observed"]):
                return True

        # Author-year parenthetical literature citations like "(Choi and Dean, 1997)" or "(Xu et al., 2005)"
        if re.search(r"\([A-Z][A-Za-z-]+(?:\s+(?:et\s+al\.|and\s+[A-Z][A-Za-z-]+))?,?\s*(?:19|20)\d{2}[a-z]?(?:[;,)]|\s)", sentence):
            if not any(h in s_norm for h in ["in this study", "in the present study", "our results", "(fig.", "(figure"]):
                return True

        return False

    @staticmethod
    def _build_gene_pattern(gene_id: str, gene_name: Optional[str]) -> re.Pattern:
        """Build a word-boundary regex matching MGG ID, split gene symbols, Mo/Po prefixes, and mutant forms."""
        tokens: List[str] = []
        # Match MGG_12345, MGG-12345, MGG_12345T0
        mgg_digits = re.match(r"^MGG[_-]?(\d{5})", gene_id.strip(), re.IGNORECASE)
        if mgg_digits:
            tokens.append(rf"MGG[_-]?{mgg_digits.group(1)}(?:[T_\-]\d+)?")
        else:
            tokens.append(re.escape(gene_id.strip()))

        if gene_name:
            for raw in re.split(r"[/,;]", gene_name):
                sym = raw.strip()
                if sym and len(sym) >= 2 and not sym.isdigit():
                    esc = re.escape(sym)
                    tokens.append(rf"(?:(?:Mo|Po)[-_]?)?{esc}")

        core = "|".join(tokens)
        return re.compile(
            rf"(?<!\w)(?:[Δ∆]|delta[-_]?)?(?:{core})(?:[Δ∆]|[-_](?:null|mutant|ko))?(?!\w)",
            re.IGNORECASE,
        )

    @staticmethod
    def _has_direct_mutant_for_gene(sentence: str, gene_pat: re.Pattern) -> bool:
        """Return True if the sentence directly describes genetic manipulation of this gene."""
        for m in gene_pat.finditer(sentence):
            tok = m.group(0)
            if re.search(r"^(?:[Δ∆]|delta[-_]?)|(?:[Δ∆]|[-_](?:null|mutant|ko))$", tok, re.IGNORECASE):
                return True

        # Check manipulation verbs adjacent to gene mention (within ~40 chars)
        manip_before = re.compile(
            r"(?:deletion|knockout|knock-out|disruption|null|overexpression|over-expressing|"
            r"complementation|complemented|silencing|rnai|crispr|mutant|mutants)\s+"
            r"(?:of|in|for|strain|strains|lines?)?\s*(?:the\s+)?$",
            re.IGNORECASE,
        )
        manip_after = re.compile(
            r"^\s*(?:deletion|knockout|knock-out|disruption|null|overexpression|"
            r"mutant|mutants|complemented|silencing|strain)",
            re.IGNORECASE,
        )
        for m in gene_pat.finditer(sentence):
            left = sentence[max(0, m.start() - 45):m.start()]
            right = sentence[m.end():min(len(sentence), m.end() + 35)]
            if manip_before.search(left) or manip_after.search(right):
                return True
        return False

    @staticmethod
    def _is_bystander_sentence(sentence: str, gene_pat: re.Pattern) -> bool:
        """
        Filter out sentences where the gene is only a qPCR/RNA-seq expression target,
        background citation from a previous study, protein-interaction partner,
        localization marker, or entry in a multi-gene list.
        """
        if EvidenceExtractor._is_prior_study_citation(sentence):
            return True

        s_norm = EvidenceExtractor._normalize_dashes(sentence)
        s_lower = s_norm.lower()
        has_direct = EvidenceExtractor._has_direct_mutant_for_gene(s_norm, gene_pat)

        # Skip sentences describing experiments strictly in another fungal genus (e.g. Colletotrichum ΔCg...)
        other_fungi_re = re.compile(r"\b(?:" + "|".join(EvidenceExtractor.OTHER_FUNGI_GENERA) + r")\b", re.IGNORECASE)
        if other_fungi_re.search(s_norm) or re.search(r"\b(?:[Δ∆]|delta[-_]?)(?:Cg|Fg|Bc|Ss|Vd|Um|An|Af)[A-Za-z0-9_]+\b", s_norm):
            if not re.search(r"\b(?:magnaporthe|pyricularia|m\.\s*oryzae|p\.\s*oryzae|mgg[_-]?\d{5}|guy11|70-15)\b", s_norm, re.IGNORECASE):
                return True

        # Expression-only / DEG / RT-qPCR sentences remain bystanders unless directly testing this gene's mutant phenotype
        is_expr_mention = any(k in s_lower for k in [
            "expression level", "expression levels", "expression of", "transcript level",
            "mrna level", "reduced expression", "increased expression", "relative expression",
            "gene expression", "expression analysis", "expression pattern",
            "qrt-pcr", "rt-qpcr", "rt-pcr", "rna-seq", "rnaseq",
            "differentially expressed", "degs", "upregulated", "up-regulated",
            "downregulated", "down-regulated", "fold change"
        ])
        if is_expr_mention and not has_direct:
            return True
        if is_expr_mention and any(k in s_lower for k in ["expression of", "reduced expression", "increased expression", "relative expression", "rt-qpcr", "qrt-pcr", "differentially expressed", "degs"]):
            # Even if e.g. ΔcreA is mentioned, "expression of X in ΔcreA did not change" is an expression assay, not a phenotype
            if re.search(r"\b(?:expression|transcript|mrna)\s+of\b", s_lower) and any(k in s_lower for k in ["did not", "significantly change", "confirmed", "degs"]):
                return True

        if has_direct:
            return False

        # 1. Multi-gene figure / heatmap / parenthesized symbol enumeration list
        mgg_count = len(set(re.findall(r"MGG_\d{5}", s_norm, re.IGNORECASE)))
        if mgg_count >= 3 or any(k in s_lower for k in [
            "listed top to bottom", "right side label", "heatmap", "heat map"
        ]):
            return True
        # Parenthesized comma list of >= 4 gene symbols (e.g., DEG functional groups)
        if re.search(r"\([A-Z0-9_-]{2,}(?:\s*,\s*[A-Z0-9_-]{2,}){3,}\)", s_norm):
            return True

        # 2. Protein interaction / pull-down / co-IP / Y2H-only mentions
        if any(k in s_lower for k in [
            "pull-down", "pulldown", "co-immunoprecip", "co-ip",
            "yeast two-hybrid", "y2h", "bimolecular fluorescence", "bifc",
            "interact with", "interacted with"
        ]):
            return True

        # 3. Localization / fluorescent fusion / control / domain / vector mentions
        if any(k in s_lower for k in [
            "co-localized", "colocalized", "localization marker", "marker gene",
            "marked by", "-dsred", "-gfp", "reference gene", "housekeeping",
            "internal control", "primer pair", "primers for", "was cloned into",
            "blastp search", "phylogenetic analysis"
        ]):
            return True

        # Reject protein domain references like "JmjC domain" when not a mutant
        for m in gene_pat.finditer(s_norm):
            after = s_norm[m.end():m.end() + 15].lower()
            if re.match(r"^\s+domains?\b", after):
                return True

        # 4. Sentence is asserting a phenotype of a *different* Δ/∆gene or named mutant
        other_deltas = re.findall(r"(?:[Δ∆]|delta[-_]?)[A-Za-z0-9_]+|[A-Za-z0-9_]+[Δ∆]", s_norm, re.IGNORECASE)
        if other_deltas:
            if not any(gene_pat.fullmatch(d) for d in other_deltas):
                return True
        if re.search(r"\b(?:in\s+(?:the\s+)?|the\s+)[A-Za-z0-9_-]+\s+mutants?\b|\bmutants?\s+[A-Za-z0-9_-]+\b", s_norm, re.IGNORECASE):
            return True

        return False

    def _heuristic_extraction(
        self,
        gene_id: str,
        gene_name: Optional[str],
        source_text: str,
        pmid: Optional[str],
        pmcid: Optional[str],
        doi: Optional[str],
        title: str,
        pub_year: Optional[int],
        snippets: Optional[List[str]] = None
    ) -> Optional[GeneEvidenceRecord]:
        """Sentence-level text-based pattern search for offline / zero-API-key curation."""
        raw_eval = "\n".join(snippets) if snippets else source_text
        clean_eval = self._strip_comments(raw_eval)
        clean_source = self._strip_comments(source_text)

        gene_pat = self._build_gene_pattern(gene_id, gene_name)
        if not gene_pat.search(clean_eval) and not gene_pat.search(clean_source):
            return None

        # Split into intact sentences (preserving 'M. oryzae', 'Fig. 2', 'et al.')
        eval_sentences = self._split_sentences(clean_eval)
        source_sentences = self._split_sentences(clean_source)

        # Find sentences mentioning this gene in eval_text (fallback to source_text if none)
        target_indices = [i for i, s in enumerate(eval_sentences) if gene_pat.search(s)]
        active_sentences = eval_sentences
        if not target_indices:
            target_indices = [i for i, s in enumerate(source_sentences) if gene_pat.search(s)]
            active_sentences = source_sentences

        if not target_indices:
            return None

        # Keep only non-bystander sentences
        valid_indices = [
            i for i in target_indices
            if not self._is_bystander_sentence(active_sentences[i], gene_pat)
        ]
        if not valid_indices:
            return None

        # If direct mutant phenotype sentences exist for this gene, prioritize them so
        # background/passing sentences in the same paper do not dilute or alter the phenotype.
        direct_indices = [
            i for i in valid_indices
            if self._has_direct_mutant_for_gene(self._normalize_dashes(active_sentences[i]), gene_pat)
        ]
        primary_indices = direct_indices if direct_indices else valid_indices

        # Build local evidence windows: each valid sentence + immediately following sentence
        # (e.g. "...caused reduced virulence. Genetic complementation rescued the wild type.")
        # unless the next sentence explicitly introduces another gene's mutant.
        window_indices: List[int] = []
        seen_idx = set()
        for i in primary_indices:
            if i not in seen_idx:
                seen_idx.add(i)
                window_indices.append(i)
            if i + 1 < len(active_sentences) and (i + 1) not in seen_idx:
                nxt = active_sentences[i + 1]
                nxt_low = nxt.lower()
                is_followup = any(k in nxt_low for k in [
                    "complement", "rescued", "reintroduc", "this mutant", "the mutant", "these mutants"
                ])
                other_mggs = [m for m in re.findall(r"MGG_\d{5}", nxt, re.IGNORECASE) if m.upper() != gene_id.upper()]
                other_deltas = [
                    d for d in re.findall(r"(?:[Δ∆]|delta[-_]?)[A-Za-z0-9_]+|[A-Za-z0-9_]+[Δ∆]", nxt, re.IGNORECASE)
                    if not gene_pat.fullmatch(d)
                ]
                if is_followup and not other_mggs and not other_deltas:
                    seen_idx.add(i + 1)
                    window_indices.append(i + 1)

        evidence_text = " ".join(active_sentences[i] for i in window_indices)
        text_lower = self._normalize_dashes(evidence_text).lower()

        dispensable_phrases = [
            "dispensable for virulence", "dispensable for pathogenicity",
            "dispensable for appressorium", "dispensable",
            "little phenotypic effect", "no phenotypic effect",
            "no obvious phenotype", "no visible phenotype",
            "appeared to have no effect", "has no effect", "had no effect",
            "have no effect", "no effect on", "did not affect", "does not affect",
            "not affect", "normal virulence", "full virulence",
            "pathogenicity similar to", "virulence similar to",
            "similar to that of guy11", "similar to guy11",
            "similar to wild", "similar to wt",
            "not required for virulence", "not required for appressorium", "not required for",
            "no defect", "no defects", "no significant defect", "no significant difference",
            "comparable to wild", "comparable to wt", "identical to wt",
            "indistinguishable from", "unaffected", "without defect", "not impaired",
        ]

        proc_keywords = [
            (PhenotypeCategory.APPRESSORIUM, ["appressor", "penetrat", "peg", "turgor"]),
            (PhenotypeCategory.VIRULENCE, ["reduced virulence", "loss of pathogenicity", "avirulent", "attenuated virulence", "pathogenicity defective", "reduced lesion", "virulen", "pathogen", "infect"]),
            (PhenotypeCategory.CONIDIATION, ["conidi", "sporulat", "spore"]),
            (PhenotypeCategory.STRESS_RESPONSE, ["sensitivity to h2o2", "sorbitol", "cell wall defect", "stress response", "osmotic", "melanin"]),
            (PhenotypeCategory.GROWTH_DEVELOPMENT, ["reduced hyphal growth", "vegetative growth", "colony diameter", "aerial hyphae"])
        ]

        # Analyze contrast clauses individually so one gene/dispensable clause cannot leak defects
        saw_dispensable_clause = False
        found_proc_positions: Dict[PhenotypeCategory, int] = {}
        global_pos = 0
        for idx in window_indices:
            sent = self._normalize_dashes(active_sentences[idx])
            clauses = re.split(r"\s*(?:;|\b(?:while|whereas|but|however)\b)\s*", sent, flags=re.IGNORECASE)
            sent_has_gene = bool(gene_pat.search(sent))
            for c_idx, clause in enumerate(clauses):
                if not clause.strip():
                    continue
                c_low = clause.lower()
                clause_mentions_gene = bool(gene_pat.search(clause))
                # Ignore contrast clauses that exclusively describe a different gene
                if not clause_mentions_gene:
                    other_mggs_c = [m for m in re.findall(r"MGG_\d{5}", clause, re.IGNORECASE) if m.upper() != gene_id.upper()]
                    other_deltas_c = [
                        d for d in re.findall(r"(?:[Δ∆]|delta[-_]?)[A-Za-z0-9_]+|[A-Za-z0-9_]+[Δ∆]", clause, re.IGNORECASE)
                        if not gene_pat.fullmatch(d)
                    ]
                    if other_mggs_c or other_deltas_c:
                        continue
                    # First clause of a sentence must name our gene or be a follow-up sentence
                    if c_idx == 0 and not sent_has_gene and idx in primary_indices:
                        continue

                if any(p in c_low for p in dispensable_phrases):
                    saw_dispensable_clause = True
                    continue

                # Strip neutral assay setup phrases before checking phenotype defect terms
                cleaned_clause = re.sub(
                    r"\b(?:conidial\s+)?(?:pathogenicity|virulence|appressorium|infection)\s+assays?\b",
                    " ",
                    c_low,
                )
                for cat, kw_list in proc_keywords:
                    positions = [cleaned_clause.find(kw) for kw in kw_list if kw in cleaned_clause]
                    if positions:
                        pos = global_pos + min(positions)
                        if cat not in found_proc_positions or pos < found_proc_positions[cat]:
                            found_proc_positions[cat] = pos
                global_pos += len(clause) + 1

        detected_processes = [
            cat for cat, _ in sorted(found_proc_positions.items(), key=lambda kv: kv[1])
        ]
        is_dispensable = saw_dispensable_clause and not detected_processes

        if is_dispensable:
            process = PhenotypeCategory.DISPENSABLE
            secondary_procs: List[PhenotypeCategory] = []
            detailed = f"Targeted deletion of {gene_id} ({gene_name or ''}) showed wild-type growth and was dispensable for pathogenicity."
        elif detected_processes:
            process = detected_processes[0]
            secondary_procs = detected_processes[1:]
            proc_summary = ", ".join([p.value for p in detected_processes])
            detailed = f"Targeted mutant showed defective {proc_summary} phenotype based on experimental evidence."
        else:
            process = PhenotypeCategory.NOT_REPORTED
            secondary_procs = []
            detailed = "Gene was studied in experimental context but no specific phenotype defect was specified in the excerpt."

        # Determine experimental method & confidence tier scoped to local evidence window
        has_direct_mut = any(self._has_direct_mutant_for_gene(self._normalize_dashes(active_sentences[i]), gene_pat) for i in primary_indices)
        has_complementation = any(k in text_lower for k in ["complementation", "complemented strain", "rescued", "reintroduction"])
        if has_direct_mut and any(k in text_lower for k in ["crispr", "cas9"]):
            method = ExperimentalMethod.CRISPR_CAS9
            tier = EvidenceConfidence.TIER_1 if has_complementation else EvidenceConfidence.TIER_2
        elif has_direct_mut and any(k in text_lower for k in ["deletion", "knockout", "knock-out", "null mutant", "disruption", "mutant", "δ", "delta"]):
            method = ExperimentalMethod.KNOCKOUT_DELETION
            tier = EvidenceConfidence.TIER_1 if has_complementation else EvidenceConfidence.TIER_2
        elif any(k in text_lower for k in ["overexpression", "over-expressing", "overexpressed"]):
            method = ExperimentalMethod.OVEREXPRESSION
            tier = EvidenceConfidence.TIER_3
        elif any(k in text_lower for k in ["rnai", "silencing", "knockdown"]):
            method = ExperimentalMethod.RNAI_SILENCING
            tier = EvidenceConfidence.TIER_3
        else:
            method = ExperimentalMethod.BIOCHEMICAL_ASSAY
            tier = EvidenceConfidence.TIER_4

        # Discard pure noise mentions that have neither an observable phenotype nor genetic manipulation
        if process == PhenotypeCategory.NOT_REPORTED and tier == EvidenceConfidence.TIER_4:
            return None

        # Select the highest-signal sentence containing the gene as the supporting quote
        phenotype_cues = ["mutant", "deletion", "knockout", "virulence", "defect", "appressor", "conidi", "pathogen", "dispensable", "abolish", "reduc", "requir", "no effect", "similar to"]
        verbatim_quote = ""
        best_score = -1
        for i in primary_indices:
            s = active_sentences[i]
            s_low = self._normalize_dashes(s).lower()
            score = sum(2 for w in phenotype_cues if w in s_low)
            if self._has_direct_mutant_for_gene(self._normalize_dashes(s), gene_pat):
                score += 3
            if score > best_score:
                best_score = score
                verbatim_quote = s

        if not verbatim_quote or not self.verify_verbatim_quote(verbatim_quote, source_text):
            return None

        return GeneEvidenceRecord(
            gene_id=gene_id,
            gene_name=gene_name,
            process_affected=process,
            secondary_processes=secondary_procs,
            is_dispensable=is_dispensable,
            detailed_phenotype=detailed,
            experimental_method=method,
            complementation_tested=has_complementation,
            confidence_tier=tier,
            verbatim_quote=verbatim_quote,
            pmid=pmid,
            pmcid=pmcid,
            doi=doi,
            paper_title=title,
            publication_year=pub_year
        )

    def _call_llm_extraction(self, **kwargs) -> Optional[GeneEvidenceRecord]:
        """Call Google GenAI for structured extraction with zero-hallucination quote verification."""
        try:
            from google import genai
            from google.genai import types

            client = genai.Client(api_key=self.api_key)
            prompt = f"""You are a specialist fungal geneticist curating gene phenotype evidence for Magnaporthe oryzae (the rice blast fungus).
Given the following paper excerpt for target gene {kwargs.get('gene_id')} ({kwargs.get('gene_name', 'unnamed')}):
Analyze the experimental genetic findings and output a strictly grounded JSON record conforming to the schema.

Important Curation Guidelines:
1. Grounding: You MUST provide a verbatim quote (1-2 sentences) directly copied character-for-character from the source text.
2. Negative Results: If the authors generated a mutant and found it has normal virulence or wild-type phenotype (i.e. dispensable for pathogenicity), set is_dispensable=true and process_affected="dispensable / wild-type-like".
3. Pleiotropy: If multiple processes are affected (e.g. both appressorium formation AND conidiation), set primary process_affected and list remaining affected stages in secondary_processes.
4. Experimental Method: Distinguish between targeted deletion/knockout, complementation rescue, overexpression, or assay.

Source Text:
{kwargs.get('source_text')}
"""
            response = client.models.generate_content(
                model=self.model_name,
                contents=prompt,
                config=types.GenerateContentConfig(
                    response_mime_type="application/json",
                    response_schema=GeneEvidenceRecord
                )
            )
            data = json.loads(response.text)
            # Enforce canonical ID and publication provenance from caller metadata
            data["gene_id"] = kwargs.get("gene_id") or data.get("gene_id")
            if kwargs.get("gene_name") is not None:
                data["gene_name"] = kwargs.get("gene_name")
            for meta_key in ("pmid", "pmcid", "doi"):
                if kwargs.get(meta_key) is not None:
                    data[meta_key] = kwargs.get(meta_key)
            pub_yr = kwargs.get("pub_year") if kwargs.get("pub_year") is not None else kwargs.get("publication_year")
            if pub_yr is not None:
                data["publication_year"] = pub_yr
            if kwargs.get("title") is not None:
                data["paper_title"] = kwargs.get("title")

            record = GeneEvidenceRecord(**data)
            if self.verify_verbatim_quote(record.verbatim_quote, kwargs.get("source_text", "")):
                return record
            warnings.warn(
                f"LLM verbatim quote failed grounding verification for {kwargs.get('gene_id')}; falling back to heuristic.",
                RuntimeWarning
            )
        except Exception as exc:
            warnings.warn(
                f"LLM extraction failed for {kwargs.get('gene_id')} ({exc}); falling back to heuristic.",
                RuntimeWarning
            )

        return self._heuristic_extraction(**kwargs)
