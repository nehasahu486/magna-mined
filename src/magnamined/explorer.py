"""
Interactive Explorer & Export Ecosystem for Magna-Mined (Phase 7).
Provides structured filtering, gene evidence dossiers, and exports to
CSV, TSV, Excel (.xlsx), JSON, and self-contained interactive HTML reports.
"""

import csv
import html
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from .database import DatabaseManager
from .reconciler import EvidenceReconciler


class CatalogExplorer:
    """Query, inspect, and export the Magna-Mined gene catalog and evidence ledger."""

    def __init__(self, db: Optional[DatabaseManager] = None, db_path: str = "magnaporthe.db"):
        self.db = db or DatabaseManager(db_path)

    def search_genes(
        self,
        query: Optional[str] = None,
        process: Optional[str] = None,
        confidence: Optional[str] = None,
        source: Optional[str] = None,
        has_extractions: bool = False,
        newly_classified: bool = False,
        knockout_only: bool = False,
        has_additional_functions: bool = False,
        conflicts_only: bool = False,
        limit: Optional[int] = 50,
    ) -> List[Dict[str, Any]]:
        """Filter genes by text query, biological process, tier, provenance source, knockout status, novelty, or conflicts."""
        clauses: List[str] = []
        params: List[Any] = []

        if query:
            like = f"%{query.strip()}%"
            clauses.append(
                "(UPPER(g.gene) LIKE UPPER(?) OR UPPER(COALESCE(g.gene_name, '')) LIKE UPPER(?) "
                "OR LOWER(COALESCE(g.blast2go_annotation, '')) LIKE LOWER(?) "
                "OR LOWER(COALESCE(g.pfam, '')) LIKE LOWER(?) "
                "OR LOWER(COALESCE(g.detailed_mutant_phenotype, '')) LIKE LOWER(?))"
            )
            params.extend([like, like, like, like, like])

        if process:
            clauses.append("LOWER(COALESCE(g.process_affected, '')) LIKE LOWER(?)")
            params.append(f"%{process.strip()}%")

        if confidence:
            clauses.append(
                "(LOWER(COALESCE(g.consensus_tier, '')) = LOWER(?) OR EXISTS ("
                "SELECT 1 FROM extractions e WHERE e.gene_id = g.id AND LOWER(e.confidence_tier) = LOWER(?)))"
            )
            params.extend([confidence.strip(), confidence.strip()])

        if source:
            clauses.append(
                "EXISTS (SELECT 1 FROM extractions e WHERE e.gene_id = g.id AND LOWER(e.evidence_source) = LOWER(?))"
            )
            params.append(source.strip())

        if has_extractions or newly_classified:
            clauses.append("EXISTS (SELECT 1 FROM extractions e WHERE e.gene_id = g.id)")

        if newly_classified:
            clauses.append("LOWER(TRIM(COALESCE(g.baseline_process, 'not reported'))) IN ('not reported', '')")

        if knockout_only:
            clauses.append("g.knockout_reported = 1")

        if has_additional_functions:
            clauses.append("g.additional_functions IS NOT NULL AND TRIM(g.additional_functions) != ''")

        if conflicts_only:
            clauses.append("g.conflict_status = 'conflict'")

        where_sql = f"WHERE {' AND '.join(clauses)}" if clauses else ""
        limit_sql = f"LIMIT {int(limit)}" if limit is not None else ""

        sql = f"""
            SELECT g.*,
                   (SELECT COUNT(*) FROM extractions e WHERE e.gene_id = g.id) AS extraction_count,
                   CASE
                       WHEN LOWER(TRIM(COALESCE(g.baseline_process, 'not reported'))) IN ('not reported', '')
                        AND EXISTS (SELECT 1 FROM extractions e WHERE e.gene_id = g.id)
                       THEN 1 ELSE 0
                   END AS newly_classified
            FROM genes g
            {where_sql}
            ORDER BY newly_classified DESC, extraction_count DESC, g.gene ASC
            {limit_sql}
        """

        with self.db.get_connection() as conn:
            rows = conn.execute(sql, params).fetchall()
            return [dict(r) for r in rows]

    def get_gene_dossier(self, gene_query: str) -> Optional[Dict[str, Any]]:
        """Retrieve a gene's master catalog record, aliases, and ranked extraction quotes."""
        gene = self.db.get_gene(gene_query)
        if not gene:
            return None

        with self.db.get_connection() as conn:
            ext_rows = conn.execute(
                "SELECT * FROM extractions WHERE gene_id = ? ORDER BY id ASC",
                (gene["id"],)
            ).fetchall()

        extractions = sorted(
            [dict(r) for r in ext_rows],
            key=EvidenceReconciler.rank_extraction,
            reverse=True,
        )
        aliases = self.db.linker.get_aliases_for_gene(gene["gene"])

        return {
            "gene": gene,
            "aliases": aliases,
            "extractions": extractions,
        }

    def _fetch_extractions_for_genes(self, gene_ids: List[int]) -> List[Dict[str, Any]]:
        if not gene_ids:
            return []
        placeholders = ",".join("?" for _ in gene_ids)
        sql = f"""
            SELECT g.gene AS mgg_id, g.gene_name, e.*
            FROM extractions e
            JOIN genes g ON e.gene_id = g.id
            WHERE e.gene_id IN ({placeholders})
            ORDER BY g.gene ASC, e.id ASC
        """
        with self.db.get_connection() as conn:
            return [dict(r) for r in conn.execute(sql, gene_ids).fetchall()]

    def export_catalog(
        self,
        output_path: str,
        fmt: Optional[str] = None,
        query: Optional[str] = None,
        process: Optional[str] = None,
        confidence: Optional[str] = None,
        source: Optional[str] = None,
        has_extractions: bool = False,
        newly_classified: bool = False,
        knockout_only: bool = False,
        has_additional_functions: bool = False,
        conflicts_only: bool = False,
    ) -> int:
        """
        Export filtered (or full) catalog to csv, tsv, xlsx, json, or html.
        Returns the number of gene rows exported.
        """
        out = Path(output_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        resolved_fmt = (fmt or out.suffix.lstrip(".")).lower()
        if resolved_fmt not in ("csv", "tsv", "xlsx", "json", "html"):
            raise ValueError(f"Unsupported export format: {resolved_fmt}")

        genes = self.search_genes(
            query=query,
            process=process,
            confidence=confidence,
            source=source,
            has_extractions=has_extractions,
            newly_classified=newly_classified,
            knockout_only=knockout_only,
            has_additional_functions=has_additional_functions,
            conflicts_only=conflicts_only,
            limit=None,
        )
        extractions = self._fetch_extractions_for_genes([g["id"] for g in genes])

        if resolved_fmt in ("csv", "tsv"):
            delim = "\t" if resolved_fmt == "tsv" else ","
            cols = [
                "gene", "gene_name", "baseline_process", "process_affected", "newly_classified",
                "knockout_reported", "knockout_phenotype",
                "additional_functions", "additional_function_references",
                "consensus_tier", "conflict_status", "extraction_count", "size_aa",
                "pmid_or_link_to_publication", "detailed_mutant_phenotype",
                "blast2go_annotation", "pfam",
            ]
            with open(out, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=cols, delimiter=delim, extrasaction="ignore")
                writer.writeheader()
                writer.writerows(genes)
            return len(genes)

        if resolved_fmt == "json":
            by_gid: Dict[int, List[Dict[str, Any]]] = {}
            for ext in extractions:
                by_gid.setdefault(ext["gene_id"], []).append(ext)
            payload = []
            for g in genes:
                item = dict(g)
                item["extractions"] = by_gid.get(g["id"], [])
                payload.append(item)
            out.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
            return len(genes)

        if resolved_fmt == "xlsx":
            import openpyxl

            wb = openpyxl.Workbook()
            ws_genes = wb.active
            ws_genes.title = "Genes"
            gene_cols = [
                "gene", "gene_name", "baseline_process", "process_affected", "newly_classified",
                "knockout_reported", "knockout_phenotype",
                "additional_functions", "additional_function_references",
                "consensus_tier", "conflict_status", "extraction_count", "size_aa",
                "pmid_or_link_to_publication", "detailed_mutant_phenotype",
                "blast2go_annotation", "pfam",
            ]
            ws_genes.append(gene_cols)
            for g in genes:
                ws_genes.append([g.get(c) for c in gene_cols])

            ws_ext = wb.create_sheet("Extractions")
            ext_cols = [
                "mgg_id", "gene_name", "pmid", "pmcid", "doi", "publication_year",
                "process_affected", "experimental_method", "complementation_tested",
                "confidence_tier", "evidence_source", "reviewed_by_human",
                "detailed_phenotype", "verbatim_quote", "paper_title",
            ]
            ws_ext.append(ext_cols)
            for ext in extractions:
                ws_ext.append([ext.get(c) for c in ext_cols])

            wb.save(str(out))
            return len(genes)

        # html interactive report
        return self._write_html_report(out, genes, extractions)

    def _write_html_report(
        self,
        out: Path,
        genes: List[Dict[str, Any]],
        extractions: List[Dict[str, Any]],
    ) -> int:
        by_gid: Dict[int, List[Dict[str, Any]]] = {}
        for ext in extractions:
            by_gid.setdefault(ext["gene_id"], []).append(ext)

        records_for_js = []
        for g in genes:
            records_for_js.append({
                "gene": g["gene"],
                "gene_name": g.get("gene_name") or "",
                "baseline_process": g.get("baseline_process") or "Not reported",
                "newly_classified": bool(g.get("newly_classified")),
                "process_affected": g.get("process_affected") or "not reported",
                "knockout_reported": bool(g.get("knockout_reported")),
                "knockout_phenotype": g.get("knockout_phenotype") or "",
                "additional_functions": g.get("additional_functions") or "",
                "additional_function_references": g.get("additional_function_references") or "",
                "consensus_tier": g.get("consensus_tier") or "",
                "conflict_status": g.get("conflict_status") or "none",
                "blast2go": g.get("blast2go_annotation") or "",
                "pfam": g.get("pfam") or "",
                "extractions": [
                    {
                        "pmid": e.get("pmid") or "",
                        "year": e.get("publication_year") or "",
                        "process": e.get("process_affected") or "",
                        "method": e.get("experimental_method") or "",
                        "tier": e.get("confidence_tier") or "",
                        "source": e.get("evidence_source") or "",
                        "quote": e.get("verbatim_quote") or "",
                        "title": e.get("paper_title") or "",
                    }
                    for e in by_gid.get(g["id"], [])
                ],
            })

        data_json = json.dumps(records_for_js).replace("</", "<\\/")
        stats = self.db.get_stats()
        total_genes = stats.get("total_genes", 0)
        baseline_unclassified = stats.get("baseline_unclassified", 0)
        prior_classified = max(0, total_genes - baseline_unclassified)
        newly_classified_db = stats.get("newly_classified_from_unclassified", 0)
        remaining_unclassified = max(0, baseline_unclassified - newly_classified_db)
        total_classified_now = prior_classified + newly_classified_db
        add_funcs_db = stats.get("genes_with_additional_functions", 0)
        ko_genes_db = stats.get("knockout_reported_genes", 0)
        papers_checked_db = stats.get("papers_triaged", 0)
        exp_papers_db = stats.get("experimental_papers", 0)
        db_summary_json = json.dumps({
            "total_genes": total_genes,
            "prior_classified": prior_classified,
            "baseline_unclassified": baseline_unclassified,
            "newly_classified": newly_classified_db,
            "remaining_unclassified": remaining_unclassified,
            "total_classified_now": total_classified_now,
            "additional_functions_genes": add_funcs_db,
            "knockout_reported_genes": ko_genes_db,
            "papers_checked": papers_checked_db,
            "experimental_papers": exp_papers_db,
        })

        html_doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Magna-Mined Curation Explorer</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; margin: 24px auto; max-width: 1200px; padding: 0 16px; color: #1f2937; background: #f8fafc; }}
  h1 {{ margin-bottom: 4px; color: #0f172a; }}
  .subtitle {{ color: #475569; margin: 0 0 18px 0; font-size: 14px; }}
  .stats {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(175px, 1fr)); gap: 12px; margin: 16px 0; }}
  .stat {{ background: #fff; border: 1px solid #e2e8f0; border-radius: 10px; padding: 12px 16px; box-shadow: 0 1px 2px rgba(15, 23, 42, 0.04); }}
  .stat span {{ font-size: 12px; color: #64748b; display: block; }}
  .stat strong {{ display: block; font-size: 1.4rem; color: #0f766e; margin-top: 4px; }}
  .stat small {{ font-size: 11px; color: #64748b; }}
  .charts-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(340px, 1fr)); gap: 16px; margin: 20px 0; }}
  .chart-card {{ background: #fff; border: 1px solid #e2e8f0; border-radius: 10px; padding: 16px; box-shadow: 0 1px 2px rgba(15, 23, 42, 0.04); }}
  .chart-card h3 {{ margin: 0 0 6px 0; font-size: 15px; color: #0f172a; }}
  .chart-card p.chart-note {{ margin: 0 0 12px 0; font-size: 12px; color: #64748b; }}
  .legend {{ display: flex; flex-wrap: wrap; gap: 10px; margin-top: 10px; font-size: 12px; color: #334155; }}
  .legend-item {{ display: inline-flex; align-items: center; gap: 6px; }}
  .swatch {{ width: 11px; height: 11px; border-radius: 3px; display: inline-block; }}
  .controls {{ display: flex; gap: 10px; margin: 16px 0; flex-wrap: wrap; }}
  input, select {{ padding: 8px 10px; border: 1px solid #cbd5e1; border-radius: 6px; font-size: 14px; background: #fff; }}
  input {{ flex: 1; min-width: 240px; }}
  table {{ width: 100%; border-collapse: collapse; background: #fff; border: 1px solid #e2e8f0; border-radius: 10px; overflow: hidden; }}
  th, td {{ text-align: left; padding: 10px 12px; border-bottom: 1px solid #f1f5f9; vertical-align: top; font-size: 13px; }}
  th {{ background: #f1f5f9; font-weight: 600; color: #334155; }}
  .badge {{ display: inline-block; padding: 2px 7px; border-radius: 999px; font-size: 11px; font-weight: 600; background: #e0f2fe; color: #0369a1; }}
  .tier_1_high {{ background: #dcfce7; color: #166534; }}
  .tier_2_high {{ background: #dbeafe; color: #1e40af; }}
  .tier_3_medium {{ background: #fef9c3; color: #854d0e; }}
  .tier_4_low {{ background: #f3f4f6; color: #4b5563; }}
  .novel {{ background: #d1fae5; color: #065f46; }}
  .conflict {{ background: #fee2e2; color: #991b1b; }}
  .quote-box {{ margin-top: 6px; padding: 8px; background: #f8fafc; border-left: 3px solid #0f766e; font-size: 12px; }}
</style>
</head>
<body>
  <h1>Magna-Mined Gene Curation Explorer</h1>
  <p class="subtitle">Automated <em>Magnaporthe oryzae</em> literature mining dashboard for MagnaGenes v2.0 maintenance.</p>
  <div class="stats">
    <div class="stat"><span>Prior to Search (Baseline Classified)</span><strong>{prior_classified:,}</strong><small>out of {total_genes:,} total genes ({baseline_unclassified:,} Not reported)</small></div>
    <div class="stat"><span>Newly Classified (Total DB)</span><strong>+{newly_classified_db:,}</strong><small>promoted from Not reported ({total_classified_now:,} now classified)</small></div>
    <div class="stat"><span>Baseline Genes Enriched</span><strong>+{add_funcs_db:,}</strong><small>additional biological functions found</small></div>
    <div class="stat"><span>Genes in This View</span><strong>{len(genes):,}</strong><small>{len(extractions):,} extractions | {ko_genes_db:,} KO genes in DB</small></div>
    <div class="stat"><span>Papers Checked</span><strong>{papers_checked_db:,}</strong><small>{exp_papers_db:,} experimental | {stats.get('conflicts_flagged', 0)} conflicts</small></div>
  </div>

  <div class="charts-grid">
    <div class="chart-card">
      <h3>Database Coverage: Before vs. After Literature Search</h3>
      <p class="chart-note">Compares genes already classified in MagnaGenes prior to search against newly classified genes and remaining unclassified genes.</p>
      <div id="chartCoverage"></div>
    </div>
    <div class="chart-card">
      <h3>Curated Biological Processes</h3>
      <p class="chart-note">Distribution of phenotypes among genes matching the active filters below.</p>
      <div id="chartProcesses"></div>
    </div>
    <div class="chart-card">
      <h3>Evidence Tiers &amp; Document Sources</h3>
      <p class="chart-note">Experimental rigor (Tier 1 = KO + complementation) and text source across filtered extractions.</p>
      <div id="chartTiers"></div>
    </div>
    <div class="chart-card">
      <h3>Extractions by Year &amp; Knockout Status</h3>
      <p class="chart-note">Publication years of supporting papers and proportion of filtered genes with reported knockouts.</p>
      <div id="chartYears"></div>
    </div>
  </div>

  <div class="controls">
    <input id="q" type="text" placeholder="Search MGG ID, gene name, annotation, or quote...">
    <select id="tier">
      <option value="">All Confidence Tiers</option>
      <option value="tier_1_high">tier_1_high</option>
      <option value="tier_2_high">tier_2_high</option>
      <option value="tier_3_medium">tier_3_medium</option>
      <option value="tier_4_low">tier_4_low</option>
    </select>
    <select id="onlyExt">
      <option value="all">All Exported Rows</option>
      <option value="novel">Newly Classified Only (Unclassified in Baseline)</option>
      <option value="ko">Knockout Reported Only</option>
      <option value="additional">Additional Functions Discovered Only</option>
      <option value="ext">With Automated Extractions Only</option>
      <option value="conflict">Conflicts Only</option>
    </select>
  </div>
  <table>
    <thead><tr><th>Gene</th><th>Baseline &amp; Status</th><th>Curated Process &amp; Additional Functions</th><th>Knockout Reported &amp; Phenotype</th><th>Tier</th><th>Evidence Quotes</th></tr></thead>
    <tbody id="tbody"></tbody>
  </table>
<script>
const DATA = {data_json};
const DB_STATS = {db_summary_json};
const qEl = document.getElementById('q'), tierEl = document.getElementById('tier'), extEl = document.getElementById('onlyExt'), tb = document.getElementById('tbody');
const esc = s => String(s || '').replace(/[&<>"]/g, c => ({{'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}}[c]));

function renderCoverageChart() {{
  const el = document.getElementById('chartCoverage');
  const total = Math.max(1, DB_STATS.total_genes);
  const segs = [
    {{ label: 'Prior Classified (Baseline)', val: DB_STATS.prior_classified, color: '#0284c7' }},
    {{ label: 'Newly Classified (Magna-Mined)', val: DB_STATS.newly_classified, color: '#10b981' }},
    {{ label: 'Remaining Not Reported', val: DB_STATS.remaining_unclassified, color: '#cbd5e1' }}
  ];
  const r = 46, cx = 64, cy = 64, circ = 2 * Math.PI * r;
  let offset = 0;
  const arcs = segs.map(s => {{
    const len = (s.val / total) * circ;
    const dash = `${{len.toFixed(2)}} ${{circ.toFixed(2)}}`;
    const elStr = `<circle cx="${{cx}}" cy="${{cy}}" r="${{r}}" fill="none" stroke="${{s.color}}" stroke-width="20" stroke-dasharray="${{dash}}" stroke-dashoffset="${{(-offset).toFixed(2)}}"><title>${{esc(s.label)}}: ${{s.val}}</title></circle>`;
    offset += len;
    return elStr;
  }}).join('');
  const pctNow = ((DB_STATS.total_classified_now / total) * 100).toFixed(1);
  const priorW = ((DB_STATS.prior_classified / total) * 100).toFixed(2);
  const newW = ((DB_STATS.newly_classified / total) * 100).toFixed(2);
  const remW = ((DB_STATS.remaining_unclassified / total) * 100).toFixed(2);

  el.innerHTML = `
    <div style="display:flex;gap:16px;align-items:center;flex-wrap:wrap;">
      <svg width="128" height="128" viewBox="0 0 128 128">
        <g transform="rotate(-90 64 64)">${{arcs}}</g>
        <text x="64" y="61" text-anchor="middle" font-size="14" font-weight="700" fill="#0f172a">${{pctNow}}%</text>
        <text x="64" y="76" text-anchor="middle" font-size="10" fill="#64748b">classified</text>
      </svg>
      <div style="flex:1;min-width:170px;font-size:12px;">
        <div style="margin-bottom:8px;">
          <div style="font-weight:600;color:#334155;">Prior to Search vs. Current</div>
          <div style="color:#475569;">Before: <strong>${{DB_STATS.prior_classified.toLocaleString()}}</strong> classified / ${{DB_STATS.baseline_unclassified.toLocaleString()}} Not reported</div>
          <div style="color:#065f46;">Now: <strong>${{DB_STATS.total_classified_now.toLocaleString()}}</strong> classified (+${{DB_STATS.newly_classified}} new)</div>
          <div style="color:#1e40af;">Enriched baseline genes (+additional role): <strong>${{DB_STATS.additional_functions_genes}}</strong></div>
        </div>
        <div style="height:10px;border-radius:999px;overflow:hidden;display:flex;background:#e2e8f0;">
          <div style="width:${{priorW}}%;background:#0284c7;" title="Prior Baseline: ${{DB_STATS.prior_classified}}"></div>
          <div style="width:${{newW}}%;background:#10b981;" title="Newly Classified: ${{DB_STATS.newly_classified}}"></div>
          <div style="width:${{remW}}%;background:#cbd5e1;" title="Remaining Not Reported: ${{DB_STATS.remaining_unclassified}}"></div>
        </div>
      </div>
    </div>
    <div class="legend">${{segs.map(s => `<span class="legend-item"><span class="swatch" style="background:${{s.color}}"></span>${{esc(s.label)}}: <strong>${{s.val.toLocaleString()}}</strong></span>`).join('')}}</div>`;
}}

function hBars(containerId, items, defaultColor) {{
  const el = document.getElementById(containerId);
  if (!items.length) {{
    el.innerHTML = '<div style="color:#94a3b8;font-size:12px;padding:12px 0;">No data in current filter view.</div>';
    return;
  }}
  const maxVal = Math.max(1, ...items.map(i => i.val));
  const rowH = 26, w = 340, labelW = 148, barMaxW = 140;
  const svgH = items.length * rowH + 6;
  const rowsSvg = items.map((it, idx) => {{
    const y = idx * rowH + 6;
    const bw = Math.max(2, (it.val / maxVal) * barMaxW);
    return `
      <text x="0" y="${{y + 13}}" font-size="11" fill="#334155">${{esc(it.label.slice(0, 22))}}</text>
      <rect x="${{labelW}}" y="${{y + 2}}" width="${{bw.toFixed(1)}}" height="14" rx="3" fill="${{it.color || defaultColor}}"><title>${{esc(it.label)}}: ${{it.val}}</title></rect>
      <text x="${{labelW + bw + 6}}" y="${{y + 13}}" font-size="11" font-weight="600" fill="#0f172a">${{it.val}}</text>`;
  }}).join('');
  el.innerHTML = `<svg width="100%" viewBox="0 0 ${{w}} ${{svgH}}" style="max-width:${{w}}px;overflow:visible;">${{rowsSvg}}</svg>`;
}}

function renderDynamicCharts(rows) {{
  // 1. Processes
  const procMap = {{}};
  rows.forEach(r => {{
    String(r.process_affected || 'not reported').split(',').forEach(p => {{
      const k = p.trim();
      if (k) procMap[k] = (procMap[k] || 0) + 1;
    }});
  }});
  const procItems = Object.entries(procMap).sort((a,b) => b[1] - a[1]).slice(0, 7).map(([k,v]) => ({{ label: k, val: v, color: '#0f766e' }}));
  hBars('chartProcesses', procItems, '#0f766e');

  // 2. Confidence tiers + sources
  const tierColors = {{ tier_1_high: '#16a34a', tier_2_high: '#2563eb', tier_3_medium: '#ca8a04', tier_4_low: '#64748b' }};
  const tierCounts = {{ tier_1_high: 0, tier_2_high: 0, tier_3_medium: 0, tier_4_low: 0 }};
  const srcCounts = {{ jats_xml: 0, pdf: 0, abstract_only: 0 }};
  const yearCounts = {{}};
  let koYes = 0, koNo = 0;

  rows.forEach(r => {{
    if (r.knockout_reported) koYes++; else koNo++;
    if (r.extractions.length) {{
      r.extractions.forEach(e => {{
        if (e.tier && tierCounts[e.tier] !== undefined) tierCounts[e.tier]++;
        if (e.source && srcCounts[e.source] !== undefined) srcCounts[e.source]++;
        if (e.year) yearCounts[e.year] = (yearCounts[e.year] || 0) + 1;
      }});
    }} else if (r.consensus_tier && tierCounts[r.consensus_tier] !== undefined) {{
      tierCounts[r.consensus_tier]++;
    }}
  }});

  const tierItems = Object.entries(tierCounts).filter(([,v]) => v > 0).map(([k,v]) => ({{ label: k, val: v, color: tierColors[k] }}));
  const srcItems = Object.entries(srcCounts).filter(([,v]) => v > 0).map(([k,v]) => ({{ label: 'src: ' + k, val: v, color: '#0284c7' }}));
  hBars('chartTiers', [...tierItems, ...srcItems], '#0284c7');

  // 3. Publication Years vertical bar chart + KO badge summary
  const yrEl = document.getElementById('chartYears');
  const yrEntries = Object.entries(yearCounts).sort((a,b) => Number(a[0]) - Number(b[0]));
  const koSummary = `<div class="legend" style="margin-top:0;margin-bottom:8px;">
    <span class="legend-item"><span class="swatch" style="background:#16a34a"></span>Knockout Reported: <strong>${{koYes}}</strong></span>
    <span class="legend-item"><span class="swatch" style="background:#94a3b8"></span>No KO Reported: <strong>${{koNo}}</strong></span>
  </div>`;
  if (!yrEntries.length) {{
    yrEl.innerHTML = koSummary + '<div style="color:#94a3b8;font-size:12px;">No dated extractions in current filter view.</div>';
  }} else {{
    const maxY = Math.max(1, ...yrEntries.map(e => e[1]));
    const bw = 34, gap = 14, chartH = 92;
    const svgW = Math.max(240, yrEntries.length * (bw + gap) + 16);
    const bars = yrEntries.map(([yr, cnt], idx) => {{
      const h = Math.max(4, (cnt / maxY) * 68);
      const x = 10 + idx * (bw + gap);
      const y = chartH - 18 - h;
      return `
        <rect x="${{x}}" y="${{y.toFixed(1)}}" width="${{bw}}" height="${{h.toFixed(1)}}" rx="4" fill="#6366f1"><title>${{yr}}: ${{cnt}} extractions</title></rect>
        <text x="${{x + bw/2}}" y="${{(y - 4).toFixed(1)}}" text-anchor="middle" font-size="10" font-weight="600" fill="#1e293b">${{cnt}}</text>
        <text x="${{x + bw/2}}" y="${{chartH - 3}}" text-anchor="middle" font-size="10" fill="#475569">${{esc(yr)}}</text>`;
    }}).join('');
    yrEl.innerHTML = koSummary + `<svg width="100%" viewBox="0 0 ${{svgW}} ${{chartH}}" style="max-width:${{svgW}}px;">${{bars}}</svg>`;
  }}
}}

function render() {{
  const q = qEl.value.trim().toLowerCase(), tier = tierEl.value, mode = extEl.value;
  let rows = DATA.filter(r => {{
    if (tier && r.consensus_tier !== tier && !r.extractions.some(e => e.tier === tier)) return false;
    if (mode === 'novel' && !r.newly_classified) return false;
    if (mode === 'ko' && !r.knockout_reported) return false;
    if (mode === 'additional' && !r.additional_functions) return false;
    if (mode === 'ext' && r.extractions.length === 0) return false;
    if (mode === 'conflict' && r.conflict_status !== 'conflict') return false;
    if (!q) return true;
    return (r.gene + ' ' + r.gene_name + ' ' + r.process_affected + ' ' + r.knockout_phenotype + ' ' + r.additional_functions + ' ' + r.blast2go + ' ' + r.extractions.map(e => e.quote).join(' ')).toLowerCase().includes(q);
  }});
  renderDynamicCharts(rows);
  tb.innerHTML = rows.slice(0, 300).map(r => {{
    const quotes = r.extractions.map(e =>
      `<div class="quote-box"><span class="badge ${{esc(e.tier)}}">${{esc(e.tier)}}</span> <span class="badge">${{esc(e.source)}}</span> ${{e.pmid ? `<a href="https://pubmed.ncbi.nlm.nih.gov/${{esc(e.pmid)}}/" target="_blank">PMID:${{esc(e.pmid)}}</a>` : ''}} (${{esc(e.year)}})<div>${{esc(e.quote)}}</div></div>`
    ).join('');
    const novelBadge = r.newly_classified ? `<span class="badge novel">Newly Classified</span>` : `<span class="badge">In Baseline</span>`;
    const conflictBadge = r.conflict_status === 'conflict' ? `<span class="badge conflict">conflict</span>` : '';
    const koHtml = r.knockout_reported
      ? `<span class="badge tier_1_high">Yes</span><div style="margin-top:4px;font-size:12px;">${{esc(r.knockout_phenotype)}}</div>`
      : `<span class="badge tier_4_low">No</span>`;
    const addHtml = r.additional_functions
      ? `<div style="margin-top:6px;"><span class="badge tier_2_high">+Additional: ${{esc(r.additional_functions)}}</span><div style="font-size:11px;color:#4b5563;margin-top:2px;">${{esc(r.additional_function_references)}}</div></div>`
      : '';
    return `<tr><td><strong>${{esc(r.gene)}}</strong><br><span style="color:#6b7280">${{esc(r.gene_name)}}</span></td><td>${{esc(r.baseline_process)}}<br>${{novelBadge}}</td><td>${{esc(r.process_affected)}}${{addHtml}}</td><td>${{koHtml}}</td><td>${{r.consensus_tier ? `<span class="badge ${{esc(r.consensus_tier)}}">${{esc(r.consensus_tier)}}</span>` : '-'}} ${{conflictBadge}}</td><td>${{quotes || '<span style="color:#9ca3af">Baseline only</span>'}}</td></tr>`;
  }}).join('');
}}
renderCoverageChart();
qEl.addEventListener('input', render); tierEl.addEventListener('change', render); extEl.addEventListener('change', render); render();
</script>
</body>
</html>"""
        out.write_text(html_doc, encoding="utf-8")
        return len(genes)
