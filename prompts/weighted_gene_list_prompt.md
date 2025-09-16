You are assisting in building a weighted gene panel representing a signaling pathway for downstream scoring on RNA‑seq (bulk and single‑cell). Produce a CSV with exactly these six columns and a header (no extra columns, no prose):

gene_symbol,gene_name,module,annotation,evidence_source,weight

Instructions:
- gene_symbol: Official gene symbol for the specified species (mouse or human), one per row.
- gene_name: Descriptive full name for the gene.
- module: High‑level role within the pathway (e.g., ligand, receptor, co‑receptor, adaptor, kinase, transcription factor).
- annotation: One short phrase capturing the specific role (e.g., “core: p75NTR receptor”, “extended adaptor: NF‑κB arm”).
- evidence_source: Curated sources that justify inclusion (e.g., “KEGG:mmu04722; Reactome: p75NTR signaling; peer‑reviewed reviews”).
- weight: Real number > 0 encoding relative importance; higher means more influence. Suggested scale: 0.3–1.0 (core nodes ≈1.0; extended ≈0.3–0.7). Keep weights within [0, 1.2].

Coverage guidance:
- Include core ligand(s), receptor(s)/co‑receptors, critical adaptors/kinases, and key transcriptional effectors.
- Prefer 10–60 genes total. If uncertain, include with a lower weight rather than omitting.

Output only CSV rows following the header above.

Provide the panel for: [describe pathway], species: [human|mouse].
