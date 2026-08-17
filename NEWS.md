# CENetwork (development)

### Liver network v1.8 becomes the reference

- New datasets `liver_1.8_network` and `liver_1.8_rwr_closest_dfr`. v1.8 is v1.7
  with a pruned pathway layer (2018 → 814 Reactome pathways), keeping only
  pathways that map to hepatotoxicity key characteristics. This removes the
  generic root pathways that acted as diffusion hubs.
- The APAP / VPA use cases (`data-raw/result_use_case_paper.R`), the ORA/GSEA
  comparison (`data-raw/comparison_ora_gsea.R`, dataset `comparison_ora_gsea`)
  and the vignette now all run on v1.8.
- `comparison_ora_gsea`: the diffusion term sets are now thresholded at
  `padj < 0.05` from `compute_diffusion_enrichment()`, instead of listing every
  term present in the subnetwork — the ORA and GSEA sets are filtered at
  `qvalue < 0.05`, so the three methods are now compared at the same threshold.
- The diffusion / ORA / GSEA supplementary table (`results/comparison_ora_gsea.xlsx`)
  is now written by `data-raw/comparison_ora_gsea.R` (it used to sit at the end of
  `data-raw/result_use_case_paper.R`, which required the not-yet-computed results).

### The Venn diagrams are replaced by a term × method dot plot

- `vignettes/img/comparison_venn.png` gives way to `comparison_dotplot.png`. On v1.8
  no panel has a triple intersection and pairwise overlaps are 0-2 terms, so the four
  Venn panels spent a page drawing empty regions. Worse, a Venn cannot separate a term
  a method *tested and rejected* from one it never tested — the dot plot shows that as
  a third state (filled disc / open circle / cross).
- `comparison_ora_gsea` gains a `full` element holding every term each method tested
  within the network background, with its qvalue. These were already computed
  (`pvalueCutoff = 1`) and simply discarded.
- No truncation: every term significant for at least one method is drawn. Ranking by
  q would have been misleading — the 39 GSEA hits on APAP are dominated by a redundant
  cluster of nested Reactome translation sets (11 tie at q = 5.2e-09), and a top-N would
  have hidden the bioactivation pathways (Biological oxidations, Phase I, Xenobiotics,
  Cytochrome P450) that carry the known APAP mechanism.
- "not tested" is reserved for the two methods that can genuinely leave a term unscored:
  the diffusion (node never visited) and GSEA (gene set outside its 10-500 gene range).
  ORA always has a verdict — `enrichGO()` / `enrichPathway()` only report gene sets
  sharing at least one gene with the signature, so a term absent from their output has a
  zero overlap, i.e. p = 1, and is drawn as tested-and-not-significant. On APAP this is
  not a detail: only 23 of 814 pathways share a gene with the 28-gene signature.
- Terms sit in their own framed band (`ggh4x::facet_nested`), which also merges the
  repeated compound x ontology strip. `ggh4x` joins `Suggests`; `ggvenn` and `patchwork`
  leave it — the Venn was their only use.

# CENetwork 1.1.0 (2026-08-05)

### Breaking-free of netOmics

- netOmics is removed from `Imports`. The package is deprecated in Bioconductor
  (removed in 3.21), which would have made CENetwork uninstallable and its
  datasets non-reproducible.
- New exported function `combine_layers()`: merges two layers of a multi-layer
  network, or connects a layer to new nodes through a table of interactions.
  Reimplementation of `netOmics::combine_layers()` and `netOmics:::merge_graphs()`
  (GPL-3), behaviour preserved as-is so networks rebuilt with it are identical
  to the published ones.
- Internal `remove_unconnected_nodes()` replaces `netOmics:::remove_unconnected_nodes()`
  in `RWR_build_complete()`.
- Dropped stale `@importFrom netOmics` tags in `get_route()` (the imported
  functions were not called).

# CENetwork 1.0.2 (2026-05-25)

### New features

- `get_route()` gains `mandatory_nodes` parameter: fixed seed nodes always
  present in the subnetwork and held constant across permutations.
- `get_route()` gains `do_permutation`, `n_permut`, and `seed_permut`
  parameters for permutation-based significance testing.
- Permutation results stored as vertex attributes `visit_count` and `z_score`
  on the output network.
- New function `compute_diffusion_enrichment()`: computes z-scores and
  BH-adjusted p-values for GO and pathway nodes from permutation results.
- New dataset `liver_1.7_network`: updated liver network (v1.7) with extended
  drug layer and refined hepatotoxicity GO term list.
- New dataset `liver_1.7_rwr_closest_dfr`: pre-calculated shortest paths
  matching `liver_1.7_network`.
- New dataset `comparison_ora_gsea`: pre-computed ORA and GSEA enrichment
  results for APAP and VPA (padj < 0.05), filtered to the network-specific
  hepatotoxicity background, for comparison with diffusion results.
- New `data-raw/result_use_case_paper.R`: reproducible script for the APAP/VPA
  use case with permutations.
- New `data-raw/comparison_ora_gsea.R`: reproducible script for ORA/GSEA
  computation and Venn diagram generation.

### Bug fixes

- `report()`: fixed missing `@importFrom` for `summarise`, `n`, `nest`,
  `unnest`, `map_dbl`, `imap_dbl`, `map`, `keep`, and `str_detect`, which
  caused errors when the package was loaded without `tidyverse` attached.
- `report()`: `complete_network` is now optional (default `NULL`); ORA and
  radar plot blocks are skipped gracefully when omitted.

### Vignette

- Updated throughout to use `liver_1.7_network` and `liver_1.7_rwr_closest_dfr`.
- Added permutation analysis and enrichment statistics sections.
- Added Diffusion vs ORA vs GSEA comparison section with 2×2 Venn diagrams
  (GO terms and Reactome pathways, APAP and VPA).

---

2022-06-28

+ First

2022-07-11

+ Fix liver network with no display_name attribute
+ Issues with nodes: P35228, P30711, P29475, P48448, P11511, O15528, P16444, 
    P11086, P16050, Q07973, Q9GZR5 which were uncorrectly assigned to 
    "drug/compound" type (proteins)
+ Add url (link attribute) to GO and side_effect nodes.
+ Add SIDER_se_id (side effect id)
+ Recomputing Closest dfr as node types attributes were changed.
