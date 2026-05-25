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
