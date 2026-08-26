# liver_1.8 — réseau foie à couche pathway hépatotoxique

Construit `liver_1.8` = `liver_1.7` dont la **couche pathway** est restreinte aux
pathways Reactome hépatotoxiques, puis en précalcule la diffusion RWR.

## Motivation
Sur `liver_1.7` la diffusion converge sur des pathways racines génériques
(*Metabolism*, *Gene expression*…) car la couche pathway (2018 nœuds) est
sélectionnée par simple expression foie, sans pertinence toxicologique.
`liver_1.8` restreint cette couche aux mécanismes hépatotoxiques.

## Source de la curation (citable)
Mapping **pathway → Key Characteristic** publié et revu par les pairs :
- **Tsai et al. 2025**, *Toxicol Sci*, doi:10.1093/toxsci/kfaf036 — Supplemental
  File 1 = `ressources/toxsci-25-0020-File008.xlsx` (matrice pathway × 34 termes KC,
  assignation par 2 scientifiques + consensus).
- On retient les pathways assignés à l'un des **11 termes umbrella des 12 Key
  Characteristics des hépatotoxicants** (**Rusyn et al. 2021**, *Hepatology*,
  doi:10.1002/hep.31999 ; domaine « HEP »). KC2 et KC3 partagent un terme → 11.

Aucun réglage anti-hub : fidélité au mapping publié (les hubs résiduels sont
laissés au z-score par permutation de la diffusion).

## Scripts (exécuter depuis la racine du package)
1. `01_build_network_1.8.R` — parse File008 → set hépatotox (noms → R-HSA via
   `ressources/ReactomePathways.txt`) → élague `liver_1.7` → `data/liver_1.8_network.rda`.
   Écrit aussi `hepatotox_pathways.tsv` (traçabilité pathway → KC).
   Résultat : couche pathway **2018 → 814**, réseau 31033 → 29829 nœuds.
   `liver_1.8` est exactement le sous-graphe induit de `liver_1.7` (aucun attribut modifié).
2. `02_rwr_1.8.R` — `pipeline_RWR()` (restart = 0.7) → `data/liver_1.8_rwr_closest_dfr.rda`
   (même structure que `liver_1.7_rwr_closest_dfr`). Calcul dense n×n (~29,8k nœuds) → plusieurs Go de RAM.

## Dépendances dev
`readxl` (parse xlsx, non ajouté aux Imports du package), + Imports du package
(`igraph`, `Matrix`, `netOmics`, tidyverse pour le build).

## Ressources
- `toxsci-25-0020-File008.xlsx` — Tsai 2025 Suppl. File 1 (mapping Reactome × KC).
- `ReactomePathways.txt` — Reactome (nom ↔ R-HSA id, Homo sapiens).
