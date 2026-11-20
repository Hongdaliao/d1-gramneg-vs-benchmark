# d1-gramneg-vs-benchmark

Reproducible code, derived data and figures for:

**Benchmarking structure-based virtual screening across Gram-negative antibacterial targets: target-specific enrichment, co-crystal ligand pose limitations and failure modes**

## Contents
- `code/` — analysis and figure-generation scripts (Track A/B, calibration gate, method-class comparison)
- `manuscript/` — manuscript markdown/PDF (for reference)
- `SI_tables_S1-S31.xlsx` — supporting tables S1–S31
- `results/` — derived per-target metrics, sensitivity analyses, gate results and method-class comparison tables (JSON/CSV)
- `figures/` — final figures (PNG/PDF, 300 dpi)
- `data/` — curated manifests, decoy protocol and ligand-class re-annotation (no raw docking poses)
- `checksums.txt` — SHA-256 of released files

## What is not included
Raw docking trajectories, receptor structure files and third-party model weights are not redistributed. The DrugCLIP model weights are licensed CC BY-NC 4.0 and must be obtained from the original authors.

## Reproduction
Scripts are Python 3.10+ with standard scientific stack (numpy, pandas, matplotlib, scipy). See `environment_note.md` for the versions used.

## Data availability
- Derived tables and code: this repository (archived on Zenodo, DOI: [to be inserted]).
- The MAES calibration manifest/report is identified by SHA-256 `e4ae5c23461474f773a016ec415e9da8e5dd6ff10bbef55eda73c726cb4fa3ca`.

## License
Code: MIT. Derived data and figures: CC BY 4.0. See `LICENSE` and `DATA_LICENSE`.

## Citation
See `CITATION.cff`.
