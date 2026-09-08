# Resources Catalog

Everything gathered for: *Network Pharmacology Approaches to Drug Repurposing — Do
Topological Drug-Target Similarity Scores Predict Clinical Repurposing Success?*

**Summary**: 74 papers downloaded (1,192 screened), 12 datasets, 2 repositories cloned,
plus a validated Python 3 implementation of the method under test and a fully assembled
benchmark of 5,797 labelled drug–indication pairs.

---

## Papers

74 open-access PDFs in `papers/` (of the top 110 ranked; 36 paywalled, abstracts retained).
Full-text XML for 5 key papers in `papers/fulltext_xml/`.
Detailed listing with abstracts: `papers/README.md`.

The ones that matter most:

| Rank | Title | Year | Cites | File | Why it matters |
|---|---|---|---|---|---|
| 12 | Network-based in silico drug efficacy screening (Guney et al.) | 2016 | 437 | `papers/012_…pdf` | **Defines the metric under test.** Reports AUC = 0.66; shows shared-target similarity (0.80) ≈ proximity similarity (0.81), P = 0.12 |
| 10 | A standard database for drug repositioning (Brown & Patel) | 2017 | 206 | `papers/010_…pdf` + XML | **Supplies the ground truth**; exists precisely to fix the assumed-negatives problem |
| 18 | Uncovering disease-disease relationships through the incomplete interactome (Menche et al.) | 2015 | 1,089 | abstract only (paywalled) | Disease-module formalism; identifiability conditions under incompleteness |
| 6 | Network-based approach to prediction and population-based validation (Cheng et al.) | 2018 | 392 | `papers/006_…pdf` | Proximity validated against patient claims data |
| 11 | Network medicine framework for COVID-19 repurposing (Morselli Gysi et al.) | 2021 | 345 | `papers/011_…pdf` | Shows network pipelines disagree with each other |
| 19 | Tissue specificity of human disease module (Kitsak et al.) | 2016 | 77 | `papers/019_…pdf` | Basis for the tissue-specificity arm |
| 52 | A reference map of the human binary protein interactome (HuRI) | 2020 | 1,083 | `papers/052_…pdf` | The study-bias-free network that makes the bias test possible |
| 66 | Genomic data integration systematically biases interactome mapping | 2018 | 31 | `papers/066_…pdf` | Network "improvements" can be circular |
| 41 | GhostBuster: literature-unbiased gene prioritization | 2025 | 0 | `papers/041_…pdf` | Literature bias inflates standard metrics |
| 58 | Benchmarking data leakage in biomedical KGE | 2026 | 2 | `papers/058_…pdf` | Leakage inflates biomedical link prediction |
| 70 | BETA: benchmark for computational drug–target prediction | 2022 | 15 | `papers/070_…pdf` | Random-shuffle CV is inadequate; use stratified splits |
| 21 | Clinical-trial and ontology-derived positive/negative benchmarks | 2026 | 1 | `papers/021_…pdf` | The negative-set problem is still open |

## Datasets

12 datasets, ~700 MB. Not committed to git; reproduce with `code/download_datasets.sh`.
Full detail in `datasets/README.md`.

| Name | Source | Size | Role | Location |
|---|---|---|---|---|
| **repoDB** | Brown & Patel 2017 (via authors' GitHub RData) | 10,800 pairs | **Ground truth** — 6,677 approved / 4,123 clinically failed | `datasets/repodb/` |
| BioGRID | thebiogrid.org | 1.29M human physical | Literature-curated interactome | `datasets/ppi/` |
| HuRI / HI-union / Lit-BM | interactome-atlas.org | 52k / 63k / 12k edges | **Study-bias-free** interactomes | `datasets/ppi/` |
| STRING v12 | stringdb-downloads.org | 233k edges @ ≥700 | Weighted multi-evidence network | `datasets/ppi/` |
| Hetionet v1.0 | github.com/hetio/hetionet | 2.25M edges | Drug targets (CbG), PPI (GiG), tissue expression (AeG) | `datasets/hetionet/` |
| DisGeNET | Zenodo 48426 (`dhimmel/disgenet`) | 427k assoc | **Disease genes, CUI-keyed** — joins repoDB directly | `datasets/disgenet/` |
| DrugCentral | drugcentral.org | 19,379 rows | Alternative drug–target source | `datasets/drugcentral/` |
| GTEx v8 | GTEx portal | 56,200 × 54 | Tissue expression | `datasets/tissue/` |
| Jensen TISSUES | download.jensenlab.org | 4.5 MB | Alternative tissue expression | `datasets/tissue/` |
| Jensen DISEASES | download.jensenlab.org | 11 MB | Alternative disease genes (DOID) | `datasets/diseases_jensen/` |
| Disease Ontology | DiseaseOntology GitHub | 7.2 MB | DOID ↔ UMLS CUI xrefs | `datasets/ontology/` |
| NCBI gene mappings | ftp.ncbi.nlm.nih.gov | 303 MB | Symbol/Ensembl → Entrez | `datasets/mappings/` |

**Assembled benchmark** (`datasets/processed/`, built by `code/build_dataset.py`):
7 interactomes × 2 disease-gene definitions × 5,797 max labelled pairs.

## Code repositories

| Name | URL | Purpose | Location | Notes |
|---|---|---|---|---|
| emreg00/toolbox | github.com/emreg00/toolbox | Reference implementation of network proximity, by the method's author | `code/emreg00_toolbox/` | **Python 2 — will not run here.** Used as the spec for `code/proximity.py` |
| repoDB | github.com/adam-sam-brown/repoDB | Source of the gold standard + its build scripts | `code/repoDB/` | `Shiny_Application/data/shiny.RData` holds the full database |

## Notes on resource gathering

### Search strategy
The paper-finder service was unavailable ("service not running at localhost:8000"), so
literature search was done directly against **Europe PMC** (the workhorse for this biomedical
topic), **Semantic Scholar**, and **arXiv** across 30 queries in two rounds. Round 1 covered
the method and its applications; round 2 deliberately targeted the *critique* literature
(study bias, data leakage, evaluation pitfalls, negative-set construction), because the
hypothesis is a methodological challenge and the field's own papers are overwhelmingly
self-confirming. 1,477 hits → 1,192 unique after title normalisation.

### Selection criteria
Ranked by a hypothesis-weighted term score plus citation and recency boosts, with a
**penalty for TCM-style "network pharmacology" papers**, which otherwise flood this query
space and are unrelated to interactome proximity. Top 110 were pursued for full text.

### Challenges encountered
- **paper-finder service down** → manual multi-source API search (documented above).
- **Semantic Scholar rate-limiting** (HTTP 429) throughout; Europe PMC carried the search.
  Roughly a third of S2 queries returned nothing even with exponential backoff.
- **PDF retrieval**: the documented Europe PMC `fullTextPdf` REST route returns 404; the
  working route is `https://europepmc.org/articles/{PMCID}?pdf=render`. OpenAlex was used
  as the DOI → OA resolver (Unpaywall was avoided since it requires submitting an email
  address to a third party).
- **repoDB has no direct download** — it is served by a Shiny app. Extracted from the authors'
  own repository RData via `pyreadr`; row counts match the paper exactly, confirming fidelity.
- **DisGeNET is now gated** behind registration. Used the archived `dhimmel/disgenet` Zenodo
  release (same underlying data, ODbL 1.0).
- **HuRI downloads 301 and drop the body** on the `www.` host; `https://interactome-atlas.org/`
  works.
- 36 of 110 papers are paywalled/publisher-blocked (403). Abstracts retained for all;
  the three most important of these (Menche 2015, Barabási 2011, Huttlin 2017) are
  well-summarised in papers we do have.
- `uv add` fails without `[tool.uv] package = false` (hatchling tries to build a wheel).

### Gaps and workarounds
- **DrugBank full data is licence-gated**; used Hetionet's DrugBank-keyed CbG edges
  (1,389 drugs, joins to repoDB directly) with DrugCentral as a sensitivity check.
- **CTD** (3.2 GB) was skipped as disproportionate; DisGeNET covers the need.
- **AACT/ClinicalTrials.gov** was not ingested directly — repoDB already encodes it, and a
  prospective post-2017 outcome analysis was pruned from scope (see `planning.md`, D8).
- Disease-gene coverage is the binding constraint: 1,234 of 2,051 repoDB CUIs have DisGeNET
  genes, which is what caps the benchmark at 5,797 of 10,800 pairs.

## Recommendations for experiment design

1. **Primary dataset**: repoDB, joined via `datasets/processed/pairs.tsv`. Use the
   `biogrid × disgenet_all` configuration (5,797 pairs) as the headline and
   `× disgenet_curated` (2,993 pairs, higher confidence) as the confirmation.
2. **Method under test**: `code/proximity.py`, `closest` z-score, n_random=1000 — faithful to
   Guney et al. Also report `shortest`, `kernel`, `centre`, and raw distance.
3. **Baselines** (all cheap, all necessary): overlap count, Jaccard, number of drug targets,
   disease module size, mean target degree. The hypothesis's second claim is specifically
   about overlap count, so that comparison needs a **DeLong test**, not just two numbers.
4. **The key control**: rerun the same pipeline with Guney-style *random* negatives. If
   AUROC jumps toward 0.66 there and sits near chance on repoDB, the published performance is
   substantially an artefact of the negative set — the most informative single result available.
5. **Network sweep**: all 7 interactomes plus random edge-removal curves (10–90%), to test the
   incompleteness/study-bias clause directly.
6. **Tissue arm**: GTEx-expressed-gene filtering per disease tissue.
7. **Statistics**: AUROC + AUPRC, bootstrap CIs resampled **by drug and by disease** (pairs are
   not independent), DeLong for paired AUC differences.
8. **Report honestly**: repoDB negatives are noisier than its positives (F1 0.55 vs 0.98) and
   ~30% of failures are administrative rather than efficacy-related. A low AUROC is partly
   attributable to these, and the write-up must say so rather than claiming the metric is
   worthless. Stratifying by `DetailedStatus` is the way to separate the two.
