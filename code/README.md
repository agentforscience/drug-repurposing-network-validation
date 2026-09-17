# Code

## Written for this project

| File | Purpose |
|---|---|
| `proximity.py` | **Python 3 port of the Guney et al. (2016) proximity measure.** The core method under test. |
| `build_dataset.py` | Assembles everything into `datasets/processed/` (networks, drug targets, disease genes, labelled pairs, coverage report). |
| `download_datasets.sh` | Raw dataset downloads. Idempotent — skips non-empty files. |
| `search_lit.py`, `search_lit2.py` | Literature search over Europe PMC / Semantic Scholar / arXiv. |
| `rank_papers.py` | Scores and ranks candidates by relevance to the hypothesis. |
| `download_papers.py`, `retry_papers.py` | Open-access PDF retrieval (Europe PMC → S2 → OpenAlex → arXiv/bioRxiv). |

### `proximity.py` — usage

```python
import sys; sys.path.insert(0, 'code')
import pandas as pd, networkx as nx
from proximity import ProximityCalculator, overlap_count

e = pd.read_csv('datasets/processed/ppi_hetionet_gig.tsv', sep='\t')
G = nx.from_pandas_edgelist(e, 'gene_a', 'gene_b')

pc = ProximityCalculator(G, bin_size=100, n_random=1000, seed=452456)
res = pc.proximity(drug_targets, disease_genes, distance='closest')
# -> {'d':…, 'z':…, 'mu':…, 'sigma':…, 'n_targets':…, 'n_genes':…}
```

Implements equations 1–4 of the paper (`closest`, `shortest`, `kernel`, `centre`) and the
degree-preserving reference distribution: degree bins grown until each holds ≥100 nodes,
1,000 random target/disease set pairs matched bin-wise, `z = (d − μ)/σ`.

**Performance.** The constructor precomputes the all-pairs shortest-path matrix once
(unweighted BFS via `scipy.sparse.csgraph`, stored `uint8`), after which every evaluation is
numpy indexing. Measured on Hetionet GiG (15,139 nodes / 147,149 edges):

- APSP + binning: **41 s**, 229 MB resident
- Proximity: **~3 ms/pair at n_random=100** → the full 5,797-pair sweep at n_random=1000 is
  minutes per network, not hours.

`overlap_count()` and `jaccard()` provide the trivial baselines the hypothesis compares against.

**Verified**: runs end-to-end on real repoDB pairs (`notes/smoke_test.md`); network is fully
connected (max hop distance 9, no unreachable pairs), 47 degree bins.

## Cloned repositories

### `emreg00_toolbox/` — github.com/emreg00/toolbox
The **authoritative reference implementation** of network proximity, by Guney (first author of
the 2016 paper). Key entry points: `wrappers.py::calculate_proximity`,
`wrappers.py::calculate_closest_distance`, `network_utilities.py::get_degree_binning`,
`network_utilities.py::pick_random_nodes_matching_selected`.

⚠️ **Python 2 only** — uses `print` statements, `iteritems()`, `xrange`. It will not run under
this project's Python 3.12 environment. It was cloned as the specification against which
`code/proximity.py` was written, and `proximity.py` is what should actually be used.
Also contains ~60 `parse_*.py` modules for the source databases (DrugBank, DisGeNET, Hetionet,
STRING, DailyMed, clinical trials …) which are useful documentation of the original data joins.

### `repoDB/` — github.com/adam-sam-brown/repoDB
Source repository for the repoDB gold standard (Brown & Patel 2017). `Reproduction_Scripts/`
holds the R code that built it from DrugCentral + AACT; `Shiny_Application/data/shiny.RData`
holds the full database, which is where `datasets/repodb/repodb_full.csv` came from.

## Environment

Isolated `uv` venv at `.venv/` (Python 3.12.8), dependencies pinned in `pyproject.toml`:
`pandas, numpy, scipy, networkx, scikit-learn, statsmodels, matplotlib, pypdf, requests,
httpx, pyreadr`.

```bash
source .venv/bin/activate
uv add <package>     # adds to pyproject.toml
```

Note `[tool.uv] package = false` in `pyproject.toml` — without it `uv` tries to build the
workspace as a wheel and fails.
