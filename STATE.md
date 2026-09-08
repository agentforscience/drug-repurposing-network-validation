# Research State

- Current phase: `experiment_runner`
- Pipeline completed: `False`

## Previous phases

resource_finder (succeeded)

## Current phase context

- Phase: `experiment_runner`
- Status: `in_progress`
- Started: `2026-09-08T18:47:04.562046Z`
- Next steps:
  - Validate the report and experimental artifacts before finalizing.

## Workspace check

- Expected: `/workspaces/network_pharmacology_approache_20260908_172336_17b06cb2`
- Actual: `/app`
- Directory usable: `True`
- Current process matches workspace: `False`

## Output validation

No phase output validation recorded yet.

## Agent notes

<!-- NEURICO_AGENT_NOTES_START -->
### resource_finder
<!-- NEURICO_AGENT_NOTES_START:resource_finder -->
### resource_finder

**Status**: COMPLETE. All expected artifacts exist on disk.

**Environment**: isolated `uv` venv at `.venv/` (Python 3.12.8), deps in `pyproject.toml`.
Requires `[tool.uv] package = false` or `uv add` fails. Activate with
`source .venv/bin/activate`. The experiment runner should reuse this venv.

**What was completed**
- Literature: paper-finder service was DOWN; fell back to Europe PMC + Semantic Scholar +
  arXiv across 30 queries in 2 rounds (round 2 deliberately targeted the critique /
  evaluation-pitfalls literature). 1,477 hits -> 1,192 unique -> top 110 pursued ->
  **74 PDFs in `papers/`** + 5 full-text XML. 36 paywalled (abstracts retained in
  `paper_search_results/ranked.jsonl`).
- Deep-read in full: Guney et al. 2016 (the method under test) and Brown & Patel 2017 (repoDB).
- Datasets: 12 sources downloaded, then assembled by `code/build_dataset.py` into
  `datasets/processed/` -> **7 interactomes x 2 disease-gene definitions x up to 5,797
  labelled repoDB pairs**.
- Implemented + validated `code/proximity.py` (Python 3 port of Guney proximity).

**Key findings driving the experiment design**
1. Guney et al. 2016 report **AUC = 0.66** for their best measure (`closest` z-score) -
   already at the hypothesis's 0.65 boundary. Other measures: kernel 0.61, separation 0.59,
   shortest 0.58, centre 0.58.
2. That 0.66 assumes **all unobserved drug-disease pairs are negatives** (stated explicitly
   in their Methods). repoDB replaces these with pairs that entered clinical trials and
   FAILED - genuinely hard negatives, and the correct reading of "clinical repurposing success".
3. Guney's own Fig. 2d: proximity-based drug similarity AUC 0.81 vs **shared-target count
   AUC 0.80, P = 0.12 (not significant)** - i.e. the hypothesis's "marginally better than
   drug-target overlap" claim appears in the originating paper.
4. Guney Table 1: drugs withdrawn for ADVERSE EFFECTS are strongly proximal (z = -5.6, -2.2),
   while efficacy failures are distant. Proximity may track pharmacological engagement, not
   therapeutic benefit. Important for interpretation.
5. repoDB label-quality asymmetry: positives F1 ~0.98 (DrugCentral/OMOP), negatives F1 ~0.55
   (NLM Medical Text Indexer). ~30% of failures are administrative (funding, low accrual),
   not efficacy. **Stratify by `DetailedStatus`; do not read a low AUROC as pure evidence
   against proximity.**

**Direction budget** (10 enumerated, scored, top 3 kept - see `planning.md` for the full
table and pruning reasons): D1 repoDB hard-negative benchmark w/ overlap baselines + DeLong
(absorbs D4 random-negative ablation and D5 failure-reason stratification); D2 network
sensitivity sweep across the study-bias spectrum + edge-removal curves; D3 tissue-specific
filtering via GTEx. Pruned: D6 deep-learning comparators, D7 single-cell interactome
construction, D8 prospective ClinicalTrials.gov validation, D9 drug combinations,
D10 chemical/side-effect similarity.

**Evidence paths**
- `planning.md` - hypothesis decomposition (C1/C2/C3), direction budget, statistical plan
- `literature_review.md` - synthesis, methods table, gaps, recommendations
- `resources.md` - full catalogue + gathering notes
- `notes/smoke_test.md` - end-to-end validation run
- `datasets/processed/coverage.json` - pairs surviving every network x gene-source join
- `code/proximity.py` - the method under test; `code/emreg00_toolbox/` - Python 2 reference

**Validated smoke test (NOT a result)**: 300 sampled pairs, Hetionet GiG x DisGeNET-curated,
n_random=100 -> AUROC 0.454 (proximity z), 0.465 (overlap count), 0.407 (raw distance),
0.401 (n targets), 0.360 (module size). Consistent with the hypothesis but far too small to
conclude anything. Runtime **~3 ms/pair**; full 5,797-pair x n_random=1000 sweep is minutes
per network after a one-off 41 s APSP build (229 MB uint8 matrix).

**Gotchas the next phase will hit**
- `code/emreg00_toolbox/` is **Python 2** and will not run. Use `code/proximity.py`.
- Europe PMC `fullTextPdf` REST route 404s; use `https://europepmc.org/articles/{PMCID}?pdf=render`.
- HuRI must be fetched from `https://interactome-atlas.org/` (the `www.` host 301s, body dropped).
- DisGeNET live download is gated; use the Zenodo `dhimmel/disgenet` copy already in `datasets/`.
- Semantic Scholar rate-limits aggressively (429) - budget backoff if searching more.
- Pairs are non-independent (one drug spans many indications): resample by drug AND by disease.

**Next phase**: `experiment_runner`. Concrete next steps:
1. Run D1: all 5 distance measures + 5 baselines over `biogrid x disgenet_all` (5,797 pairs),
   n_random=1000. AUROC/AUPRC + bootstrap CIs + DeLong(proximity vs overlap).
2. Run the random-negative ablation (D4) - the single most informative control.
3. Run D2 network sweep across all 7 interactomes + edge-removal curves.
4. Run D3 tissue filtering with GTEx.
5. Sanity check that degree correction works here: corr(z, target degree) should be ~0 while
   corr(d, degree) should be strongly negative (Guney: -0.01 vs -0.46). If not, degree bias
   is uncontrolled and results are not interpretable.

**Unresolved uncertainty**
- Disease-gene coverage caps the benchmark at 5,797/10,800 raw pairs. The retained set is 64%
  positive vs 62% raw - close, but **verify per-configuration that dropped pairs are not
  label-biased** before reporting.
- DisGeNET "all" includes text-mined associations, which carry the same literature bias the
  hypothesis blames. Running both `curated` and `all` is the intended control, but neither is
  bias-free; this limits how cleanly C3 can be attributed to the *network* alone.
<!-- NEURICO_AGENT_NOTES_END:resource_finder -->

### experiment_runner
<!-- NEURICO_AGENT_NOTES_START:experiment_runner -->
Update this section at the end of the `experiment_runner` phase.
<!-- NEURICO_AGENT_NOTES_END:experiment_runner -->

<!-- NEURICO_AGENT_NOTES_END -->
