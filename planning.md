# Planning: Do Topological Drug-Target Similarity Scores Predict Clinical Repurposing Success?

## Motivation & Novelty Assessment

### Why This Research Matters

Network proximity between a drug's protein targets and a disease's gene module
(Guney et al., *Nat. Commun.* 2016; 437 citations) has become the default
quantitative primitive of computational drug repurposing. It underpins COVID-19
repurposing pipelines (Morselli Gysi et al. 2021), drug-combination prediction
(Cheng et al. 2019), and dozens of derived tools. Trials cost $10-100M each, so
if the score is used to triage candidates, its true discriminative power for
*clinical* outcomes is a decision-relevant quantity. Yet essentially every
published evaluation scores it against a benchmark in which the negatives are
*unobserved* drug-disease pairs -- pairs nobody has ever tried. That is the wrong
question. The decision a repurposing programme actually faces is: given a
biologically plausible candidate, will it succeed in the clinic?

### Gap in Existing Work

From `literature_review.md`: the originating paper reports AUC = 0.66 under the
explicit assumption that all 18,162 unobserved drug-disease pairs are negatives.
repoDB (Brown & Patel 2017) exists precisely to fix that assumption -- it supplies
4,123 drug-indication pairs that entered clinical trials and *failed* -- but the
proximity literature has not been re-evaluated on it. Separately, the critique
literature (Lit-BM/HuRI study bias; "Genomic data integration systematically biases
interactome mapping", 2018; GhostBuster 2025; BETA 2022) argues that
literature-curated interactomes and text-mined disease genes carry a shared
study bias that inflates network-based metrics, but nobody has held the pair set
fixed and swept the interactome across the bias spectrum to measure it.

### Our Novel Contribution

Three things that, to our reading of the 74 papers gathered, have not been done
together:

1. **A hard-negative clinical benchmark.** Score proximity on repoDB, where the
   negatives are trial failures rather than untried pairs, against deliberately
   trivial baselines (target/disease-gene overlap count), with DeLong tests.
2. **A same-everything-but-the-negatives ablation.** Hold the positives, network,
   genes, measure and code path fixed and swap only the negative set
   (clinical failures vs Guney-style random unobserved pairs). This directly
   quantifies how much of the published AUC is a property of the benchmark rather
   than of the metric.
3. **Mechanism tests for the two named confounds.** A seven-interactome sweep
   spanning literature-curated to systematically-mapped (study-bias-free)
   networks, controlled random edge deletion, and GTEx tissue-restricted
   interactomes with a *shuffled-tissue control* that separates "tissue relevance"
   from "any sparsification".

### Experiment Justification

- **D1 (headline benchmark)** -- tests C1 and C2 directly. Without it there is no
  clinically grounded estimate of AUROC at all, and the hypothesis's two
  quantitative claims are untestable.
- **D4 (negative-set ablation, folded into D1)** -- the single most informative
  control available. If AUROC recovers toward 0.66 on random negatives while
  sitting near chance on clinical failures, the discrepancy is attributable to the
  benchmark, not to our implementation, our networks or our gene sets. Without
  this control a low AUROC is uninterpretable.
- **D5 (failure-reason stratification, folded into D1)** -- ~30% of repoDB
  failures are administrative (funding, accrual). If proximity were informative
  about *efficacy* it should do better against efficacy failures than against
  administrative ones. This is the honest way to avoid over-claiming from a noisy
  negative set.
- **D2 (network sweep + edge removal)** -- tests the "incompleteness" half of the
  hypothesised mechanism. A real biological signal should be robust across
  interactomes and degrade smoothly under edge deletion; a study-bias artefact
  should collapse on systematically-mapped networks.
- **D3 (tissue specificity)** -- tests the "tissue-specificity" half. The shuffled
  tissue control is what makes the result interpretable: without it, any change
  could be explained by sparsification alone.
- **Supervised comparator (GroupKFold by drug and by disease)** -- asks whether
  proximity adds anything at all on top of the trivial features, which is a
  stronger form of C2 than a single pairwise AUC comparison.


## Hypothesis under test

> Network proximity between drug targets and disease modules predicts successful drug
> repurposing with **AUROC below 0.65**, performing only **marginally better than
> drug-target overlap counts alone**, because PPI network incompleteness and
> tissue-specificity confound topological metrics.

Three separable claims:

| # | Claim | Falsifiable form |
|---|-------|------------------|
| C1 | Absolute performance is weak | AUROC(proximity) < 0.65 on a clinically grounded benchmark |
| C2 | No real gain over a trivial baseline | AUROC(proximity) − AUROC(overlap count) is small / not significant (DeLong) |
| C3 | Mechanism of the weakness | Performance varies systematically with network incompleteness / study bias and with tissue-specific filtering |

## Why the benchmark choice decides the answer

This is the single most important design decision, and it is what makes the
hypothesis non-trivial.

Guney et al. (2016), who introduced the proximity measure, report **AUC = 0.66**
for the best measure (`closest`, z-score) — already at the hypothesis boundary.
But that number comes from treating **all 18,162 unobserved drug–disease pairs as
negatives**; the authors state this assumption explicitly. Unobserved pairs are
mostly pharmacologically absurd (e.g. gliclazide for AML) and are trivially easy
to reject.

repoDB (Brown & Patel 2017) replaces that assumption: its negatives are drug–indication
pairs that **entered clinical trials and failed** (Terminated / Withdrawn / Suspended).
These are hard negatives — a human expert already judged them plausible enough to fund
a trial. Measuring on repoDB is therefore the honest operationalisation of
"predicts clinical repurposing **success**", and is expected to be much harder than
the original evaluation.

A 300-pair smoke test (Hetionet PPI × DisGeNET-curated, `closest`, n_random=100)
already gave AUROC ≈ 0.45 for proximity and ≈ 0.47 for overlap count — consistent
with C1 and C2, though far too small a sample to be a result. See
`notes/smoke_test.md`.

**Guney's own Fig. 2d also supports C2 directly**: drug–drug similarity by
*proximity* scored AUC 0.81 vs *shared targets* AUC 0.80, difference not
significant (P = 0.12).

## Direction budget: candidates considered

Ten plausible directions were enumerated and scored 1–5 on evidence support (E),
relevance to the hypothesis (R), expected information gain (I), and feasibility with
the data actually in hand (F).

| # | Direction | E | R | I | F | Total | Verdict |
|---|-----------|---|---|---|---|-------|---------|
| D1 | repoDB hard-negative benchmark: proximity (5 measures) vs overlap-count and degree baselines, with DeLong tests | 5 | 5 | 5 | 5 | **20** | **KEEP** |
| D2 | Network-sensitivity sweep: same pairs across 7 interactomes spanning literature-curated → systematic (bias/incompleteness axis) | 5 | 5 | 5 | 5 | **20** | **KEEP** |
| D3 | Tissue-specificity: GTEx/Hetionet-anatomy-filtered subnetworks; does tissue conditioning rescue or degrade proximity? | 4 | 5 | 4 | 4 | **17** | **KEEP** |
| D4 | Negative-set ablation: rerun with Guney-style random negatives to quantify how much of published AUC is benchmark artefact | 5 | 4 | 4 | 5 | 18 | Folded into D1 |
| D5 | Failure-reason stratification: efficacy failures vs administrative failures (funding/accrual) in repoDB `DetailedStatus` | 4 | 4 | 4 | 5 | 17 | Folded into D1 |
| D6 | Deep-learning repurposing models (deepDR, KG embeddings) as extra comparators | 3 | 2 | 3 | 2 | 10 | Pruned — hypothesis is about *topological proximity*, not SOTA-chasing; heavy compute; leakage-prone |
| D7 | Build a new tissue-specific interactome from single-cell atlases | 3 | 3 | 3 | 1 | 10 | Pruned — large data/compute cost, GTEx bulk suffices for D3 |
| D8 | Prospective validation against post-2017 ClinicalTrials.gov outcomes | 4 | 4 | 4 | 1 | 13 | Pruned — needs AACT ingest + long temporal follow-up; out of scope for this phase |
| D9 | Extend to drug *combinations* proximity | 3 | 1 | 2 | 3 | 9 | Pruned — different hypothesis |
| D10 | Chemical-structure / side-effect similarity baselines | 3 | 2 | 3 | 2 | 10 | Pruned — DrugBank structures gated; adds a non-topological axis the hypothesis does not claim |

### Retained directions (top 3)

**D1 — Hard-negative clinical benchmark (tests C1 + C2).**
Compute all five Guney distance measures (`closest`, `shortest`, `kernel`, `centre`,
plus raw un-normalised distance) and the degree-preserving z-scores over the 5,797
repoDB pairs. Compare against deliberately trivial baselines: drug-target/disease-gene
**overlap count**, Jaccard, number of targets, disease-module size, and drug/disease
degree. Report AUROC + AUPRC with bootstrap CIs and **DeLong** tests for paired
AUC differences. Includes D4 (swap in random negatives to isolate the benchmark
effect) and D5 (stratify failures by reason).

**D2 — Network sensitivity (tests C3, incompleteness/study-bias half).**
Hold pairs fixed; vary the interactome across the bias spectrum already assembled:
BioGRID full (976k edges, literature) → BioGRID low-throughput (139k, most
study-biased) → Lit-BM (12k, literature binary) → HuRI (52k, systematic Y2H,
bias-free by construction) → HI-union (63k) → STRING≥700 (233k) → Hetionet GiG
(147k). If proximity is a real biological signal it should be *robust*; if it
tracks literature bias it should collapse on HuRI. Also run edge-removal curves
(randomly delete 10–90% of edges) to measure sensitivity to incompleteness directly.

**D3 — Tissue specificity (tests C3, tissue half).**
Restrict the interactome to genes expressed in the tissue relevant to each disease
(GTEx v8 median TPM; Hetionet Anatomy–expresses–Gene as a second source) and
recompute proximity. Kitsak et al. (2016) argue disease modules are tissue-specific;
the hypothesis predicts tissue-naïve topology is confounded. Either outcome is
informative: improvement supports "tissue-specificity is a real confound",
no-improvement supports "topology is uninformative regardless".

### Pruning rule
Directions D6–D10 are not to be revisited unless D1–D3 produce evidence that
invalidates the ranking. Any such change must be recorded in STATE.md with its
justification.

## Statistical plan

- Primary metric **AUROC**, with **AUPRC** reported alongside (class balance ≈ 64/36).
- 2,000-resample bootstrap CIs, resampling **by drug** and **by disease** as well as
  by pair, since pairs are non-independent (one drug appears in many pairs).
- **DeLong** test for the paired proximity-vs-overlap comparison (this is C2's
  decision rule).
- Grouped cross-validation (`GroupKFold` on drug, and separately on disease) for any
  fitted comparator, to avoid the leakage that inflates published numbers.
- Pre-declared: C1 is supported if the upper bootstrap CI bound for AUROC is below
  0.65; C2 is supported if the DeLong P-value for the difference exceeds 0.05 or the
  point difference is below 0.02.

## Known confounds to control explicitly

1. **Degree / hub bias** — drug targets are better-studied and higher-degree
   (Guney: mean target degree 28.6 vs interactome 21.2). The z-score's
   degree-preserving randomisation is the intended correction; verify it works here
   by checking correlation of z with target degree (Guney reports ρ = −0.01 for z
   vs ρ = −0.46 for raw distance).
2. **Non-independence of pairs** — same drug across many indications.
3. **Disease-module size** — DisGeNET "all" modules are large (median 30 genes),
   which alone can drive distance.
4. **repoDB failure semantics** — ~30% of failures carry administrative reasons
   (funding, low accrual), not lack of efficacy. Stratify; do not silently treat
   them as pharmacological negatives.
5. **Cross-source ID loss** — CUI→DrugBank→Entrez joins lose ~46% of raw repoDB
   pairs; report coverage, and check that lost pairs are not systematically
   different in label.
