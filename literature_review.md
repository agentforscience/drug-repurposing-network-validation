# Literature Review

**Topic**: Network pharmacology approaches to drug repurposing — do topological
drug–target similarity scores predict clinical repurposing success?

**Hypothesis**: Network proximity between drug targets and disease modules predicts
successful drug repurposing with AUROC below 0.65, performing only marginally better
than drug-target overlap counts alone, because PPI network incompleteness and
tissue-specificity confound topological metrics.

74 open-access PDFs downloaded (of 110 top-ranked candidates; 1,192 unique papers
screened by abstract). Full catalogue in `papers/README.md`.

---

## 1. Research area overview

Network medicine holds that disease-associated genes cluster into connected
"disease modules" in the human interactome, and that a drug works if its targets sit
inside or adjacent to that module. This yields a concrete repurposing score:
the **network proximity** between a drug's target set and a disease's gene set.
The field is large and largely self-confirming — most papers propose a proximity
variant, apply it to one disease (very often COVID-19), and validate on a handful of
literature-supported examples.

The hypothesis under test targets the field's weak point: **how proximity is
evaluated**, not how it is computed. Two facts from the primary literature make it
sharply testable.

1. **The canonical result is already borderline.** Guney et al. (2016), who
   introduced the measure, report **AUC = 0.66** for their best variant.
2. **The canonical evaluation assumes its negatives.** That 0.66 comes from
   treating all 18,162 unobserved drug–disease pairs as negatives — an assumption the
   authors state explicitly and that Brown & Patel (2017) built repoDB to remove.

---

## 2. Key papers

### 2.1 Guney, Menche, Vidal & Barabási (2016) — *Network-based in silico drug efficacy screening*
`papers/012_Network_based_in_silico_drug_efficacy_screening.pdf` · Nat Commun 7:10331 · 437 cites · **read in full**

The foundational paper; defines the metric this project evaluates.

- **Method**: five distance measures between drug target set *T* and disease gene set *S* —
  *closest* (mean over targets of the minimum distance to any disease gene), *shortest*,
  *kernel* (exponentially down-weighted paths), *centre*, *separation*. Each is converted to a
  z-score against a reference distribution from **degree-preserving random node sets**
  (1,000 permutations; degree bins built so each bin holds ≥100 nodes).
- **Data**: interactome of 141,150 interactions / 13,329 proteins (Menche et al. 2015);
  disease genes from OMIM + GWAS catalog, diseases restricted to ≥20 genes;
  drug targets from DrugBank ("Targets" only, excluding enzymes/carriers/transporters);
  indications from MEDI-HPS. Final set: **238 drugs × 78 diseases = 402 known pairs**.
- **Headline results**:
  - AUC by measure: **closest 0.66**, kernel 0.61, separation 0.59, shortest 0.58, centre 0.58.
  - Raw distance *d*<sub>c</sub> correlates with target degree (ρ = −0.46); the z-score does
    not (ρ = −0.01). The degree-preserving randomisation is therefore doing real work.
  - Drug targets are higher degree than the network average (28.6 vs 21.2) — the authors
    attribute this to **literature bias toward drug targets**.
  - Only 15.4% of known pairs have any drug-target/disease-gene overlap.
  - 165 of 402 known pairs are *distant* — the interactome cannot explain them.
- **Directly relevant to the hypothesis**: in their own Fig. 2d, drug–drug similarity by
  **proximity scores AUC 0.81 vs shared-target count AUC 0.80, P = 0.12 (not significant)**.
  This is the "only marginally better than drug-target overlap" claim, in the originating paper.
- **Also relevant**: Table 1 shows drugs that failed for *lack of efficacy* have z ≈ 0.0–1.8
  (distant, as the theory predicts), but drugs withdrawn for *adverse effects* are strongly
  proximal (semagacestat z = −5.6, terfenadine z = −2.2). Proximity appears to measure
  pharmacological engagement, which is not the same thing as clinical success.
- **Code**: `github.com/emreg00/toolbox` (cloned to `code/emreg00_toolbox/`; Python 2).

### 2.2 Brown & Patel (2017) — *A standard database for drug repositioning* (repoDB)
`papers/010_A_standard_database_for_drug_repositioning.pdf` + `papers/fulltext_xml/PMC5349249.xml` · Sci Data 4:170029 · 206 cites · **read in full**

The paper that supplies this project's ground truth, and its rationale is the hypothesis.

- **Motivation**: repositioning papers "typically assume that all other drug-indication pairs
  are false… This assumption is unsatisfying [and] suggests that all novel repositioning
  predictions are false."
- **Positives**: 6,677 approved drug–indication pairs from DrugCentral (UMLS-mapped from FDA
  labels via the OMOP pipeline, **F1 ≈ 0.98**).
- **Negatives**: 4,123 pairs from AACT/ClinicalTrials.gov with phase 0–3 and status
  *suspended / terminated / withdrawn*, mapped via NLM Medical Text Indexer (**F1 ≈ 0.55**).
- Pairs later approved were removed from the failed set, so labels do not conflict.
- Scale: 1,571 drugs × 2,051 UMLS indications.
- **Caveat this project must carry**: the positive and negative labels are not annotated to
  the same quality (F1 0.98 vs 0.55). Noisier negatives depress attainable AUROC
  independently of whether proximity works. Any low AUROC must be reported with this caveat,
  not presented purely as evidence against proximity.

### 2.3 Menche et al. (2015) — *Uncovering disease-disease relationships through the incomplete interactome*
Science 347:1257601 · 1,089 cites · paywalled (abstract retained)

Origin of the disease-module formalism and the *separation* measure S<sub>AB</sub>; supplies the
interactome Guney used. Derives **mathematical conditions for when a disease module is
identifiable at all** given interactome incompleteness — the theoretical backing for the
hypothesis's "incompleteness confounds topological metrics" clause. Modules are only
detectable when enough disease genes are known and enough of the interactome is mapped.

### 2.4 Cheng et al. (2018) — *Network-based approach to prediction and population-based validation of in silico drug repurposing*
`papers/006_Network_based_approach_to_prediction_and_population_based_validati.pdf` · Nat Commun · 392 cites

Applies proximity at scale and adds validation against **patient-level insurance claims** rather
than literature annotation — the strongest external-validity design in this literature, and the
closest existing analogue to asking whether proximity tracks clinical reality.

### 2.5 Morselli Gysi et al. (2021) — *Network medicine framework for identifying drug-repurposing opportunities for COVID-19*
`papers/011_Network_medicine_framework_for_identifying_drug_repurposing_opport.pdf` · PNAS · 345 cites

Important as a **cautionary result**: the authors run multiple network/AI pipelines and find
the methods disagree substantially with one another, resolving this by ensembling and
validating against outcomes. Evidence that individual proximity rankings are not stable.

### 2.6 Kitsak et al. (2016) — *Tissue specificity of human disease module*
`papers/019_Tissue_Specificity_of_Human_Disease_Module.pdf` · Sci Rep · 77 cites

The basis for the tissue-specificity arm. Genes expressed in the same tissue localise in the
same interactome neighbourhood; disease manifestation depends on the **integrity and
completeness of the disease module's expression in that tissue**. Implies a tissue-naive
interactome mixes together modules that are never co-expressed — precisely the confound named
in the hypothesis.

### 2.7 Luck et al. (2020) — *A reference map of the human binary protein interactome* (HuRI)
`papers/052_A_reference_map_of_the_human_binary_protein_interactome.pdf` · Nature · 1,083 cites

HuRI is ~53,000 systematically screened binary interactions — testing all pairs rather than
following literature attention. This makes it **free of study bias by construction**, which is
what turns "PPI incompleteness/bias confounds proximity" from an assertion into an experiment:
run identical pairs on literature-curated vs systematic networks and compare.

### 2.8 Study-bias and evaluation-methodology papers
- **Stacey et al. (2018), *Genomic data integration systematically biases interactome mapping***
  (`papers/066_...pdf`) — integrating genomic data raises apparent functional coherence but
  *reduces* power to find novel interactions, and its novel predictions are **no more likely to
  be experimentally confirmed** than predictions without it. A concrete demonstration that
  network "quality" improvements can be circular.
- **GhostBuster (2025)** (`papers/041_...pdf`) — literature bias means well-studied genes
  dominate; models trained on literature-derived labels score higher while being *less* able
  to recover newly discovered annotations. Direct evidence that standard metrics
  "overestimate biological relevance by overfitting literature-derived patterns."
- **Benchmarking data leakage in biomedical KGE (2026)** (`papers/058_...pdf`) — train/test
  redundancy inflates biomedical link prediction; performance drops substantially on an
  independent Orphanet test set. Notably, they find **no evidence that node degree alone drives
  predictions**, a useful counterpoint to assume-the-worst degree-bias arguments.
- **BETA (2022)** (`papers/070_...pdf`) — random-shuffle cross-validation is inadequate for
  drug–target prediction; proposes connectivity- and category-stratified splits. Motivates
  grouped CV by drug and by disease here.
- **IxIDN/ORDON (2026)** (`papers/021_...pdf`) — very recent, independent effort building
  clinical-trial-derived positives and ontology-derived negatives for rare disease
  repurposing. Confirms the negative-set problem is still live nine years after repoDB.

---

## 3. Common methodology in this literature

| Component | Usual choice | Notes for us |
|---|---|---|
| Interactome | Literature-curated union (Menche/Cheng ~140–350k edges) | Carries study bias; almost never varied as a sensitivity analysis |
| Disease genes | OMIM + GWAS, or DisGeNET; often ≥20-gene cut | Module size varies hugely and is itself predictive |
| Drug targets | DrugBank pharmacological targets only | Excluding enzymes/carriers/transporters matters |
| Score | z-score of `closest` distance, 1,000 degree-preserving permutations | Guney's `closest` is the field default |
| Negatives | **Unobserved pairs assumed false** | The central methodological weakness |
| Metric | AUROC, occasionally AUPRC | Rarely with CIs; DeLong tests rare |

## 4. Standard baselines
- Drug-target/disease-gene **overlap count** — the trivial baseline the hypothesis names.
- **Shared-target** drug–drug similarity (Guney: AUC 0.80).
- Chemical structure similarity (0.78), GO-term similarity (0.71), shortest-path target
  similarity (0.71), LINCS expression similarity (0.65), side-effect similarity (0.81, low coverage).
- Gene-expression-signature correlation (DvD): AUC 0.53 — barely above chance.
- **Not usually reported, and needed here**: number of drug targets, disease module size, and
  target degree, each of which can produce apparent signal on its own.

## 5. Evaluation metrics
AUROC is the field standard and the hypothesis is stated in AUROC, so it is primary.
Add AUPRC (repoDB is ~64% positive), bootstrap CIs resampled by drug and by disease (pairs are
not independent), and **DeLong** paired tests for proximity-vs-overlap — that comparison is the
hypothesis's second claim and needs a significance test, not a point estimate.

## 6. Datasets used in the literature
- Interactomes: Menche/Cheng curated union, STRING, BioGRID, HuRI/HI-union, Lit-BM.
- Disease genes: OMIM, GWAS Catalog, DisGeNET, DISEASES (Jensen lab).
- Drug targets: DrugBank, DrugCentral, ChEMBL, Hetionet.
- Gold standards: repoDB, PREDICT/Gottlieb, Fdataset/Cdataset, Hetionet CtD, ClinicalTrials.gov.

## 7. Gaps and opportunities

1. **The negative-set gap.** Proximity's headline AUC is measured against negatives assumed
   to be false. repoDB's clinically failed pairs are far harder — and are the correct
   operationalisation of "clinical repurposing success". Almost no proximity paper uses them.
2. **The network-sensitivity gap.** Proximity is nearly always reported on one interactome.
   Guney does compare against STRING and binary screens in a supplementary table (finding
   their curated network best), but the systematic sweep across the bias spectrum is missing.
3. **The trivial-baseline gap.** Overlap count is rarely reported as a standalone comparator,
   despite Guney's own data showing shared-target similarity is statistically
   indistinguishable from proximity-based similarity.
4. **Failure-reason conflation.** repoDB negatives mix efficacy failures with administrative
   ones (funding, low accrual). No proximity study stratifies these.
5. **Tissue-naivety.** Kitsak shows modules are tissue-specific; proximity is computed on a
   tissue-agnostic network anyway.

## 8. Recommendations for our experiment

- **Ground truth**: repoDB — 5,797 pairs (3,708 approved / 2,089 failed) survive the join to
  drug targets and disease genes. Report the coverage loss and check it is not label-dependent.
- **Primary score**: Guney `closest` z-score, 1,000 degree-preserving permutations, faithful to
  the original. Also compute `shortest`, `kernel`, `centre`, and raw un-normalised distance.
- **Baselines**: overlap count (primary comparator), Jaccard, number of targets, module size,
  mean target degree.
- **Networks**: sweep BioGRID (976k edges) → BioGRID low-throughput (139k) → STRING≥700 (233k)
  → Hetionet GiG (147k) → HI-union (63k) → HuRI (52k) → Lit-BM (12k). Add random edge-removal
  curves to probe incompleteness directly.
- **Tissue arm**: restrict to GTEx-expressed genes per disease tissue and recompute.
- **Ablation**: rerun with Guney-style random negatives to quantify how much of the published
  AUC is an artefact of the negative set. This is the single most informative control.
- **Statistics**: AUROC + AUPRC, bootstrap CIs by drug and by disease, DeLong for the
  proximity-vs-overlap difference.

### Methodological cautions
- Do not read a low AUROC as pure evidence against proximity: repoDB negatives are noisier
  (F1 ≈ 0.55) than positives (F1 ≈ 0.98), and this alone depresses attainable performance.
- Verify the degree-correction works in *our* setup (Guney: ρ(z, degree) ≈ −0.01 vs
  ρ(d, degree) ≈ −0.46). If it does not, degree bias is uncontrolled.
- Pairs are not independent — one drug spans many indications. Grouped resampling is required.
- Proximity may track pharmacological engagement rather than therapeutic benefit: Guney's
  adverse-effect withdrawals were *strongly proximal*. Consider this when interpreting.
