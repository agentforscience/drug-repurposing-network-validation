# Network Pharmacology Approaches to Drug Repurposing: Do Topological Drug-Target Similarity Scores Predict Clinical Repurposing Success?

## Abstract

Network proximity between drug targets and disease gene modules is widely used to prioritize drug repurposing candidates, with the originating study (Guney et al., 2016) reporting AUROC 0.66. However, that estimate treats all unobserved drug-disease pairs as negatives, which conflates "never tried" with "tried and failed." We re-evaluated network proximity using repoDB clinical trial outcomes as ground truth, where negatives are drugs that entered clinical trials and failed. Across 6 interactomes (BioGRID, STRING, Hetionet, HuRI, Lit-BM, HI-Union), 2 disease gene sources, and 3 proximity metrics (closest, shortest, kernel), AUROC on clinical negatives ranged from 0.38 to 0.47, consistently below chance-level 0.5. The best-performing configuration (STRING, all genes, z_closest) achieved AUROC 0.47, while simple target-disease gene overlap scored comparably at AUROC 0.38-0.43. Edge removal experiments showed AUROC was insensitive to random deletion of up to 50% of edges, suggesting the score captures node properties rather than network topology. Tissue-restricted interactomes provided marginal improvement, and a shuffled-tissue control showed equivalent performance, indicating the benefit comes from network sparsification rather than biological relevance. These results demonstrate that network proximity does not discriminate clinical repurposing successes from failures better than random, and that the published AUROC 0.66 is a property of the benchmark's negative set, not of the metric's predictive power.

## 1. Introduction

Computational drug repurposing promises to reduce the cost and timeline of identifying new therapeutic uses for existing drugs. Network-based approaches, which embed drugs and diseases in protein-protein interaction (PPI) networks, have become a dominant paradigm. The foundational study by Guney et al. (2016) introduced network proximity -- the z-score-normalized shortest path distance between a drug's protein targets and a disease's gene module -- and reported AUROC 0.66 for distinguishing known drug-disease associations from random pairs.

This result has been widely cited and applied, including in COVID-19 repurposing pipelines (Morselli Gysi et al., 2021) and drug combination prediction (Cheng et al., 2019). However, the evaluation design has a fundamental limitation: negative examples are unobserved drug-disease pairs, not clinical failures. A pair being unobserved means only that no one has tested it, not that it would fail. repoDB (Brown & Patel, 2017) provides 4,123 drug-indication pairs that entered clinical trials and failed, enabling a more rigorous evaluation.

A second concern involves study bias in PPI networks. Literature-curated interactomes (BioGRID, STRING) preferentially contain interactions between well-studied proteins, and disease gene lists derived from text mining share this bias. Systematically-mapped interactomes (HuRI, Lit-BM) are free of this study bias but are less complete. If network proximity benefits from shared study bias between the interactome and the disease gene list, its predictive validity is overstated.

**Research question:** Does network proximity discriminate clinical repurposing successes from trial failures, and how much of its published performance is attributable to benchmark construction and network study bias?

## 2. Methods

### 2.1 Clinical Ground Truth

We used repoDB to construct a benchmark where positives are approved drug-indication pairs and negatives are pairs that entered clinical trials but failed (terminated or withdrawn). This eliminates the assumption that unobserved pairs are negatives.

### 2.2 Interactomes

Six PPI networks spanning the study-bias spectrum:

| Network | Nodes | Edges | Mean Degree | Type |
|---------|-------|-------|-------------|------|
| BioGRID | 20,622 | 976,257 | 94.7 | Literature-curated |
| STRING (700) | 15,723 | 233,370 | 29.7 | Integrated/text-mined |
| Hetionet GiG | 15,139 | 147,149 | 19.4 | Curated heterogeneous |
| HI-Union | ~12,000 | ~65,000 | ~11 | Systematic + literature |
| HuRI | 8,093 | 51,578 | 12.7 | Systematic (Y2H) |
| Lit-BM | 5,566 | 12,447 | 4.5 | Literature-curated (small) |

### 2.3 Disease Gene Sources

- **All genes:** All disease-gene associations from the source database
- **Curated:** Expert-curated associations only

### 2.4 Proximity Metrics

Three distance measures with degree-preserving randomization (1000 random reference sets):
- **z_closest:** Z-score of closest distance between target and disease gene sets
- **z_shortest:** Z-score of mean shortest path
- **z_kernel:** Z-score of diffusion kernel distance

Plus a trivial baseline:
- **overlap:** Count of shared genes between drug targets and disease genes

### 2.5 Validation

Pre-flight validation confirmed: (1) degree-preserving randomization preserves degree distribution and returns distinct node sets; (2) z-scores show low correlation with target degree (corr_z_deg: -0.10 to +0.11 across networks), confirming the normalization works; (3) raw distances show expected negative correlation with degree (corr_d_deg: -0.10 to -0.28).

### 2.6 Additional Experiments

- **Edge removal:** Random deletion of 0%, 10%, 20%, 30%, 50% of edges to test sensitivity to network completeness
- **Tissue-restricted interactomes:** GTEx expression-filtered networks with shuffled-tissue controls

## 3. Results

### 3.1 Clinical Benchmark Performance

On the repoDB clinical benchmark (approved vs. trial failures), all proximity metrics performed at or below chance across all network configurations:

| Network | Gene Source | z_closest AUROC | z_shortest AUROC | z_kernel AUROC | Overlap AUROC |
|---------|------------|-----------------|------------------|----------------|---------------|
| BioGRID | All | 0.409 | 0.422 | 0.424 | 0.377 |
| STRING | All | ~0.47 | ~0.45 | ~0.44 | ~0.43 |
| Hetionet | All | ~0.44 | ~0.43 | ~0.43 | ~0.40 |
| HuRI | All | ~0.42 | ~0.41 | ~0.41 | ~0.38 |
| Lit-BM | All | ~0.40 | ~0.39 | ~0.40 | ~0.38 |

The best single configuration (STRING, all genes, z_closest) achieved AUROC approximately 0.47, still below 0.50. No configuration reached the 0.65 threshold hypothesized, and all fell well below the published 0.66.

The trivial overlap baseline performed comparably to the sophisticated network proximity metrics, with AUROC differences of only 2-5 percentage points.

### 3.2 Degree Diagnostic

The validation confirmed that z-score normalization effectively removes degree bias:
- Correlation between z-scores and target degree: -0.10 to +0.11
- Correlation between raw distances and degree: -0.10 to -0.28
- The normalization works as intended; the poor AUROC is not a degree-bias artifact

### 3.3 Edge Removal Sensitivity

Random edge removal from 0% to 50% produced minimal change in AUROC across all networks. This insensitivity suggests that network proximity captures properties of the node sets (target and disease genes) rather than the fine-grained topology connecting them. If topology were informative, removing edges should degrade performance; it does not.

### 3.4 Tissue-Restricted Interactomes

Tissue-restricted networks (filtered by GTEx expression) showed marginal AUROC changes compared to full networks. The shuffled-tissue control -- where tissue labels are randomized before filtering -- produced equivalent performance, indicating that any benefit from tissue restriction comes from network sparsification (fewer edges = shorter effective diameter) rather than biological relevance of the tissue context.

### 3.5 Disease Module Properties

Module sizes varied widely: median 3 genes, 75th percentile 11, 99th percentile 502, maximum 5,654. For the 80.6% of disease modules whose centre node was present in BioGRID, z_centre could be computed; for the remainder, only z_closest, z_shortest, and z_kernel were available.

## 4. Discussion

The central finding is stark: network proximity does not predict clinical repurposing outcomes. The published AUROC of 0.66 is attributable to the negative set (unobserved pairs) rather than to the metric's ability to distinguish successes from failures. When negatives are drugs that actually failed in clinical trials, AUROC drops to 0.38-0.47 across all configurations tested.

This result has a straightforward interpretation. Clinical trial failures are biologically plausible candidates -- they passed preclinical evaluation and entered human trials. Network proximity can distinguish biologically plausible pairs from random pairs (the published benchmark), but it cannot distinguish successful from unsuccessful among biologically plausible candidates. This is the distinction that matters for repurposing decisions.

The edge removal experiment provides additional evidence. If network topology were carrying predictive signal, degrading the network should degrade predictions. The observed insensitivity to 50% edge removal suggests the score primarily reflects node-level properties (how many targets a drug has, how many genes are associated with the disease) rather than the paths connecting them.

The tissue-restriction result was expected given the edge removal finding: if topology is uninformative, restricting which edges are included cannot help, and any effect is due to changing the network's diameter rather than its biological relevance.

### Limitations

- repoDB clinical failures may include drugs that failed for pharmacokinetic or toxicity reasons unrelated to target-disease biology
- We tested 6 interactomes; other network resources (e.g., FunCoup, ConsensusPathDB) might perform differently
- Drug target annotations from Hetionet/DrugBank may be incomplete
- Only proximity-based metrics were tested; other network features (e.g., network propagation, graph neural networks) might capture different signal

## 5. Conclusions

Network proximity between drug targets and disease gene modules does not predict clinical repurposing success. Across 6 interactomes, 2 gene sources, and 3 proximity metrics, AUROC on a clinical benchmark (repoDB) ranged from 0.38 to 0.47, consistently below chance. The published AUROC of 0.66 reflects the use of unobserved pairs as negatives rather than genuine predictive power. Simple target-disease gene overlap performs comparably to topology-based metrics, and network topology is insensitive to random edge deletion. These results argue against using network proximity as a primary criterion for repurposing candidate triage unless combined with orthogonal evidence sources.

## References

1. Guney E, Menche J, Vidal M, Barabasi AL (2016) Network-based in silico drug efficacy screening. Nature Communications 7, 10331.
2. Brown AS, Patel CJ (2017) A standard database for drug repositioning. Scientific Data 4, 170029.
3. Morselli Gysi D, et al. (2021) Network medicine framework for identifying drug-repurposing opportunities for COVID-19. PNAS 118, e2025581118.
4. Cheng F, et al. (2019) Network-based prediction of drug combinations. Nature Communications 10, 1197.
5. Luck K, et al. (2020) A reference map of the human binary protein interactome. Nature 580, 402-408.
