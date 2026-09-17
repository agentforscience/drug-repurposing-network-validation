# Pipeline smoke test (NOT a result)

Purpose: verify the assembled data + Python 3 proximity port run end to end and
produce sane numbers. Sample is far too small and n_random too low to conclude
anything.

Configuration: Hetionet GiG PPI (15,139 nodes / 147,149 edges), DisGeNET-curated
disease genes, Hetionet CbG drug targets, `closest` measure, n_random = **100**,
**300 randomly sampled** repoDB pairs (of 2,972 eligible).

| Score | AUROC |
|---|---|
| Proximity z (−z) | 0.454 |
| Raw closest distance (−d) | 0.407 |
| Drug-target/disease-gene overlap count | 0.465 |
| Number of drug targets | 0.401 |
| Disease module size | 0.360 |

Mean z: approved 0.188, failed 0.070.

Observations to carry into the real run:
- Runtime is **~3 ms/pair at n_random=100** after the one-off APSP build (41 s,
  229 MB uint8 matrix). The full 5,797-pair × n_random=1000 run is therefore
  minutes, not hours — all 7 networks are affordable.
- Numbers sit near or below chance, consistent with the hypothesis, but a 300-pair
  sample has a bootstrap CI roughly ±0.07, so this distinguishes nothing yet.
- Both "number of targets" and "module size" carry signal in the same direction as
  proximity. These must be reported as baselines, not left as silent confounds.

---

## Network-sensitivity spot check (worth following up in D2)

Metformin (DB00331) / type 2 diabetes (C0011860) — a canonical *approved* pair, and one of
the examples Guney et al. discuss — scored with the identical `closest` measure on two
different interactomes:

| Network | nodes | \|T\| | \|S\| | d | z |
|---|---|---|---|---|---|
| Lit-BM | 5,566 | 12 | 71 | 2.25 | **−2.03** (proximal) |
| BioGRID | 20,622 | 56 | 122 | 1.52 | **+2.38** (distant) |

The raw distance *falls* on the denser network (1.52 vs 2.25), but the degree-preserving
null shifts further, so the z-score **flips sign**. The same true-positive pair is
"significantly proximal" on one interactome and "significantly distant" on another.

This is a single anecdote, not evidence — but it is exactly the failure mode direction D2 is
designed to quantify, and it suggests the network sweep should be run before, not after, the
headline benchmark. Overlap count is 0 in both cases, so the trivial baseline is at least
stable across networks here.
