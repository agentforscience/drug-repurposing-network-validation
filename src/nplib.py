#!/usr/bin/env python
"""
nplib -- core library for the network-proximity / clinical-repurposing benchmark.

Contents
--------
1. Data loading      : interactomes, drug targets, disease genes, repoDB pairs
2. ProximityEngine   : Guney et al. (2016) network proximity -- all four distance
                       measures, degree-preserving z-scores, fully vectorised
                       over the random reference sets
3. Baselines         : the deliberately-trivial comparators named in the hypothesis
4. Evaluation        : AUROC / AUPRC, clustered bootstrap CIs, fast DeLong test

Fidelity to the reference implementation (code/emreg00_toolbox, Python 2):
  * distance measures closest / shortest / kernel / centre  (Guney eqs 1-4)
  * centre = mean distance to *all* tied minimum-eccentricity seed nodes
  * degree binning: consecutive degrees merged until a bin holds >= 100 nodes
  * random node sets drawn from the same degree bin, with up to 20 rejection
    retries to keep the drawn nodes distinct
  * z = (d - mu_rand) / sigma_rand, n_random = 1000, z := 0 when sigma = 0
Documented deviation: the reference collapses its random set with `set()` after
the retries, so a random set can end up SMALLER than the observed set (which
biases the null for large seed sets). Here the set size is held fixed, retaining
any residual duplicate. See REPORT.md.

All randomness flows through explicit seeds so every number in results/ is
reproducible.
"""
from __future__ import annotations

import time

import networkx as nx
import numpy as np
import pandas as pd
import scipy.stats as st
from scipy.sparse.csgraph import shortest_path
from sklearn.metrics import average_precision_score, roc_auc_score

SEED = 42
PROC = "datasets/processed"
UNREACH = 255  # uint8 sentinel; never occurs inside an LCC


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


def set_seed(seed: int = SEED):
    import random

    random.seed(seed)
    np.random.seed(seed)


# ====================================================================== loading

NETWORKS = [
    "biogrid",        # 976k edges, literature-curated  -> study-biased, dense
    "biogrid_lt",     # low-throughput subset           -> most study-biased
    "lit_bm",         # literature binary               -> study-biased, sparse
    "huri",           # systematic Y2H                  -> study-bias-free
    "hi_union",       # union of systematic screens     -> study-bias-free
    "string700",      # multi-evidence, score >= 700
    "hetionet_gig",   # the PPI used by Rephetio
]
GENE_SOURCES = ["all", "curated"]  # DisGeNET


def load_network(name: str) -> nx.Graph:
    """Load an interactome edge list (Entrez ids, already reduced to its LCC)."""
    df = pd.read_csv(f"{PROC}/ppi_{name}.tsv", sep="\t")
    G = nx.Graph()
    G.add_edges_from(zip(df.gene_a.values, df.gene_b.values))
    return G


def load_drug_targets(source: str = "hetionet"):
    """(dict drug key -> set(entrez), key column name).

    hetionet is keyed by DrugBank id, drugcentral by lower-cased drug name;
    repoDB carries both keys."""
    if source == "hetionet":
        df = pd.read_csv(f"{PROC}/drug_targets_hetionet.tsv", sep="\t")
        key = "drug_id"
    elif source == "drugcentral":
        df = pd.read_csv(f"{PROC}/drug_targets_drugcentral.tsv", sep="\t")
        key = "drug_name_lc"
    else:
        raise ValueError(source)
    return {k: set(g.entrez.values) for k, g in df.groupby(key)}, key


def load_disease_genes(source: str = "all", min_score: float = 0.0) -> dict:
    """UMLS CUI -> set(entrez)."""
    df = pd.read_csv(f"{PROC}/disease_genes_disgenet_{source}.tsv", sep="\t")
    if min_score > 0:
        df = df[df.score >= min_score]
    return {k: set(g.entrez.values) for k, g in df.groupby("cui")}


def load_pairs() -> pd.DataFrame:
    """repoDB drug-indication pairs. label = 1 approved, 0 clinically failed."""
    return pd.read_csv(f"{PROC}/pairs.tsv", sep="\t")


# --- repoDB failure-reason taxonomy (planning.md confound #4) ----------------
# Order matters: safety > efficacy > administrative, since a free-text reason can
# mention several. ~30% of repoDB failures are administrative, not pharmacological.
SAFETY_PATTERNS = [
    "safety", "adverse", "toxic", "side effect", "death", "mortality",
    "serious event", " sae", "tolerab",
]
EFFICACY_PATTERNS = [
    "efficac", "futility", "futile", "lack of effect", "no benefit",
    "did not meet", "failed to meet", "interim analysis", "ineffective",
    "lack of response", "no significant", "primary endpoint", "lack of benefit",
    "not effective", "no effect",
]
ADMIN_PATTERNS = [
    "fund", "accrual", "recruit", "enroll", "business", "sponsor", "financial",
    "administrat", "budget", "supply", "logistic", "personnel", "investigator left",
    "slow", "competing", "study never", "never activated", "not activated",
    "resources", "strategic", "company decision", "portfolio", "priorit",
    "merger", "no participants", "insufficient particip", "lack of particip",
    "written up", "staff", "feasibil", "revisions needed", "documentation",
]


def classify_failure(detailed) -> str:
    """Coarse taxonomy of a repoDB failure reason from its free-text field."""
    if not isinstance(detailed, str) or not detailed.strip():
        return "unstated"
    s = " " + detailed.lower()
    if any(p in s for p in SAFETY_PATTERNS):
        return "safety"
    if any(p in s for p in EFFICACY_PATTERNS):
        return "efficacy"
    if any(p in s for p in ADMIN_PATTERNS):
        return "administrative"
    return "other"


# ============================================================ proximity engine


class ProximityEngine:
    """Guney et al. (2016) network-based proximity.

    d_measure(T, S) is the raw distance between drug-target set T and disease-gene
    set S; z = (d - mu_rand)/sigma_rand with the reference distribution built from
    `n_random` degree-preserving random node sets of the same sizes.

    The APSP hop-distance matrix is precomputed once per network (uint8) so a
    proximity evaluation is a pure numpy gather + reduce; all `n_random` reference
    sets are evaluated in vectorised blocks rather than a Python loop.
    """

    MAX_BLOCK = 60_000_000  # ceiling on a gathered (chunk, |T|, |S|) block

    def __init__(self, G: nx.Graph, bin_size: int = 100, seed: int = SEED,
                 verbose: bool = True):
        self.G = G
        self.nodes = sorted(G.nodes())
        self.idx = {n: i for i, n in enumerate(self.nodes)}
        t0 = time.time()
        A = nx.to_scipy_sparse_array(G, nodelist=self.nodes, format="csr",
                                     dtype=np.int8)
        D = shortest_path(A, method="D", unweighted=True, directed=False)
        D[np.isinf(D)] = UNREACH
        self.D = np.ascontiguousarray(D.astype(np.uint8))
        self.n_nodes = len(self.nodes)
        if verbose:
            log(f"  APSP {self.n_nodes} nodes, {G.number_of_edges()} edges in "
                f"{time.time() - t0:.1f}s ({self.D.nbytes / 1e6:.0f} MB)")
        self.deg = np.array([G.degree(n) for n in self.nodes], dtype=np.int32)
        self.bins = self._degree_binning(bin_size)
        self.node2bin = np.empty(self.n_nodes, dtype=np.int32)
        for bi, members in enumerate(self.bins):
            self.node2bin[members] = bi
        # flat (concatenated) bin membership, so a degree-matched draw for an
        # arbitrary vector of bin ids is a single vectorised gather
        self.bin_flat = np.concatenate(self.bins).astype(np.int32)
        self.bin_len = np.array([len(b) for b in self.bins], dtype=np.int64)
        self.bin_off = np.concatenate([[0], np.cumsum(self.bin_len)[:-1]])
        self.rng = np.random.default_rng(seed)
        # exp(-(d+1)) lookup so the kernel measure never materialises a float
        # copy of the whole distance block
        self.KLUT = np.exp(-(np.arange(256, dtype=np.float32) + 1.0))
        self._T_cache: dict = {}

    # ------------------------------------------------------------ degree bins
    def _degree_binning(self, bin_size: int):
        """Consecutive degree values merged until each bin holds >= bin_size nodes
        (same rule as network_utilities.get_degree_binning)."""
        deg2nodes: dict[int, list] = {}
        for i, d in enumerate(self.deg):
            deg2nodes.setdefault(int(d), []).append(i)
        values = sorted(deg2nodes)
        bins, i = [], 0
        while i < len(values):
            val = list(deg2nodes[values[i]])
            while len(val) < bin_size:
                i += 1
                if i == len(values):
                    break
                val.extend(deg2nodes[values[i]])
            if i == len(values):
                i -= 1
            i += 1
            if len(val) < bin_size and bins:
                bins[-1] = bins[-1] + val
            else:
                bins.append(val)
        return [np.array(b, dtype=np.int32) for b in bins]

    # ------------------------------------------------------------ helpers
    def to_idx(self, genes) -> np.ndarray:
        """Gene ids -> node indices, dropping genes absent from this interactome."""
        return np.array(sorted(self.idx[g] for g in genes if g in self.idx),
                        dtype=np.int32)

    #: sets larger than this skip duplicate-rejection (see `random_sets`)
    DEDUP_MAX_K = 500

    def _draw(self, binids) -> np.ndarray:
        """Vectorised degree-matched draw: one random member of bin b per entry."""
        b = np.asarray(binids)
        off = self.bin_off[b] + (self.rng.random(b.shape) * self.bin_len[b]
                                 ).astype(np.int64)
        return self.bin_flat[off]

    def random_sets(self, gene_idx: np.ndarray, n_random: int,
                    retries: int = 20) -> np.ndarray:
        """(n_random, |gene_idx|) degree-preserving random node-index sets.

        Each slot is drawn from its node's degree bin. Up to `retries` vectorised
        rejection rounds then remove within-row duplicates, matching the reference
        implementation's behaviour.

        Duplicate rejection is skipped for sets larger than DEDUP_MAX_K: bins hold
        >= 100 nodes, so a several-thousand-gene module cannot be made distinct at
        all (the reference collapses such a set with `set()`, shrinking it far
        below the observed size, which biases the null). At that size a residual
        duplicate shifts a mean- or min-based statistic negligibly, whereas the
        rejection loop would dominate runtime.
        """
        binids = self.node2bin[gene_idx]
        k = len(gene_idx)
        out = self._draw(np.broadcast_to(binids, (n_random, k))).astype(np.int32)
        if 1 < k <= self.DEDUP_MAX_K:
            prev = None
            for _ in range(retries):
                order = np.argsort(out, axis=1, kind="stable")
                srt = np.take_along_axis(out, order, axis=1)
                dup_sorted = np.zeros_like(srt, dtype=bool)
                dup_sorted[:, 1:] = srt[:, 1:] == srt[:, :-1]
                n_dup = int(dup_sorted.sum())
                # stop when clean, or when rejection has stopped making progress
                # (a bin can be too small to supply distinct nodes for every slot)
                if n_dup == 0 or (prev is not None and n_dup > 0.9 * prev):
                    break
                prev = n_dup
                dup = np.zeros_like(dup_sorted)
                np.put_along_axis(dup, order, dup_sorted, axis=1)
                rows, cols = np.nonzero(dup)
                out[rows, cols] = self._draw(binids[cols])
        return out

    def targets_random(self, key, T: np.ndarray, n_random: int) -> np.ndarray:
        """Cached random draws for a drug (same drug recurs across indications)."""
        ck = (key, len(T), n_random)
        if ck not in self._T_cache:
            self._T_cache[ck] = self.random_sets(T, n_random)
        return self._T_cache[ck]

    def clear_cache(self):
        self._T_cache.clear()

    # ------------------------------------------------------------ measures
    def _block_measures(self, Tr: np.ndarray, Sr: np.ndarray, n_S: int,
                        centre_mask: np.ndarray | None):
        """All measures for a batch of (target-set, seed-set) index pairs.

        Tr : (n, |T|) node indices    Sr : (n, |S|) node indices
        centre_mask : (n, |S|) bool marking the topological centre(s) of each Sr,
                      or None to skip the centre measure.
        Returns dict measure -> (n,) float array.
        """
        n, nT = Tr.shape
        chunk = max(1, min(n, self.MAX_BLOCK // max(1, nT * n_S)))
        acc = {k: [] for k in ("closest", "shortest", "kernel")}
        if centre_mask is not None:
            acc["centre"] = []
        for lo in range(0, n, chunk):
            hi = min(lo + chunk, n)
            blk = self.D[Tr[lo:hi][:, :, None], Sr[lo:hi][:, None, :]]  # uint8
            acc["closest"].append(blk.min(axis=2).mean(axis=1, dtype=np.float64))
            acc["shortest"].append(
                blk.sum(axis=(1, 2), dtype=np.int64) / (nT * n_S))
            acc["kernel"].append(
                -np.log(self.KLUT[blk].sum(axis=2) / n_S).mean(axis=1))
            if centre_mask is not None:
                cm = centre_mask[lo:hi]
                num = np.einsum("ijk,ik->i", blk.astype(np.float32),
                                cm.astype(np.float32))
                acc["centre"].append(num / (nT * cm.sum(axis=1)))
        return {k: np.concatenate(v) for k, v in acc.items()}

    def centre_mask(self, Sr: np.ndarray, max_S: int = 600):
        """(n, |S|) bool: the minimum-total-distance node(s) of each seed set.

        Ties are all retained and averaged over, as in get_center_of_subnetwork.
        Returns None when |S| exceeds `max_S` (the O(|S|^2) block is skipped)."""
        n, nS = Sr.shape
        if nS > max_S:
            return None
        chunk = max(1, min(n, self.MAX_BLOCK // max(1, nS * nS)))
        out = np.empty((n, nS), dtype=bool)
        for lo in range(0, n, chunk):
            hi = min(lo + chunk, n)
            blk = self.D[Sr[lo:hi][:, :, None], Sr[lo:hi][:, None, :]]
            sums = blk.sum(axis=2, dtype=np.int64)
            out[lo:hi] = sums == sums.min(axis=1, keepdims=True)
        return out

    # ------------------------------------------------------------ main entry
    def proximity(self, T: np.ndarray, S: np.ndarray, n_random: int = 1000,
                  drug_key=None, Sr: np.ndarray | None = None,
                  S_centre_mask=None, Sr_centre_mask=None) -> dict:
        """Raw distances + z-scores for every measure.

        Sr / Sr_centre_mask / S_centre_mask may be supplied by the caller so that
        the per-disease work (random seed sets and their centres) is done once and
        reused across all drugs tested against that disease.
        """
        nT, nS = len(T), len(S)
        if nT == 0 or nS == 0:
            return {}
        if Sr is None:
            Sr = self.random_sets(S, n_random)
            Sr_centre_mask = self.centre_mask(Sr)
            S_centre_mask = None if Sr_centre_mask is None else \
                self.centre_mask(S[None, :])
        Tr = self.targets_random(drug_key, T, n_random) if drug_key is not None \
            else self.random_sets(T, n_random)

        obs = self._block_measures(T[None, :], S[None, :], nS, S_centre_mask)
        rnd = self._block_measures(Tr, Sr, nS, Sr_centre_mask)

        out = {}
        for k in ("closest", "shortest", "kernel", "centre"):
            if k not in obs:
                out[f"d_{k}"] = np.nan
                out[f"z_{k}"] = np.nan
                continue
            d = float(obs[k][0])
            mu, sd = float(rnd[k].mean()), float(rnd[k].std())
            out[f"d_{k}"] = d
            out[f"z_{k}"] = 0.0 if sd == 0 else (d - mu) / sd
            out[f"mu_{k}"] = mu
            out[f"sd_{k}"] = sd
        out["n_targets"] = nT
        out["n_genes"] = nS
        out["mean_target_degree"] = float(self.deg[T].mean())
        out["mean_gene_degree"] = float(self.deg[S].mean())
        return out


# ================================================================== baselines


def baseline_features(T_genes: set, S_genes: set) -> dict:
    """The trivial comparators. `overlap` is the one the hypothesis names."""
    inter = len(T_genes & S_genes)
    union = len(T_genes | S_genes)
    return {
        "overlap": inter,
        "overlap_bin": int(inter > 0),
        "jaccard": inter / union if union else 0.0,
        "overlap_frac_targets": inter / len(T_genes) if T_genes else 0.0,
    }


#: score column -> +1 if larger means "more likely approved", -1 if smaller does.
SCORE_DIRECTION = {
    "z_closest": -1, "z_shortest": -1, "z_kernel": -1, "z_centre": -1,
    "d_closest": -1, "d_shortest": -1, "d_kernel": -1, "d_centre": -1,
    "overlap": +1, "overlap_bin": +1, "jaccard": +1, "overlap_frac_targets": +1,
    "n_targets": +1, "n_genes": +1,
    "mean_target_degree": +1, "mean_gene_degree": +1,
}
PROXIMITY_SCORES = ["z_closest", "z_shortest", "z_kernel", "z_centre",
                    "d_closest", "d_shortest", "d_kernel", "d_centre"]
BASELINE_SCORES = ["overlap", "overlap_bin", "jaccard", "overlap_frac_targets",
                   "n_targets", "n_genes", "mean_target_degree",
                   "mean_gene_degree"]
PRETTY = {
    "z_closest": "proximity z (closest)", "z_shortest": "proximity z (shortest)",
    "z_kernel": "proximity z (kernel)", "z_centre": "proximity z (centre)",
    "d_closest": "raw distance (closest)", "d_shortest": "raw distance (shortest)",
    "d_kernel": "raw distance (kernel)", "d_centre": "raw distance (centre)",
    "overlap": "target/disease-gene overlap count", "overlap_bin": "any overlap (0/1)",
    "jaccard": "Jaccard(targets, disease genes)",
    "overlap_frac_targets": "fraction of targets in module",
    "n_targets": "number of drug targets", "n_genes": "disease module size",
    "mean_target_degree": "mean target degree", "mean_gene_degree": "mean module degree",
}


def oriented(df: pd.DataFrame, col: str) -> np.ndarray:
    """Score vector oriented so that higher = predicted approved."""
    return SCORE_DIRECTION.get(col, 1) * df[col].values.astype(float)


# ================================================================ evaluation


def auroc_auprc(y, s):
    y, s = np.asarray(y), np.asarray(s, dtype=float)
    ok = np.isfinite(s)
    if ok.sum() < 20 or len(np.unique(y[ok])) < 2:
        return np.nan, np.nan, int(ok.sum())
    return (float(roc_auc_score(y[ok], s[ok])),
            float(average_precision_score(y[ok], s[ok])),
            int(ok.sum()))


def bootstrap_auc(y, s, groups=None, n_boot=2000, seed=SEED):
    """Percentile bootstrap CI for AUROC.

    groups=None  -> resample pairs.
    groups=array -> cluster bootstrap resampling whole groups (drug or disease),
                    because repoDB pairs are not independent (one drug spans many
                    indications).
    """
    y, s = np.asarray(y), np.asarray(s, dtype=float)
    ok = np.isfinite(s)
    y, s = y[ok], s[ok]
    if groups is not None:
        groups = np.asarray(groups)[ok]
    if len(y) < 20 or len(np.unique(y)) < 2:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    vals = []
    if groups is None:
        n = len(y)
        for _ in range(n_boot):
            i = rng.integers(0, n, n)
            if len(np.unique(y[i])) < 2:
                continue
            vals.append(roc_auc_score(y[i], s[i]))
    else:
        uniq = np.unique(groups)
        member = {g: np.where(groups == g)[0] for g in uniq}
        for _ in range(n_boot):
            gs = rng.choice(uniq, len(uniq), replace=True)
            i = np.concatenate([member[g] for g in gs])
            if len(np.unique(y[i])) < 2:
                continue
            vals.append(roc_auc_score(y[i], s[i]))
    if not vals:
        return np.nan, np.nan
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


# ---------------------------------------------------------------- fast DeLong
def _midrank(x):
    J = np.argsort(x)
    Z = x[J]
    N = len(x)
    T = np.zeros(N, dtype=float)
    i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    T2 = np.empty(N, dtype=float)
    T2[J] = T
    return T2


def _fast_delong(preds_sorted, m):
    """DeLong (1988) covariance via Sun & Xu (2014).
    preds_sorted: (k, n) with the first m columns the positives."""
    n = preds_sorted.shape[1] - m
    k = preds_sorted.shape[0]
    pos, neg = preds_sorted[:, :m], preds_sorted[:, m:]
    tx = np.empty([k, m]); ty = np.empty([k, n]); tz = np.empty([k, m + n])
    for r in range(k):
        tx[r] = _midrank(pos[r])
        ty[r] = _midrank(neg[r])
        tz[r] = _midrank(preds_sorted[r])
    aucs = tz[:, :m].sum(axis=1) / m / n - (m + 1.0) / 2.0 / n
    v01 = (tz[:, :m] - tx) / n
    v10 = 1.0 - (tz[:, m:] - ty) / m
    cov = np.cov(v01).reshape(k, k) / m + np.cov(v10).reshape(k, k) / n
    return aucs, cov


def delong_test(y, s1, s2):
    """Two-sided DeLong test for two correlated AUROCs.

    Returns (auc1, auc2, diff, p, se). Non-finite values in either score drop the
    pair from BOTH, keeping the comparison paired."""
    s1, s2 = np.asarray(s1, dtype=float), np.asarray(s2, dtype=float)
    ok = np.isfinite(s1) & np.isfinite(s2)
    y = np.asarray(y)[ok].astype(int)
    preds = np.vstack([s1[ok], s2[ok]])
    order = np.argsort(-y, kind="mergesort")
    y, preds = y[order], preds[:, order]
    m = int(y.sum())
    if m < 5 or len(y) - m < 5:
        return (np.nan,) * 5
    aucs, cov = _fast_delong(preds, m)
    l = np.array([[1.0, -1.0]])
    se = float(np.sqrt(max(float(l @ cov @ l.T), 0.0)))
    diff = float(aucs[0] - aucs[1])
    p = (1.0 if diff == 0 else 0.0) if se == 0 else \
        float(2 * st.norm.sf(abs(diff / se)))
    return float(aucs[0]), float(aucs[1]), diff, p, se


# ================================================================== assembly


def build_config_table(net_name: str, gene_source: str, target_source="hetionet",
                       pairs: pd.DataFrame | None = None, G=None):
    """(G, kept pairs, drug->targets, disease->genes, drug key column).

    A pair is kept only if the drug has >=1 target and the disease >=1 gene
    present in this interactome. Coverage differences between networks are a
    finding in their own right, so they are recorded, not hidden.
    """
    if G is None:
        G = load_network(net_name)
    nodes = set(G.nodes())
    tmap, dkey = load_drug_targets(target_source)
    gmap = load_disease_genes(gene_source)
    if pairs is None:
        pairs = load_pairs()
    tmap = {k: v for k, v in ((k, v & nodes) for k, v in tmap.items()) if v}
    gmap = {k: v for k, v in ((k, v & nodes) for k, v in gmap.items()) if v}
    keep = pairs[pairs[dkey].isin(tmap) & pairs.ind_id.isin(gmap)].copy()
    return G, keep, tmap, gmap, dkey


def score_pairs(G, pairs, tmap, gmap, dkey, n_random=1000, seed=SEED,
                verbose=True, tag="", eng=None, centre_max_S=600):
    """Score every pair: proximity (all measures) + baselines.

    The loop is grouped by disease so that the expensive per-disease work (random
    seed sets and their topological centres) is done once and shared across all
    drugs tested against that disease.
    """
    if eng is None:
        eng = ProximityEngine(G, seed=seed, verbose=verbose)
    rows = []
    t0 = time.time()
    done = 0
    for cui, grp in pairs.groupby("ind_id", sort=False):
        S = eng.to_idx(gmap[cui])
        if len(S) == 0:
            continue
        Sr = eng.random_sets(S, n_random)
        Sr_cm = eng.centre_mask(Sr, max_S=centre_max_S)
        S_cm = None if Sr_cm is None else eng.centre_mask(S[None, :],
                                                          max_S=centre_max_S)
        for _, r in grp.iterrows():
            dk = r[dkey]
            T = eng.to_idx(tmap[dk])
            res = eng.proximity(T, S, n_random=n_random, drug_key=dk,
                                Sr=Sr, S_centre_mask=S_cm, Sr_centre_mask=Sr_cm)
            if not res:
                continue
            res.update(baseline_features(tmap[dk], gmap[cui]))
            res["drug_key"] = dk
            res["drug_id"] = r.drug_id
            res["drug_name"] = r.drug_name
            res["ind_id"] = cui
            res["ind_name"] = r.ind_name
            res["label"] = int(r.label)
            res["status"] = r.status
            res["sem_type"] = r.sem_type
            res["DetailedStatus"] = r.get("DetailedStatus", "")
            rows.append(res)
            done += 1
            if verbose and done % 1000 == 0:
                el = time.time() - t0
                log(f"  {tag} {done}/{len(pairs)} pairs "
                    f"({el:.0f}s, {1000 * el / done:.1f} ms/pair)")
    df = pd.DataFrame(rows)
    if verbose:
        log(f"  {tag} scored {len(df)} pairs in {time.time() - t0:.0f}s")
    return df
