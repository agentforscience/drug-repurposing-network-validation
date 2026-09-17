#!/usr/bin/env python
"""
Python 3 port of the Guney et al. (2016) network-based drug-disease proximity.
Reference implementation (Python 2): code/emreg00_toolbox/wrappers.py :: calculate_proximity

Faithful to the paper:
  - distance measures: closest / shortest / kernel / centre  (eqs 1-4)
  - reference distribution from degree-preserving random node sets
  - degree binning with >=100 nodes per bin
  - z = (d - mu_rand) / sigma_rand,  n_random = 1000

Acceleration: precompute the all-pairs shortest-path matrix once per network
(unweighted BFS, uint8), so every proximity evaluation is pure numpy indexing.
"""
import numpy as np, networkx as nx, scipy.sparse as sp
from scipy.sparse.csgraph import shortest_path

UNREACH = 255  # uint8 sentinel


# ------------------------------------------------------------------ APSP
def build_apsp(G):
    """Return (nodes, index_map, D) with D[i,j] = hop distance, uint8."""
    nodes = sorted(G.nodes())
    idx = {n: i for i, n in enumerate(nodes)}
    A = nx.to_scipy_sparse_array(G, nodelist=nodes, format='csr', dtype=np.int8)
    D = shortest_path(A, method='D', unweighted=True, directed=False)
    D[np.isinf(D)] = UNREACH
    return nodes, idx, D.astype(np.uint8)


# ------------------------------------------------------------------ degree binning
def get_degree_binning(G, bin_size=100):
    """Guney et al. binning: consecutive degrees merged until bin has >=bin_size nodes."""
    deg2nodes = {}
    for n, d in G.degree():
        deg2nodes.setdefault(d, []).append(n)
    values = sorted(deg2nodes)
    bins, i = [], 0
    while i < len(values):
        low = values[i]
        val = list(deg2nodes[values[i]])
        while len(val) < bin_size:
            i += 1
            if i == len(values): break
            val.extend(deg2nodes[values[i]])
        if i == len(values): i -= 1
        high = values[i]; i += 1
        if len(val) < bin_size and bins:
            low_, high_, val_ = bins[-1]
            bins[-1] = (low_, high, val_ + val)
        else:
            bins.append((low, high, val))
    return bins


def bin_lookup(G, bins):
    """node -> array of degree-equivalent node indices (as node ids)."""
    node2bin = {}
    for bi, (lo, hi, nodes) in enumerate(bins):
        for n in nodes:
            node2bin[n] = bi
    return node2bin


# ------------------------------------------------------------------ distances
def d_closest(D, T, S):
    """eq 1: mean over targets of min distance to any disease gene."""
    return D[np.ix_(T, S)].min(axis=1).mean()

def d_shortest(D, T, S):
    """eq 2: mean over all target x disease-gene pairs."""
    return D[np.ix_(T, S)].mean()

def d_kernel(D, T, S):
    """eq 3: exponentially down-weighted paths."""
    sub = D[np.ix_(T, S)].astype(np.float64)
    return -np.mean(np.log(np.sum(np.exp(-(sub + 1)), axis=1) / len(S)))

def d_centre(D, S_all, T, S):
    """eq 4: distance from targets to the topological centre of S."""
    sub = D[np.ix_(S, S)].astype(np.float64).sum(axis=1)
    centres = np.array(S)[sub == sub.min()]
    return D[np.ix_(T, centres)].mean()

DIST_FUNCS = {'closest': d_closest, 'shortest': d_shortest, 'kernel': d_kernel}


# ------------------------------------------------------------------ proximity
class ProximityCalculator:
    def __init__(self, G, bin_size=100, n_random=1000, seed=452456):
        self.G = G
        self.nodes, self.idx, self.D = build_apsp(G)
        self.bins = get_degree_binning(G, bin_size)
        self.node2bin = bin_lookup(G, self.bins)
        # bin -> array of node *indices*
        self.bin_members = [np.array([self.idx[n] for n in nodes], dtype=np.int32)
                            for _, _, nodes in self.bins]
        self.n_random = n_random
        self.rng = np.random.default_rng(seed)
        self._cache = {}

    def to_idx(self, genes):
        return np.array([self.idx[g] for g in genes if g in self.idx], dtype=np.int32)

    def random_sets(self, gene_idx, n=None):
        """n degree-preserving random node-index sets matching gene_idx."""
        n = n or self.n_random
        binids = [self.node2bin[self.nodes[i]] for i in gene_idx]
        out = np.empty((n, len(gene_idx)), dtype=np.int32)
        for k, b in enumerate(binids):
            members = self.bin_members[b]
            out[:, k] = members[self.rng.integers(0, len(members), size=n)]
        return out

    def proximity(self, targets, disease_genes, distance='closest', n_random=None):
        """Returns dict(d=, z=, mu=, sigma=, n_targets=, n_genes=) or None."""
        T = self.to_idx(targets); S = self.to_idx(disease_genes)
        if len(T) == 0 or len(S) == 0: return None
        n_random = n_random or self.n_random
        f = DIST_FUNCS[distance]
        d = f(self.D, T, S)
        Trand = self.random_sets(T, n_random)
        Srand = self.random_sets(S, n_random)
        vals = np.empty(n_random)
        for i in range(n_random):
            vals[i] = f(self.D, Trand[i], Srand[i])
        mu, sd = vals.mean(), vals.std()
        z = 0.0 if sd == 0 else (d - mu) / sd
        return dict(d=float(d), z=float(z), mu=float(mu), sigma=float(sd),
                    n_targets=int(len(T)), n_genes=int(len(S)))


# ------------------------------------------------------------------ baselines
def overlap_count(targets, disease_genes):
    """The 'drug-target overlap count' baseline named in the hypothesis."""
    return len(set(targets) & set(disease_genes))

def jaccard(targets, disease_genes):
    a, b = set(targets), set(disease_genes)
    return len(a & b) / len(a | b) if a | b else 0.0
