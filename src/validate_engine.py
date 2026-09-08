#!/usr/bin/env python
"""
Pre-flight validation of the proximity engine, run before any experiment.

  V1  distance measures reproduce hand-computed values on a toy graph
  V2  degree-preserving randomisation preserves degree and returns distinct nodes
  V3  the vectorised block path agrees with a naive per-set loop (incl. centre)
  V4  degree-bias diagnostic on real interactomes -- Guney et al. report
      corr(z, target degree) ~ -0.01 vs corr(d, degree) ~ -0.46 on their network;
      reported here for every interactome, since it depends on network density
  V5  timing / module-size profile, so the full sweep can be budgeted

Writes results/validation.json.
"""
import json
import sys
import time

import networkx as nx
import numpy as np
import pandas as pd

sys.path.insert(0, "src")
import nplib as N

N.set_seed()
FAIL, REPORT = [], {}


def check(name, cond, detail=""):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name} {detail}")
    if not cond:
        FAIL.append(name)


# ---------------------------------------------------------------- V1 toy graph
print("V1 toy-graph distance measures")
# path 0-1-2-3-4 with a spur 2-5
Gt = nx.Graph([(0, 1), (1, 2), (2, 3), (3, 4), (2, 5)])
eng = N.ProximityEngine(Gt, bin_size=1, verbose=False)
T = eng.to_idx([0, 4])

S = eng.to_idx([5])                      # d(0,5)=3, d(4,5)=3
m = eng._block_measures(T[None, :], S[None, :], 1, None)
check("closest |S|=1", abs(m["closest"][0] - 3.0) < 1e-9, f"= {m['closest'][0]}")
check("shortest |S|=1", abs(m["shortest"][0] - 3.0) < 1e-9, f"= {m['shortest'][0]}")
# |S|=1: kernel = -mean(log(exp(-(d+1)))) = mean(d+1) = 4
check("kernel |S|=1", abs(m["kernel"][0] - 4.0) < 1e-5, f"= {m['kernel'][0]}")

S2 = eng.to_idx([3, 5])                  # d(0,*)=3,3 ; d(4,*)=1,3
m2 = eng._block_measures(T[None, :], S2[None, :], 2, None)
check("closest |S|=2", abs(m2["closest"][0] - (3 + 1) / 2) < 1e-9, f"= {m2['closest'][0]}")
check("shortest |S|=2", abs(m2["shortest"][0] - (3 + 3 + 1 + 3) / 4) < 1e-9,
      f"= {m2['shortest'][0]}")
kexp = -np.mean([np.log((np.exp(-4) + np.exp(-4)) / 2),
                 np.log((np.exp(-2) + np.exp(-4)) / 2)])
check("kernel |S|=2", abs(m2["kernel"][0] - kexp) < 1e-5, f"= {m2['kernel'][0]:.4f}")

# centre of S={1,3}: both have total intra-set distance 2 -> tie, both are centres
S3 = eng.to_idx([1, 3])
cm = eng.centre_mask(S3[None, :])
check("centre ties retained", cm.sum() == 2, f"(n centres = {cm.sum()})")
m3 = eng._block_measures(T[None, :], S3[None, :], 2, cm)
# mean over T of mean over centres: node0 -> d(0,1)=1,d(0,3)=3 ; node4 -> 3,1
check("centre distance", abs(m3["centre"][0] - 2.0) < 1e-6, f"= {m3['centre'][0]}")

# ---------------------------------------------------------------- V2/V3
print("V2/V3 randomisation and vectorisation on a real interactome (lit_bm)")
G = N.load_network("lit_bm")
eng = N.ProximityEngine(G, seed=7, verbose=False)
rng = np.random.default_rng(0)
gi = np.sort(rng.choice(eng.n_nodes, 25, replace=False)).astype(np.int32)
R = eng.random_sets(gi, 500)
check("degree-preserving draw",
      all((eng.node2bin[R[:, k]] == eng.node2bin[gi[k]]).all() for k in range(len(gi))))
n_dup = sum(len(np.unique(r)) < len(r) for r in R)
check("rejection sampling yields distinct nodes", n_dup / 500 < 0.02,
      f"({n_dup}/500 rows still hold a duplicate)")
check("mean degree preserved to <15%",
      abs(eng.deg[R].mean() - eng.deg[gi].mean()) / eng.deg[gi].mean() < 0.15,
      f"({eng.deg[R].mean():.2f} vs {eng.deg[gi].mean():.2f})")

S = np.sort(rng.choice(eng.n_nodes, 40, replace=False)).astype(np.int32)
Rs = eng.random_sets(S, 50)
Rt = R[:50]
cm = eng.centre_mask(Rs)
vec = eng._block_measures(Rt, Rs, len(S), cm)
naive = {k: [] for k in vec}
for i in range(50):
    sub = eng.D[np.ix_(Rt[i], Rs[i])].astype(np.float64)
    naive["closest"].append(sub.min(axis=1).mean())
    naive["shortest"].append(sub.mean())
    naive["kernel"].append(-np.mean(np.log(np.exp(-(sub + 1)).sum(axis=1) / len(S))))
    intra = eng.D[np.ix_(Rs[i], Rs[i])].astype(np.int64).sum(axis=1)
    centres = np.where(intra == intra.min())[0]
    naive["centre"].append(sub[:, centres].mean())
for k in vec:
    check(f"vectorised == naive ({k})", np.allclose(vec[k], naive[k], atol=1e-4))

# ---------------------------------------------------------------- V4 degree bias
print("V4 degree-bias diagnostic (400 sampled pairs, n_random=200, DisGeNET-all)")
rows = []
for net in ["biogrid", "string700", "hetionet_gig", "huri", "lit_bm"]:
    G, pairs, tmap, gmap, dkey = N.build_config_table(net, "all")
    sub = pairs.sample(min(400, len(pairs)), random_state=0)
    t0 = time.time()
    df = N.score_pairs(G, sub, tmap, gmap, dkey, n_random=200, verbose=False)
    el = time.time() - t0
    r = dict(
        network=net, nodes=G.number_of_nodes(), edges=G.number_of_edges(),
        mean_degree=round(2 * G.number_of_edges() / G.number_of_nodes(), 1),
        corr_z_deg=round(float(np.corrcoef(df.z_closest, df.mean_target_degree)[0, 1]), 3),
        corr_d_deg=round(float(np.corrcoef(df.d_closest, df.mean_target_degree)[0, 1]), 3),
        corr_z_nT=round(float(np.corrcoef(df.z_closest, df.n_targets)[0, 1]), 3),
        corr_d_nT=round(float(np.corrcoef(df.d_closest, df.n_targets)[0, 1]), 3),
        mean_d=round(float(df.d_closest.mean()), 2),
        ms_per_pair=round(1000 * el / len(sub), 1))
    rows.append(r)
    print("   ", r)
deg_tbl = pd.DataFrame(rows)
REPORT["degree_diagnostic"] = rows
# the correction must at least reduce the degree dependence on every network
ok = all(abs(r["corr_z_deg"]) <= abs(r["corr_d_deg"]) + 0.05 for r in rows)
check("z-score does not amplify degree dependence vs raw distance", ok)

# ---------------------------------------------------------------- V5 profile
print("V5 profile")
gmap_all = N.load_disease_genes("all")
sizes = pd.Series({k: len(v) for k, v in gmap_all.items()})
prof = {f"p{q}": int(sizes.quantile(q / 100)) for q in (25, 50, 75, 90, 99)}
prof["max"] = int(sizes.max())
print("    DisGeNET-all module sizes:", prof)
G, pairs, tmap, gmap, dkey = N.build_config_table("biogrid", "all")
cov = float((pairs.ind_id.map(lambda c: len(gmap.get(c, ()))) <= 600).mean())
print(f"    pairs with |S|<=600 (centre-measure coverage): {cov:.1%}")
REPORT["module_sizes"] = prof
REPORT["centre_coverage_biogrid_all"] = cov
REPORT["failures"] = FAIL

with open("results/validation.json", "w") as f:
    json.dump(REPORT, f, indent=2)
print("\nFAILURES:", FAIL if FAIL else "none")
sys.exit(1 if FAIL else 0)
