#!/usr/bin/env python
"""
Experiment D2(b) -- interactome incompleteness, tested directly.

The hypothesis blames PPI incompleteness for confounding topological metrics.
Rather than only comparing networks that differ in many ways at once, this
degrades a single network in a controlled way: randomly delete a fraction f of
its edges, keep the largest connected component, and rescore a FIXED subsample
of drug-indication pairs.

  - if proximity carries real biological signal, AUROC should decay smoothly and
    the full network should be best
  - if the metric is unstable under incompleteness, AUROC should wander with no
    consistent relation to f

3 deletion seeds per fraction. Writes results/edge_removal.tsv.gz
"""
import argparse
import os
import sys
import time
from multiprocessing import Pool

import networkx as nx
import numpy as np
import pandas as pd

sys.path.insert(0, "src")
import nplib as N

OUT = "results"


def one_run(args):
    net, gene_source, frac, rep, n_pairs, n_random, seed = args
    G0 = N.load_network(net)
    rng = np.random.default_rng(1000 * rep + int(frac * 100) + seed)
    G = G0
    if frac > 0:
        edges = np.array(G0.edges())
        keep = rng.random(len(edges)) >= frac
        G = nx.Graph()
        G.add_nodes_from(G0.nodes())
        G.add_edges_from(map(tuple, edges[keep]))
        lcc = max(nx.connected_components(G), key=len)
        G = G.subgraph(lcc).copy()

    # the pair subsample is fixed across every (frac, rep) so the comparison is
    # like-for-like; pairs that lose all targets or all genes are dropped and
    # counted
    _, pairs_full, tmap0, gmap0, dkey = N.build_config_table(net, gene_source, G=G0)
    sub = pairs_full.sample(min(n_pairs, len(pairs_full)), random_state=N.SEED)
    nodes = set(G.nodes())
    tmap = {k: v & nodes for k, v in tmap0.items()}
    gmap = {k: v & nodes for k, v in gmap0.items()}
    keep_rows = sub[[bool(tmap.get(d)) and bool(gmap.get(s))
                     for d, s in zip(sub[dkey], sub.ind_id)]]
    if len(keep_rows) < 100:
        return None
    df = N.score_pairs(G, keep_rows, tmap, gmap, dkey, n_random=n_random,
                       verbose=False)
    rows = []
    for col in ["z_closest", "z_shortest", "z_kernel", "overlap",
                "d_closest", "n_targets", "n_genes"]:
        auc, ap, n = N.auroc_auprc(df.label.values, N.oriented(df, col))
        rows.append(dict(network=net, gene_source=gene_source, frac_removed=frac,
                         rep=rep, score=col, auroc=auc, auprc=ap, n=n,
                         nodes=G.number_of_nodes(), edges=G.number_of_edges(),
                         pairs_retained=len(keep_rows), pairs_attempted=len(sub)))
    N.log(f"  {net} f={frac} rep={rep}: {G.number_of_nodes()}n/{G.number_of_edges()}e "
          f"{len(keep_rows)} pairs, AUROC(z_closest)={rows[0]['auroc']:.3f}")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nets", default="biogrid,huri")
    ap.add_argument("--genes", default="all")
    ap.add_argument("--fracs", default="0,0.1,0.3,0.5,0.7,0.9")
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--n-pairs", type=int, default=1500)
    ap.add_argument("--n-random", type=int, default=300)
    ap.add_argument("--procs", type=int, default=14)
    a = ap.parse_args()
    N.set_seed()
    fracs = [float(x) for x in a.fracs.split(",")]
    jobs = [(net, a.genes, f, rep, a.n_pairs, a.n_random, N.SEED)
            for net in a.nets.split(",") for f in fracs
            for rep in range(1 if f == 0 else a.reps)]
    N.log(f"{len(jobs)} edge-removal runs")
    t0 = time.time()
    with Pool(min(a.procs, len(jobs))) as p:
        res = p.map(one_run, jobs)
    rows = [r for rs in res if rs for r in rs]
    df = pd.DataFrame(rows)
    os.makedirs(OUT, exist_ok=True)
    df.to_csv(f"{OUT}/edge_removal.tsv.gz", sep="\t", index=False)
    N.log(f"wrote {OUT}/edge_removal.tsv.gz ({len(df)} rows) in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
