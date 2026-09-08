#!/usr/bin/env python
"""
Experiment D1 + D2(a): score every repoDB drug-indication pair on every
interactome x disease-gene-source configuration.

One worker process per interactome (the APSP matrix is the expensive object, and
both disease-gene sources reuse it). Writes one TSV per configuration to
results/scores/<net>__<gene_source>.tsv.gz plus a manifest.

Usage:  python src/exp_score.py [--n-random 1000] [--nets a,b] [--genes all,curated]
"""
import argparse
import json
import os
import sys
import time
from multiprocessing import Pool

import pandas as pd

sys.path.insert(0, "src")
import nplib as N

OUT = "results/scores"


def run_network(args):
    net, gene_sources, n_random, target_source = args
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    G = N.load_network(net)
    eng = N.ProximityEngine(G, seed=N.SEED, verbose=True)
    out = []
    for gs in gene_sources:
        path = f"{OUT}/{net}__{gs}__{target_source}.tsv.gz"
        if os.path.exists(path):
            N.log(f"{net} x {gs}: cached, skipping")
            out.append(path)
            continue
        _, pairs, tmap, gmap, dkey = N.build_config_table(
            net, gs, target_source=target_source, G=G)
        N.log(f"{net} x {gs} ({target_source}): {len(pairs)} pairs "
              f"({int(pairs.label.sum())}+/{int((1 - pairs.label).sum())}-)")
        eng.clear_cache()
        df = N.score_pairs(G, pairs, tmap, gmap, dkey, n_random=n_random,
                           eng=eng, tag=f"{net}x{gs}")
        df["network"] = net
        df["gene_source"] = gs
        df["target_source"] = target_source
        df.to_csv(path, sep="\t", index=False)
        N.log(f"{net} x {gs}: wrote {path}")
        out.append(path)
    N.log(f"{net} complete in {time.time() - t0:.0f}s")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-random", type=int, default=1000)
    ap.add_argument("--nets", default=",".join(N.NETWORKS))
    ap.add_argument("--genes", default="all,curated")
    ap.add_argument("--target-source", default="hetionet")
    ap.add_argument("--procs", type=int, default=7)
    a = ap.parse_args()

    nets = a.nets.split(",")
    genes = a.genes.split(",")
    N.set_seed()
    jobs = [(net, genes, a.n_random, a.target_source) for net in nets]
    t0 = time.time()
    with Pool(min(a.procs, len(jobs))) as p:
        paths = p.map(run_network, jobs)
    manifest = dict(n_random=a.n_random, nets=nets, gene_sources=genes,
                    target_source=a.target_source, seed=N.SEED,
                    files=[x for xs in paths for x in xs],
                    wall_seconds=round(time.time() - t0, 1))
    os.makedirs("results", exist_ok=True)
    key = f"manifest_scores_{a.target_source}.json"
    json.dump(manifest, open(f"results/{key}", "w"), indent=2)
    N.log(f"ALL DONE in {manifest['wall_seconds']}s -> results/{key}")


if __name__ == "__main__":
    main()
