#!/usr/bin/env python
"""
Experiment D4 -- the negative-set ablation, the single most informative control.

Guney et al. (2016) evaluate proximity by treating EVERY unobserved drug-disease
pair as a negative. repoDB instead uses pairs that entered clinical trials and
failed. This script holds the positives, the network, the genes, the measure and
the code path completely fixed and swaps ONLY the negatives:

  hard    : repoDB clinical failures  (Terminated / Withdrawn / Suspended)
  random  : drug-disease pairs absent from repoDB entirely, sampled uniformly
            from the same drug and disease pools  (Guney-style)
  random_matched : as `random`, but the disease is sampled with the same
            marginal frequency as in the hard-negative set, so the two negative
            sets have the same disease-mix and differ only in plausibility

Writes results/scores_negatives/<net>__<gene>__<negset>.tsv.gz
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, "src")
import nplib as N

OUT = "results/scores_negatives"


def sample_random_negatives(pairs_all, pool_drugs, pool_dis, n, rng,
                            disease_weights=None):
    """Unobserved (drug, disease) pairs -- absent from repoDB under any status."""
    observed = set(zip(pairs_all.drug_id, pairs_all.ind_id))
    drugs = np.array(sorted(pool_drugs))
    dis = np.array(sorted(pool_dis))
    if disease_weights is not None:
        w = np.array([disease_weights.get(d, 0.0) for d in dis], dtype=float)
        w = w / w.sum()
    else:
        w = None
    out, tries = [], 0
    seen = set()
    while len(out) < n and tries < 200 * n:
        tries += 1
        d = drugs[rng.integers(0, len(drugs))]
        s = dis[rng.integers(0, len(dis))] if w is None else rng.choice(dis, p=w)
        if (d, s) in observed or (d, s) in seen:
            continue
        seen.add((d, s))
        out.append((d, s))
    return pd.DataFrame(out, columns=["drug_id", "ind_id"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", default="biogrid")
    ap.add_argument("--genes", default="all")
    ap.add_argument("--n-random", type=int, default=1000)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    N.set_seed()
    rng = np.random.default_rng(N.SEED)

    G, pairs, tmap, gmap, dkey = N.build_config_table(a.net, a.genes)
    pos = pairs[pairs.label == 1]
    hard = pairs[pairs.label == 0]
    N.log(f"{a.net} x {a.genes}: {len(pos)} approved, {len(hard)} clinical failures")

    all_pairs = N.load_pairs()
    dis_w = hard.ind_id.value_counts().to_dict()

    negsets = {"hard": hard}
    rnd = sample_random_negatives(all_pairs, set(pos[dkey]) | set(hard[dkey]),
                                  set(pairs.ind_id), len(hard), rng)
    rnd_m = sample_random_negatives(all_pairs, set(pos[dkey]) | set(hard[dkey]),
                                    set(pairs.ind_id), len(hard), rng,
                                    disease_weights=dis_w)
    meta = all_pairs.drop_duplicates("drug_id")[["drug_id", "drug_name"]]
    dis_meta = all_pairs.drop_duplicates("ind_id")[["ind_id", "ind_name", "sem_type"]]
    for nm, df in (("random", rnd), ("random_matched", rnd_m)):
        df = df.merge(meta, on="drug_id", how="left").merge(dis_meta, on="ind_id",
                                                            how="left")
        df["label"] = 0
        df["status"] = f"Unobserved({nm})"
        df["DetailedStatus"] = ""
        df["drug_name_lc"] = df.drug_name.str.lower()
        negsets[nm] = df
        N.log(f"  sampled {len(df)} {nm} negatives")

    eng = N.ProximityEngine(G, seed=N.SEED)
    # positives are scored once and reused across all three negative sets
    t0 = time.time()
    pos_path = f"{OUT}/{a.net}__{a.genes}__positives.tsv.gz"
    if os.path.exists(pos_path):
        dpos = pd.read_csv(pos_path, sep="\t")
    else:
        dpos = N.score_pairs(G, pos, tmap, gmap, dkey, n_random=a.n_random,
                             eng=eng, tag="pos")
        dpos.to_csv(pos_path, sep="\t", index=False)
    for nm, df in negsets.items():
        path = f"{OUT}/{a.net}__{a.genes}__{nm}.tsv.gz"
        if os.path.exists(path):
            N.log(f"  {nm}: cached")
            continue
        dneg = N.score_pairs(G, df, tmap, gmap, dkey, n_random=a.n_random,
                             eng=eng, tag=nm)
        out = pd.concat([dpos, dneg], ignore_index=True)
        out["negative_set"] = nm
        out.to_csv(path, sep="\t", index=False)
        N.log(f"  wrote {path}  ({len(dpos)}+ / {len(dneg)}-)")
    N.log(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
