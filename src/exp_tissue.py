#!/usr/bin/env python
"""
Experiment D3 -- tissue specificity.

Kitsak et al. (2016) argue disease modules are tissue-specific, and the
hypothesis names tissue-specificity as a confound of tissue-naive topology.
Test: restrict the interactome to genes expressed in the disease's tissue and
recompute proximity for exactly the same pairs.

Disease -> tissue assignment is data driven, so no manual curation enters:
each disease is assigned the GTEx v8 tissue with the highest mean
log2(TPM+1) over its own disease genes.

Three arms, all scored on the identical pair set:
  global   : the full interactome (tissue-naive) -- the reference
  tissue   : subgraph induced on genes with TPM >= `--tpm` in the assigned tissue
  shuffled : subgraph for a RANDOMLY assigned tissue -- the control that
             separates "tissue relevance" from "any sparsification helps"

Writes results/tissue_scores.tsv.gz and results/tissue_assignment.tsv
"""
import argparse
import gzip
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


# ---------------------------------------------------------------- GTEx loading
def load_gtex(entrez_ok: set):
    """(matrix DataFrame indexed by Entrez id, columns = 54 GTEx tissues).

    Ensembl -> Entrez via NCBI gene2ensembl; when several Ensembl ids map to the
    same Entrez gene the maximum median TPM per tissue is kept."""
    ens2ent = {}
    with gzip.open("datasets/mappings/gene2ensembl.gz", "rt") as f:
        hdr = f.readline().lstrip("#").rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(hdr)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[ix["tax_id"]] != "9606":
                continue
            ens2ent.setdefault(p[ix["Ensembl_gene_identifier"]], int(p[ix["GeneID"]]))
    df = pd.read_csv("datasets/tissue/GTEx_v8_gene_median_tpm.gct.gz", sep="\t",
                     skiprows=2)
    df["ens"] = df.Name.str.split(".").str[0]
    df["entrez"] = df.ens.map(ens2ent)
    df = df[df.entrez.notna()]
    df["entrez"] = df.entrez.astype(int)
    df = df[df.entrez.isin(entrez_ok)]
    tissues = [c for c in df.columns if c not in ("Name", "Description", "ens", "entrez")]
    mat = df.groupby("entrez")[tissues].max()
    return mat


def assign_tissues(mat: pd.DataFrame, gmap: dict, cuis) -> pd.DataFrame:
    """Each disease -> the tissue with the highest mean log2(TPM+1) of its genes."""
    lg = np.log2(mat + 1.0)
    rows = []
    for c in cuis:
        genes = [g for g in gmap[c] if g in lg.index]
        if not genes:
            continue
        prof = lg.loc[genes].mean(axis=0)
        rows.append(dict(ind_id=c, tissue=prof.idxmax(), score=float(prof.max()),
                         n_genes_with_expr=len(genes),
                         contrast=float(prof.max() - prof.mean())))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- worker
def one_tissue(args):
    (tissue, cuis_by_arm, net, gene_source, tpm, n_random) = args
    G0 = N.load_network(net)
    mat = MAT
    expressed = set(mat.index[(mat[tissue] >= tpm).values])
    sub = G0.subgraph([n for n in G0.nodes() if n in expressed])
    if sub.number_of_nodes() < 500:
        return None
    lcc = max(nx.connected_components(sub), key=len)
    G = G0.subgraph(lcc).copy()
    _, pairs_full, tmap0, gmap0, dkey = N.build_config_table(net, gene_source, G=G0)
    nodes = set(G.nodes())
    tmap = {k: v & nodes for k, v in tmap0.items()}
    gmap = {k: v & nodes for k, v in gmap0.items()}
    eng = N.ProximityEngine(G, seed=N.SEED, verbose=False)
    out = []
    for arm, cuis in cuis_by_arm.items():
        sel = pairs_full[pairs_full.ind_id.isin(cuis)]
        sel = sel[[bool(tmap.get(d)) and bool(gmap.get(s))
                   for d, s in zip(sel[dkey], sel.ind_id)]]
        if len(sel) == 0:
            continue
        df = N.score_pairs(G, sel, tmap, gmap, dkey, n_random=n_random, eng=eng,
                           verbose=False)
        df["arm"] = arm
        df["tissue"] = tissue
        df["tissue_nodes"] = G.number_of_nodes()
        df["tissue_edges"] = G.number_of_edges()
        out.append(df)
    N.log(f"  {tissue}: {G.number_of_nodes()}n/{G.number_of_edges()}e, "
          f"{sum(len(o) for o in out)} pair-scores")
    return pd.concat(out, ignore_index=True) if out else None


MAT = None


def init_worker(mat):
    global MAT
    MAT = mat


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", default="biogrid")
    ap.add_argument("--genes", default="all")
    ap.add_argument("--tpm", type=float, default=1.0)
    ap.add_argument("--n-random", type=int, default=500)
    ap.add_argument("--procs", type=int, default=12)
    a = ap.parse_args()
    N.set_seed()
    t0 = time.time()

    G0, pairs, tmap, gmap, dkey = N.build_config_table(a.net, a.genes)
    mat = load_gtex(set(G0.nodes()))
    N.log(f"GTEx: {mat.shape[0]} genes x {mat.shape[1]} tissues "
          f"(intersected with {a.net})")

    assign = assign_tissues(mat, gmap, pairs.ind_id.unique())
    rng = np.random.default_rng(N.SEED)
    tissues = sorted(mat.columns)
    assign["tissue_shuffled"] = rng.choice(tissues, len(assign))
    assign.to_csv(f"{OUT}/tissue_assignment.tsv", sep="\t", index=False)
    N.log(f"assigned tissues for {len(assign)} diseases; "
          f"{assign.tissue.nunique()} distinct")

    jobs = []
    for t in sorted(set(assign.tissue) | set(assign.tissue_shuffled)):
        arms = {}
        real = assign.ind_id[assign.tissue == t].tolist()
        shuf = assign.ind_id[assign.tissue_shuffled == t].tolist()
        if real:
            arms["tissue"] = real
        if shuf:
            arms["shuffled"] = shuf
        if arms:
            jobs.append((t, arms, a.net, a.genes, a.tpm, a.n_random))
    N.log(f"{len(jobs)} tissue subgraphs to build")

    with Pool(min(a.procs, len(jobs)), initializer=init_worker,
              initargs=(mat,)) as p:
        res = [r for r in p.map(one_tissue, jobs) if r is not None]
    tis = pd.concat(res, ignore_index=True)

    # the global (tissue-naive) arm, on the union of all pairs that survived
    keep = set(zip(tis.drug_key, tis.ind_id))
    sel = pairs[[(d, s) in keep for d, s in zip(pairs[dkey], pairs.ind_id)]]
    N.log(f"global arm: {len(sel)} pairs")
    glob = N.score_pairs(G0, sel, tmap, gmap, dkey, n_random=a.n_random,
                         tag="global")
    glob["arm"] = "global"
    glob["tissue"] = "ALL"
    glob["tissue_nodes"] = G0.number_of_nodes()
    glob["tissue_edges"] = G0.number_of_edges()

    out = pd.concat([tis, glob], ignore_index=True)
    out["network"] = a.net
    out["gene_source"] = a.genes
    out.to_csv(f"{OUT}/tissue_scores.tsv.gz", sep="\t", index=False)
    N.log(f"wrote {OUT}/tissue_scores.tsv.gz ({len(out)} rows) "
          f"in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
