#!/usr/bin/env python
"""
Analysis for all experiments. Produces results/*.csv, results/summary.json and
figures/*.png.

Sections
  A  coverage / label-bias check on the pairs lost to ID joins
  B  headline benchmark (D1): AUROC/AUPRC + clustered bootstrap CIs, every score
  C  C2 decision: DeLong proximity vs overlap count (and every other baseline)
  D  supervised comparator: does proximity add anything on top of the baselines?
  E  failure-reason and disease-type stratification (D5)
  F  negative-set ablation (D4)
  G  network sweep (D2a) and edge-removal curves (D2b)
  H  tissue specificity (D3)
  I  figures

Pre-declared decision rules (planning.md):
  C1 supported if the upper bootstrap CI bound for AUROC(proximity) < 0.65
  C2 supported if DeLong P > 0.05 OR |AUROC difference| < 0.02
"""
import glob
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, roc_curve

sys.path.insert(0, "src")
import nplib as N

RES, FIG = "results", "figures"
os.makedirs(FIG, exist_ok=True)

# validated categorical palette (see references/palette.md; validator run in log)
C_PROX, C_OVER, C_OTHER, C_ALT = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
C_GRID, C_INK, C_INK2 = "#d8d8d4", "#0b0b0b", "#52514e"
plt.rcParams.update({
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
    "axes.edgecolor": C_INK2, "axes.labelcolor": C_INK, "text.color": C_INK,
    "xtick.color": C_INK2, "ytick.color": C_INK2, "font.size": 9,
    "axes.grid": True, "grid.color": C_GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 170,
    "savefig.bbox": "tight", "legend.frameon": False,
})

HEAD_NET, HEAD_GENE = "biogrid", "all"
PRIMARY = "z_closest"          # Guney's best measure
BASELINE = "overlap"           # the comparator the hypothesis names
SUMMARY = {}


def load_scores(net=None, gene=None, target="hetionet"):
    pat = f"{RES}/scores/{net or '*'}__{gene or '*'}__{target}.tsv.gz"
    files = sorted(glob.glob(pat))
    if not files:
        return pd.DataFrame()
    return pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)


def eval_scores(df, scores, n_boot=2000, groups=True):
    """AUROC/AUPRC with pair, by-drug and by-disease bootstrap CIs."""
    rows = []
    y = df.label.values
    for col in scores:
        if col not in df.columns:
            continue
        s = N.oriented(df, col)
        auc, ap, n = N.auroc_auprc(y, s)
        r = dict(score=col, label=N.PRETTY.get(col, col), auroc=auc, auprc=ap, n=n)
        if groups and np.isfinite(auc):
            r["lo_pair"], r["hi_pair"] = N.bootstrap_auc(y, s, None, n_boot)
            r["lo_drug"], r["hi_drug"] = N.bootstrap_auc(y, s, df.drug_key.values,
                                                         n_boot)
            r["lo_dis"], r["hi_dis"] = N.bootstrap_auc(y, s, df.ind_id.values,
                                                       n_boot)
        rows.append(r)
    return pd.DataFrame(rows)


# ====================================================== A. coverage / label bias
def section_A():
    print("\n=== A. coverage and label-bias check ===")
    allp = N.load_pairs()
    rows = []
    for f in sorted(glob.glob(f"{RES}/scores/*__hetionet.tsv.gz")):
        d = pd.read_csv(f, sep="\t", usecols=["drug_id", "ind_id", "label",
                                              "network", "gene_source"])
        kept = set(zip(d.drug_id, d.ind_id))
        mask = [(a, b) in kept for a, b in zip(allp.drug_id, allp.ind_id)]
        inn, out = allp[mask], allp[[not m for m in mask]]
        # two-proportion z-test on the approved rate, kept vs dropped
        p1, p2 = inn.label.mean(), out.label.mean()
        n1, n2 = len(inn), len(out)
        pp = (inn.label.sum() + out.label.sum()) / (n1 + n2)
        se = np.sqrt(pp * (1 - pp) * (1 / n1 + 1 / n2))
        import scipy.stats as st
        z = (p1 - p2) / se
        rows.append(dict(network=d.network.iloc[0], gene_source=d.gene_source.iloc[0],
                         n_kept=n1, n_dropped=n2, frac_kept=round(n1 / len(allp), 3),
                         approved_kept=round(p1, 3), approved_dropped=round(p2, 3),
                         z=round(float(z), 2), p=float(2 * st.norm.sf(abs(z)))))
    cov = pd.DataFrame(rows).sort_values(["gene_source", "network"])
    cov.to_csv(f"{RES}/coverage_labelbias.csv", index=False)
    print(cov.to_string(index=False))
    SUMMARY["coverage"] = cov.to_dict("records")
    return cov


# ================================================ B/C. headline benchmark + DeLong
def section_BC():
    print("\n=== B. headline benchmark ===")
    df = load_scores(HEAD_NET, HEAD_GENE)
    assert len(df), "no headline scores found"
    print(f"{HEAD_NET} x DisGeNET-{HEAD_GENE}: {len(df)} pairs, "
          f"{df.label.sum()} approved / {(1 - df.label).sum()} failed, "
          f"{df.drug_key.nunique()} drugs, {df.ind_id.nunique()} diseases")
    tbl = eval_scores(df, N.PROXIMITY_SCORES + N.BASELINE_SCORES)
    tbl.insert(0, "network", HEAD_NET)
    tbl.insert(1, "gene_source", HEAD_GENE)
    tbl.to_csv(f"{RES}/headline_auc.csv", index=False)
    print(tbl.round(4).to_string(index=False))

    print("\n=== C. DeLong tests vs the trivial baselines ===")
    y = df.label.values
    rows = []
    prox_cols = [c for c in N.PROXIMITY_SCORES if df[c].notna().any()]
    for a in prox_cols:
        for b in N.BASELINE_SCORES:
            r1, r2, diff, p, se = N.delong_test(y, N.oriented(df, a),
                                                N.oriented(df, b))
            rows.append(dict(a=a, b=b, auc_a=r1, auc_b=r2, diff=diff, se=se, p=p))
    dl = pd.DataFrame(rows)
    # Benjamini-Hochberg across the whole family of comparisons
    from statsmodels.stats.multitest import multipletests
    ok = dl.p.notna()
    dl.loc[ok, "p_bh"] = multipletests(dl.p[ok], method="fdr_bh")[1]
    dl.to_csv(f"{RES}/delong_tests.csv", index=False)
    key = dl[(dl.a == PRIMARY) & (dl.b == BASELINE)].iloc[0]
    print(dl[dl.a == PRIMARY].round(4).to_string(index=False))

    # pre-declared decision rules
    h = tbl[tbl.score == PRIMARY].iloc[0]
    o = tbl[tbl.score == BASELINE].iloc[0]
    c1 = bool(h.hi_pair < 0.65)
    c1_drug = bool(h.hi_drug < 0.65)
    c2 = bool(key.p > 0.05 or abs(key["diff"]) < 0.02)
    SUMMARY["headline"] = dict(
        network=HEAD_NET, gene_source=HEAD_GENE, n_pairs=int(len(df)),
        n_pos=int(df.label.sum()), n_neg=int((1 - df.label).sum()),
        n_drugs=int(df.drug_key.nunique()), n_diseases=int(df.ind_id.nunique()),
        proximity_auroc=float(h.auroc),
        proximity_ci_pair=[float(h.lo_pair), float(h.hi_pair)],
        proximity_ci_drug=[float(h.lo_drug), float(h.hi_drug)],
        proximity_ci_disease=[float(h.lo_dis), float(h.hi_dis)],
        proximity_auprc=float(h.auprc), positive_rate=float(df.label.mean()),
        overlap_auroc=float(o.auroc), overlap_ci_pair=[float(o.lo_pair),
                                                       float(o.hi_pair)],
        delong_diff=float(key["diff"]), delong_p=float(key.p),
        C1_supported=c1, C1_supported_drug_cluster=c1_drug, C2_supported=c2)
    print(f"\n  C1 (AUROC < 0.65): upper CI = {h.hi_pair:.3f} (pair), "
          f"{h.hi_drug:.3f} (by-drug) -> {'SUPPORTED' if c1 else 'NOT supported'}")
    print(f"  C2 (no gain over overlap): diff = {key['diff']:+.4f}, "
          f"P = {key.p:.3g} -> {'SUPPORTED' if c2 else 'NOT supported'}")
    return df, tbl, dl


# ============================================ D. supervised incremental value
def section_D(df):
    print("\n=== D. supervised comparator (GroupKFold by drug and by disease) ===")
    sets = {
        "baselines only": N.BASELINE_SCORES,
        "proximity only": ["z_closest", "z_shortest", "z_kernel"],
        "proximity + baselines": ["z_closest", "z_shortest", "z_kernel"]
                                 + N.BASELINE_SCORES,
    }
    rows = []
    y = df.label.values
    for gname, groups in (("drug", df.drug_key.values), ("disease", df.ind_id.values)):
        for name, cols in sets.items():
            X = df[cols].values.astype(float)
            X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
            oof = np.full(len(y), np.nan)
            gkf = GroupKFold(n_splits=5)
            for tr, te in gkf.split(X, y, groups):
                m = make_pipeline(StandardScaler(),
                                  LogisticRegression(max_iter=2000, C=1.0))
                m.fit(X[tr], y[tr])
                oof[te] = m.predict_proba(X[te])[:, 1]
            auc = roc_auc_score(y, oof)
            lo, hi = N.bootstrap_auc(y, oof, groups, n_boot=1000)
            rows.append(dict(grouping=gname, feature_set=name, n_features=len(cols),
                             auroc=auc, lo=lo, hi=hi))
            print(f"  {gname:8s} {name:24s} AUROC = {auc:.3f} [{lo:.3f}, {hi:.3f}]")
    out = pd.DataFrame(rows)
    out.to_csv(f"{RES}/supervised.csv", index=False)
    SUMMARY["supervised"] = out.to_dict("records")
    return out


def section_D2(df):
    """Is whatever signal the z-score carries just set-size in disguise?

    The degree-preserving null fixes |T| and |S| exactly, but z is a standardised
    effect size: sigma_rand shrinks as the sets grow, so |z| inflates with set
    size. If the z-score's AUROC survives regressing out log|T| and log|S| it is
    carrying topology; if it collapses to chance it was carrying set size.
    """
    print("\n=== D2. is the z-score just set size? ===")
    d = df.dropna(subset=[PRIMARY]).copy()
    corr = {c: {k: float(np.corrcoef(d[c], np.log10(d[k]))[0, 1])
                for k in ("n_targets", "n_genes")}
            for c in ["z_closest", "z_shortest", "z_kernel", "d_closest"]
            if c in d.columns}
    for c, v in corr.items():
        print(f"  corr({c}, log10|T|) = {v['n_targets']:+.3f}   "
              f"corr({c}, log10|S|) = {v['n_genes']:+.3f}")
    X = np.column_stack([np.log10(d.n_targets), np.log10(d.n_genes),
                         np.ones(len(d))])
    rows = []
    for c in [PRIMARY, "z_shortest", "z_kernel", "d_closest", BASELINE]:
        if c not in d.columns:
            continue
        yv = d[c].values.astype(float)
        beta, *_ = np.linalg.lstsq(X, yv, rcond=None)
        resid = yv - X @ beta
        a0, _, _ = N.auroc_auprc(d.label.values, N.SCORE_DIRECTION[c] * yv)
        a1, _, _ = N.auroc_auprc(d.label.values, N.SCORE_DIRECTION[c] * resid)
        lo, hi = N.bootstrap_auc(d.label.values, N.SCORE_DIRECTION[c] * resid,
                                 d.drug_key.values, n_boot=1000)
        rows.append(dict(score=c, auroc_raw=a0, auroc_residualised=a1,
                         lo=lo, hi=hi, delta=a1 - a0,
                         corr_log_nT=corr.get(c, {}).get("n_targets"),
                         corr_log_nS=corr.get(c, {}).get("n_genes")))
        print(f"  {c:12s} AUROC {a0:.3f} -> {a1:.3f} after removing "
              f"log|T|, log|S|  [{lo:.3f}, {hi:.3f}]")
    r = pd.DataFrame(rows)
    r.to_csv(f"{RES}/residualised_auc.csv", index=False)
    SUMMARY["residualised"] = r.to_dict("records")
    SUMMARY["setsize_correlations"] = corr
    return r


# ================================================ E. stratifications
def section_E(df):
    print("\n=== E. stratification ===")
    df = df.copy()
    df["failure_reason"] = df.DetailedStatus.map(N.classify_failure)
    df.loc[df.label == 1, "failure_reason"] = "n/a (approved)"
    counts = df[df.label == 0].failure_reason.value_counts()
    print("  failure reasons:", counts.to_dict())

    rows = []
    pos = df[df.label == 1]
    for reason in ["efficacy", "safety", "administrative", "unstated", "other"]:
        sub = pd.concat([pos, df[(df.label == 0) & (df.failure_reason == reason)]])
        if (sub.label == 0).sum() < 30:
            continue
        for col in [PRIMARY, BASELINE, "d_closest"]:
            auc, ap, n = N.auroc_auprc(sub.label.values, N.oriented(sub, col))
            lo, hi = N.bootstrap_auc(sub.label.values, N.oriented(sub, col),
                                     sub.drug_key.values, n_boot=1000)
            rows.append(dict(stratum=f"failure: {reason}", score=col,
                             n_neg=int((sub.label == 0).sum()), n=n,
                             auroc=auc, lo=lo, hi=hi))
    # disease type
    df["is_cancer"] = df.sem_type.eq("Neoplastic Process")
    for nm, sub in (("cancer indications", df[df.is_cancer]),
                    ("non-cancer indications", df[~df.is_cancer])):
        for col in [PRIMARY, BASELINE, "d_closest"]:
            auc, ap, n = N.auroc_auprc(sub.label.values, N.oriented(sub, col))
            lo, hi = N.bootstrap_auc(sub.label.values, N.oriented(sub, col),
                                     sub.drug_key.values, n_boot=1000)
            rows.append(dict(stratum=nm, score=col,
                             n_neg=int((sub.label == 0).sum()), n=n,
                             auroc=auc, lo=lo, hi=hi))
    # module size terciles -- proximity is known to depend on |S|
    q = df.n_genes.quantile([1 / 3, 2 / 3]).values
    df["size_bin"] = np.where(df.n_genes <= q[0], "small module",
                     np.where(df.n_genes <= q[1], "medium module", "large module"))
    for nm, sub in df.groupby("size_bin"):
        for col in [PRIMARY, BASELINE]:
            auc, ap, n = N.auroc_auprc(sub.label.values, N.oriented(sub, col))
            lo, hi = N.bootstrap_auc(sub.label.values, N.oriented(sub, col),
                                     sub.drug_key.values, n_boot=1000)
            rows.append(dict(stratum=nm, score=col,
                             n_neg=int((sub.label == 0).sum()), n=n,
                             auroc=auc, lo=lo, hi=hi))
    strat = pd.DataFrame(rows)
    strat.to_csv(f"{RES}/stratified_auc.csv", index=False)
    print(strat.round(3).to_string(index=False))
    SUMMARY["stratified"] = strat.to_dict("records")

    # effect size of the z-score difference between approved and failed
    import scipy.stats as st
    a = df.loc[df.label == 1, PRIMARY].dropna()
    b = df.loc[df.label == 0, PRIMARY].dropna()
    sp = np.sqrt(((len(a) - 1) * a.var() + (len(b) - 1) * b.var()) /
                 (len(a) + len(b) - 2))
    d_cohen = (a.mean() - b.mean()) / sp
    u = st.mannwhitneyu(a, b, alternative="two-sided")
    SUMMARY["effect_size"] = dict(
        mean_z_approved=float(a.mean()), mean_z_failed=float(b.mean()),
        cohens_d=float(d_cohen), mannwhitney_U=float(u.statistic),
        mannwhitney_p=float(u.pvalue))
    print(f"  mean {PRIMARY}: approved {a.mean():.3f} vs failed {b.mean():.3f}, "
          f"Cohen's d = {d_cohen:.3f}, Mann-Whitney P = {u.pvalue:.3g}")
    return df, strat


def stratified_auc(df, by, col, min_each=1):
    """Conditional (within-stratum) AUROC, pooled as a weighted mean of the
    per-stratum Mann-Whitney statistics with weights n_pos * n_neg.

    Comparing only pairs that share a disease removes every disease-level
    confound (module size, how well studied the disease is, how many failed
    oncology trials it accumulated) and asks the question a repurposing
    programme actually faces: among candidates for THIS disease, does proximity
    rank the ones that succeed above the ones that fail?
    """
    num = den = 0.0
    per, s = [], N.oriented(df, col)
    ok = np.isfinite(s)
    d = df[ok].copy()
    d["_s"] = s[ok]
    for g, sub in d.groupby(by):
        pos, neg = sub._s[sub.label == 1].values, sub._s[sub.label == 0].values
        if len(pos) < min_each or len(neg) < min_each:
            continue
        cmp = (pos[:, None] > neg[None, :]).sum() + \
              0.5 * (pos[:, None] == neg[None, :]).sum()
        w = len(pos) * len(neg)
        num += cmp
        den += w
        per.append(dict(stratum=g, auc=cmp / w, n_pos=len(pos), n_neg=len(neg)))
    return (num / den if den else np.nan), den, pd.DataFrame(per)


def section_E2(df):
    print("\n=== E2. conditional (within-disease / within-drug) AUROC ===")
    rows = []
    for by, nm in (("ind_id", "within disease"), ("drug_key", "within drug")):
        for col in [PRIMARY, "z_shortest", "z_kernel", "d_closest", BASELINE,
                    "n_targets", "n_genes"]:
            if col not in df.columns:
                continue
            auc, w, per = stratified_auc(df, by, col)
            # bootstrap over strata for a CI
            rng = np.random.default_rng(N.SEED)
            if len(per) > 5:
                vals = []
                wts = (per.n_pos * per.n_neg).values
                aucs = per.auc.values
                for _ in range(2000):
                    i = rng.integers(0, len(per), len(per))
                    vals.append(float(np.average(aucs[i], weights=wts[i])))
                lo, hi = np.percentile(vals, [2.5, 97.5])
            else:
                lo = hi = np.nan
            rows.append(dict(conditioning=nm, score=col, auroc=auc,
                             n_strata=len(per), n_comparisons=int(w),
                             lo=float(lo), hi=float(hi)))
            print(f"  {nm:15s} {col:12s} AUROC = {auc:.3f} "
                  f"[{lo:.3f}, {hi:.3f}]  ({len(per)} strata, "
                  f"{int(w)} pos-neg comparisons)")
    out = pd.DataFrame(rows)
    out.to_csv(f"{RES}/conditional_auc.csv", index=False)
    SUMMARY["conditional_auc"] = out.to_dict("records")
    return out


# ==================================================== F. negative-set ablation
def section_F():
    print("\n=== F. negative-set ablation (D4) ===")
    rows = []
    for negset in ["hard", "random", "random_matched"]:
        f = f"{RES}/scores_negatives/{HEAD_NET}__{HEAD_GENE}__{negset}.tsv.gz"
        if not os.path.exists(f):
            print(f"  missing {f}")
            continue
        d = pd.read_csv(f, sep="\t")
        for col in [PRIMARY, "z_shortest", "z_kernel", "d_closest", BASELINE,
                    "n_targets", "n_genes"]:
            auc, ap, n = N.auroc_auprc(d.label.values, N.oriented(d, col))
            lo, hi = N.bootstrap_auc(d.label.values, N.oriented(d, col),
                                     d.drug_key.values, n_boot=1000)
            rows.append(dict(negative_set=negset, score=col, auroc=auc, auprc=ap,
                             lo=lo, hi=hi, n=n,
                             n_pos=int(d.label.sum()),
                             n_neg=int((1 - d.label).sum())))
    if not rows:
        return pd.DataFrame()
    ab = pd.DataFrame(rows)
    ab.to_csv(f"{RES}/negative_set_ablation.csv", index=False)
    print(ab.round(3).to_string(index=False))
    SUMMARY["negative_ablation"] = ab.to_dict("records")
    return ab


# ==================================================== G. network sweep + removal
def section_G():
    print("\n=== G. network sweep (D2a) ===")
    rows = []
    for f in sorted(glob.glob(f"{RES}/scores/*__hetionet.tsv.gz")):
        d = pd.read_csv(f, sep="\t")
        net, gs = d.network.iloc[0], d.gene_source.iloc[0]
        for col in [PRIMARY, "z_shortest", "z_kernel", "z_centre", "d_closest",
                    BASELINE, "jaccard", "n_targets", "n_genes"]:
            auc, ap, n = N.auroc_auprc(d.label.values, N.oriented(d, col))
            lo, hi = N.bootstrap_auc(d.label.values, N.oriented(d, col),
                                     d.drug_key.values, n_boot=1000)
            rows.append(dict(network=net, gene_source=gs, score=col, auroc=auc,
                             auprc=ap, lo=lo, hi=hi, n=n, n_pairs=len(d),
                             pos_rate=float(d.label.mean())))
        # per-network DeLong for the C2 comparison
        r1, r2, diff, p, se = N.delong_test(d.label.values, N.oriented(d, PRIMARY),
                                            N.oriented(d, BASELINE))
        rows.append(dict(network=net, gene_source=gs, score="__delong__",
                         auroc=r1, auprc=r2, lo=diff, hi=p, n=len(d),
                         n_pairs=len(d), pos_rate=float(d.label.mean())))
    sweep = pd.DataFrame(rows)
    sweep.to_csv(f"{RES}/network_sweep.csv", index=False)
    piv = sweep[sweep.score.isin([PRIMARY, BASELINE, "d_closest"])].pivot_table(
        index=["gene_source", "network"], columns="score", values="auroc")
    print(piv.round(3).to_string())
    SUMMARY["network_sweep"] = sweep[sweep.score != "__delong__"].to_dict("records")
    SUMMARY["network_delong"] = sweep[sweep.score == "__delong__"][
        ["network", "gene_source", "auroc", "lo", "hi"]].rename(
        columns={"auroc": "auroc_proximity", "lo": "diff", "hi": "p"}
    ).to_dict("records")

    er = None
    if os.path.exists(f"{RES}/edge_removal.tsv.gz"):
        print("\n=== G2. edge-removal curves (D2b) ===")
        er = pd.read_csv(f"{RES}/edge_removal.tsv.gz", sep="\t")
        g = er.groupby(["network", "score", "frac_removed"]).agg(
            auroc_mean=("auroc", "mean"), auroc_sd=("auroc", "std"),
            edges=("edges", "mean"), pairs=("pairs_retained", "mean"),
            reps=("auroc", "size")).reset_index()
        g.to_csv(f"{RES}/edge_removal_summary.csv", index=False)
        print(g[g.score.isin([PRIMARY, BASELINE])].round(3).to_string(index=False))
        SUMMARY["edge_removal"] = g.to_dict("records")
    return sweep, er


# ==================================================== H. tissue specificity
def section_H():
    f = f"{RES}/tissue_scores.tsv.gz"
    if not os.path.exists(f):
        print("\n=== H. tissue arm: not available ===")
        return None
    print("\n=== H. tissue specificity (D3) ===")
    d = pd.read_csv(f, sep="\t")
    # restrict to pairs present in ALL THREE arms so the comparison is paired
    key = d.drug_key.astype(str) + "|" + d.ind_id.astype(str)
    d["key"] = key
    arms = d.arm.unique()
    common = set.intersection(*[set(d.key[d.arm == a]) for a in arms])
    d = d[d.key.isin(common)]
    print(f"  {len(common)} pairs scored in all {len(arms)} arms")
    rows = []
    for arm, sub in d.groupby("arm"):
        for col in [PRIMARY, "z_shortest", "d_closest", BASELINE]:
            auc, ap, n = N.auroc_auprc(sub.label.values, N.oriented(sub, col))
            lo, hi = N.bootstrap_auc(sub.label.values, N.oriented(sub, col),
                                     sub.drug_key.values, n_boot=1000)
            rows.append(dict(arm=arm, score=col, auroc=auc, auprc=ap, lo=lo, hi=hi,
                             n=n, mean_nodes=float(sub.tissue_nodes.mean())))
    tis = pd.DataFrame(rows)
    # paired DeLong: tissue vs global, and tissue vs shuffled, on identical pairs
    piv = d.pivot_table(index="key", columns="arm", values=PRIMARY)
    lab = d.drop_duplicates("key").set_index("key").label
    piv = piv.dropna()
    yv = lab.loc[piv.index].values
    comps = []
    for a, b in [("tissue", "global"), ("tissue", "shuffled"),
                 ("shuffled", "global")]:
        if a in piv.columns and b in piv.columns:
            r1, r2, diff, p, se = N.delong_test(yv, -piv[a].values, -piv[b].values)
            comps.append(dict(a=a, b=b, auc_a=r1, auc_b=r2, diff=diff, p=p, n=len(yv)))
    tis.to_csv(f"{RES}/tissue_auc.csv", index=False)
    pd.DataFrame(comps).to_csv(f"{RES}/tissue_delong.csv", index=False)
    print(tis.round(3).to_string(index=False))
    print(pd.DataFrame(comps).round(4).to_string(index=False))
    SUMMARY["tissue"] = tis.to_dict("records")
    SUMMARY["tissue_delong"] = comps
    return tis


# ==================================================== I. figures
def fig_forest(tbl, df):
    order = ([c for c in N.PROXIMITY_SCORES if c in set(tbl.score)]
             + N.BASELINE_SCORES)
    t = tbl.set_index("score").reindex(order).dropna(subset=["auroc"])
    fig, ax = plt.subplots(figsize=(7.4, 0.34 * len(t) + 1.5))
    ypos = np.arange(len(t))[::-1]
    cols = [C_PROX if s in N.PROXIMITY_SCORES else C_OVER if s == BASELINE
            else C_OTHER for s in t.index]
    for y, (_, r), c in zip(ypos, t.iterrows(), cols):
        ax.plot([r.lo_drug, r.hi_drug], [y, y], color=c, lw=2, solid_capstyle="round")
        ax.plot([r.auroc], [y], "o", color=c, ms=7, mec="#fcfcfb", mew=1.5, zorder=3)
        ax.text(max(r.hi_drug, r.auroc) + 0.006, y, f"{r.auroc:.3f}", va="center",
                fontsize=8, color=C_INK2)
    ax.axvline(0.5, color=C_INK2, lw=1, ls="--")
    ax.axvline(0.65, color=C_ALT, lw=1.2, ls=":")
    ax.text(0.65, len(t) - 0.2, " hypothesis bound 0.65", color=C_ALT, fontsize=8,
            va="top")
    ax.text(0.5, len(t) - 0.2, " chance", color=C_INK2, fontsize=8, va="top")
    ax.set_yticks(ypos)
    ax.set_yticklabels([N.PRETTY.get(s, s) for s in t.index], fontsize=8)
    ax.set_xlabel("AUROC (bar = 95% CI, bootstrap resampled by drug)")
    ax.set_title(f"Predicting clinical repurposing success on repoDB\n"
                 f"{HEAD_NET} interactome x DisGeNET-{HEAD_GENE}, "
                 f"{len(df)} pairs ({df.label.sum()} approved / "
                 f"{(1 - df.label).sum()} clinically failed)", fontsize=10, loc="left")
    ax.set_xlim(0.3, 0.75)
    ax.grid(axis="y", visible=False)
    fig.savefig(f"{FIG}/fig1_headline_auroc.png")
    plt.close(fig)


def fig_dist_roc(df):
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.5))
    ax = axes[0]
    bins = np.linspace(-12, 6, 46)
    for lab, name, c in ((1, "approved", C_PROX), (0, "clinically failed", C_OVER)):
        v = df.loc[df.label == lab, PRIMARY].dropna()
        ax.hist(v, bins=bins, density=True, histtype="step", lw=2, color=c,
                label=f"{name} (n={len(v)})")
        ax.axvline(v.mean(), color=c, lw=1, ls=":")
    ax.set_xlabel("proximity z-score (closest)")
    ax.set_ylabel("density")
    ax.set_title("a  Proximity z by clinical outcome", loc="left", fontsize=10)
    ax.legend(fontsize=8)

    ax = axes[1]
    for col, c, nm in ((PRIMARY, C_PROX, "proximity z (closest)"),
                       (BASELINE, C_OVER, "target/disease-gene overlap"),
                       ("d_closest", C_OTHER, "raw closest distance")):
        s = N.oriented(df, col)
        ok = np.isfinite(s)
        fpr, tpr, _ = roc_curve(df.label.values[ok], s[ok])
        ax.plot(fpr, tpr, color=c, lw=2,
                label=f"{nm} (AUC {roc_auc_score(df.label.values[ok], s[ok]):.3f})")
    ax.plot([0, 1], [0, 1], color=C_INK2, lw=1, ls="--")
    ax.set_xlabel("false positive rate"); ax.set_ylabel("true positive rate")
    ax.set_title("b  ROC", loc="left", fontsize=10)
    ax.legend(fontsize=8, loc="lower right")

    ax = axes[2]
    # is z proximal-vs-distant informative at all? approval rate by z decile
    d = df.dropna(subset=[PRIMARY]).copy()
    d["dec"] = pd.qcut(d[PRIMARY], 10, labels=False, duplicates="drop")
    g = d.groupby("dec").agg(rate=("label", "mean"), n=("label", "size"),
                             zmid=(PRIMARY, "median"))
    ax.bar(g.index, g.rate, color=C_PROX, width=0.72)
    ax.axhline(df.label.mean(), color=C_OVER, lw=1.5, ls="--",
               label=f"base rate {df.label.mean():.2f}")
    for i, r in g.iterrows():
        ax.text(i, r.rate + 0.012, f"{r.rate:.2f}", ha="center", fontsize=7,
                color=C_INK2)
    ax.set_xlabel("proximity z decile (0 = most proximal)")
    ax.set_ylabel("fraction approved")
    ax.set_ylim(0, 1)
    ax.set_title("c  Approval rate by proximity decile", loc="left", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{FIG}/fig2_distributions_roc.png")
    plt.close(fig)


def fig_negatives(ab):
    if ab is None or not len(ab):
        return
    sub = ab[ab.score.isin([PRIMARY, BASELINE, "d_closest"])]
    negs = ["hard", "random_matched", "random"]
    names = {"hard": "repoDB clinical\nfailures (hard)",
             "random_matched": "random unobserved\n(disease-matched)",
             "random": "random unobserved\n(Guney-style)"}
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    w, cols = 0.26, {PRIMARY: C_PROX, BASELINE: C_OVER, "d_closest": C_OTHER}
    for i, col in enumerate([PRIMARY, BASELINE, "d_closest"]):
        vals = [sub[(sub.negative_set == n) & (sub.score == col)] for n in negs]
        x = np.arange(len(negs)) + (i - 1) * w
        h = [float(v.auroc.iloc[0]) if len(v) else np.nan for v in vals]
        lo = [float(v.lo.iloc[0]) if len(v) else np.nan for v in vals]
        hi = [float(v.hi.iloc[0]) if len(v) else np.nan for v in vals]
        ax.bar(x, h, w * 0.92, color=cols[col], label=N.PRETTY[col])
        ax.errorbar(x, h, yerr=[np.array(h) - lo, np.array(hi) - np.array(h)],
                    fmt="none", ecolor="#fcfcfb", elinewidth=1.6, capsize=0)
        ax.errorbar(x, h, yerr=[np.array(h) - lo, np.array(hi) - np.array(h)],
                    fmt="none", ecolor=C_INK2, elinewidth=0.9, capsize=2)
        for xi, hi_ in zip(x, h):
            ax.text(xi, hi_ + 0.015, f"{hi_:.2f}", ha="center", fontsize=7.5,
                    color=C_INK2)
    ax.axhline(0.5, color=C_INK2, lw=1, ls="--")
    ax.axhline(0.66, color=C_ALT, lw=1.2, ls=":")
    ax.text(2.42, 0.665, "Guney et al.\nreported 0.66", color=C_ALT, fontsize=7.5,
            va="bottom", ha="right")
    ax.set_xticks(np.arange(len(negs)))
    ax.set_xticklabels([names[n] for n in negs], fontsize=8.5)
    ax.set_ylabel("AUROC (95% CI, by drug)")
    ax.set_ylim(0.3, 0.85)
    ax.set_title("Same positives, same network, same code — only the negatives change",
                 loc="left", fontsize=10)
    ax.legend(fontsize=8, ncol=3, loc="upper left")
    ax.grid(axis="x", visible=False)
    fig.savefig(f"{FIG}/fig3_negative_set_ablation.png")
    plt.close(fig)


def fig_sweep(sweep):
    s = sweep[(sweep.score.isin([PRIMARY, BASELINE])) & (sweep.gene_source == "all")]
    nets = ["biogrid", "string700", "hetionet_gig", "biogrid_lt", "lit_bm",
            "hi_union", "huri"]
    nets = [n for n in nets if n in set(s.network)]
    fig, ax = plt.subplots(figsize=(8.2, 4.0))
    w = 0.36
    for i, (col, c) in enumerate([(PRIMARY, C_PROX), (BASELINE, C_OVER)]):
        sub = s[s.score == col].set_index("network").reindex(nets)
        x = np.arange(len(nets)) + (i - 0.5) * w
        ax.bar(x, sub.auroc, w * 0.92, color=c, label=N.PRETTY[col])
        ax.errorbar(x, sub.auroc,
                    yerr=[sub.auroc - sub.lo, sub.hi - sub.auroc], fmt="none",
                    ecolor=C_INK2, elinewidth=0.9, capsize=2)
        for xi, v in zip(x, sub.auroc):
            ax.text(xi, v + 0.012, f"{v:.2f}", ha="center", fontsize=7.5,
                    color=C_INK2)
    ax.axhline(0.5, color=C_INK2, lw=1, ls="--")
    npairs = sweep[(sweep.gene_source == "all") &
                   (sweep.score == PRIMARY)].set_index("network").reindex(nets)
    ax.set_xticks(np.arange(len(nets)))
    ax.set_xticklabels([f"{n}\n{int(v)} pairs" for n, v in
                        zip(nets, npairs.n_pairs)], fontsize=8)
    ax.set_ylabel("AUROC (95% CI, by drug)")
    ax.set_ylim(0.3, 0.72)
    ax.set_title("Across seven interactomes (DisGeNET-all disease genes)\n"
                 "left group = literature-curated, right group = systematic "
                 "(study-bias-free)", loc="left", fontsize=10)
    ax.axvline(2.5, color=C_GRID, lw=8, zorder=0)
    ax.legend(fontsize=8, ncol=2, loc="upper right")
    ax.grid(axis="x", visible=False)
    fig.savefig(f"{FIG}/fig4_network_sweep.png")
    plt.close(fig)


def fig_edge_removal(er):
    if er is None or not len(er):
        return
    g = er.groupby(["network", "score", "frac_removed"]).auroc.agg(["mean", "std"]) \
          .reset_index()
    nets = sorted(g.network.unique())
    fig, axes = plt.subplots(1, len(nets), figsize=(5.0 * len(nets), 3.6),
                             squeeze=False)
    for ax, net in zip(axes[0], nets):
        for col, c in [(PRIMARY, C_PROX), (BASELINE, C_OVER), ("d_closest", C_OTHER)]:
            sub = g[(g.network == net) & (g.score == col)].sort_values("frac_removed")
            ax.plot(sub.frac_removed, sub["mean"], "-o", color=c, lw=2, ms=5,
                    mec="#fcfcfb", mew=1.2, label=N.PRETTY[col])
            ax.fill_between(sub.frac_removed, sub["mean"] - sub["std"].fillna(0),
                            sub["mean"] + sub["std"].fillna(0), color=c, alpha=0.15,
                            lw=0)
        ax.axhline(0.5, color=C_INK2, lw=1, ls="--")
        ax.set_xlabel("fraction of edges randomly deleted")
        ax.set_ylabel("AUROC")
        ax.set_title(f"{net}", loc="left", fontsize=10)
        ax.set_ylim(0.35, 0.65)
    axes[0][0].legend(fontsize=8)
    fig.suptitle("Sensitivity to interactome incompleteness (3 deletion seeds; "
                 "band = ±1 SD)", fontsize=10, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(f"{FIG}/fig5_edge_removal.png")
    plt.close(fig)


def fig_strat_tissue(strat, tis):
    n = 2 if tis is not None else 1
    fig, axes = plt.subplots(1, n, figsize=(6.0 * n, 4.0), squeeze=False)
    ax = axes[0][0]
    s = strat[strat.score.isin([PRIMARY, BASELINE])]
    strata = [x for x in ["failure: efficacy", "failure: safety",
                          "failure: administrative", "failure: unstated",
                          "cancer indications", "non-cancer indications",
                          "small module", "medium module", "large module"]
              if x in set(s.stratum)]
    w = 0.36
    for i, (col, c) in enumerate([(PRIMARY, C_PROX), (BASELINE, C_OVER)]):
        sub = s[s.score == col].set_index("stratum").reindex(strata)
        y = np.arange(len(strata))[::-1] + (0.5 - i) * w
        ax.barh(y, sub.auroc, w * 0.92, color=c, label=N.PRETTY[col])
        ax.errorbar(sub.auroc, y, xerr=[sub.auroc - sub.lo, sub.hi - sub.auroc],
                    fmt="none", ecolor=C_INK2, elinewidth=0.9, capsize=2)
        for yi, v in zip(y, sub.auroc):
            if np.isfinite(v):
                ax.text(v + 0.008, yi, f"{v:.2f}", va="center", fontsize=7.5,
                        color=C_INK2)
    ax.axvline(0.5, color=C_INK2, lw=1, ls="--")
    ax.set_yticks(np.arange(len(strata))[::-1])
    ax.set_yticklabels(strata, fontsize=8)
    ax.set_xlim(0.3, 0.75)
    ax.set_xlabel("AUROC (95% CI, by drug)")
    ax.set_title("a  Stratified performance", loc="left", fontsize=10)
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(axis="y", visible=False)

    if tis is not None:
        ax = axes[0][1]
        arms = [a for a in ["global", "tissue", "shuffled"] if a in set(tis.arm)]
        names = {"global": "tissue-naive\n(full interactome)",
                 "tissue": "tissue-restricted\n(disease's own tissue)",
                 "shuffled": "tissue-restricted\n(random tissue, control)"}
        for i, (col, c) in enumerate([(PRIMARY, C_PROX), (BASELINE, C_OVER)]):
            sub = tis[tis.score == col].set_index("arm").reindex(arms)
            x = np.arange(len(arms)) + (i - 0.5) * w
            ax.bar(x, sub.auroc, w * 0.92, color=c, label=N.PRETTY[col])
            ax.errorbar(x, sub.auroc, yerr=[sub.auroc - sub.lo, sub.hi - sub.auroc],
                        fmt="none", ecolor=C_INK2, elinewidth=0.9, capsize=2)
            for xi, v in zip(x, sub.auroc):
                if np.isfinite(v):
                    ax.text(xi, v + 0.01, f"{v:.2f}", ha="center", fontsize=7.5,
                            color=C_INK2)
        ax.axhline(0.5, color=C_INK2, lw=1, ls="--")
        ax.set_xticks(np.arange(len(arms)))
        ax.set_xticklabels([names[a] for a in arms], fontsize=8)
        ax.set_ylabel("AUROC (95% CI, by drug)")
        ax.set_ylim(0.3, 0.7)
        ax.set_title("b  Tissue-specific interactomes (paired pairs)", loc="left",
                     fontsize=10)
        ax.legend(fontsize=8)
        ax.grid(axis="x", visible=False)
    fig.tight_layout()
    fig.savefig(f"{FIG}/fig6_stratified_tissue.png")
    plt.close(fig)


# ==================================================================== main
def main():
    N.set_seed()
    section_A()
    df, tbl, dl = section_BC()
    section_D(df)
    section_D2(df)
    dfs, strat = section_E(df)
    section_E2(df)
    ab = section_F()
    sweep, er = section_G()
    tis = section_H()

    print("\n=== I. figures ===")
    fig_forest(tbl, df); print("  fig1")
    fig_dist_roc(df); print("  fig2")
    fig_negatives(ab); print("  fig3")
    fig_sweep(sweep); print("  fig4")
    fig_edge_removal(er); print("  fig5")
    fig_strat_tissue(strat, tis); print("  fig6")

    if os.path.exists(f"{RES}/validation.json"):
        SUMMARY["validation"] = json.load(open(f"{RES}/validation.json"))
    json.dump(SUMMARY, open(f"{RES}/summary.json", "w"), indent=2, default=float)
    print(f"\nwrote {RES}/summary.json")


if __name__ == "__main__":
    main()
