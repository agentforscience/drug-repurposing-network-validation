#!/usr/bin/env python
"""Score & rank candidate papers for relevance to the hypothesis."""
import json, re, sys, collections

HYP_TERMS = {
    # core mechanism terms (high weight)
    'network proximity': 6, 'network-based': 4, 'network medicine': 5,
    'disease module': 5, 'drug repurposing': 4, 'drug repositioning': 4,
    'interactome': 4, 'proximity': 3, 'protein-protein interaction': 3,
    'ppi': 2, 'drug target': 3, 'drug-target': 3,
    # evaluation / critique terms (high value for this hypothesis)
    'auroc': 4, 'auc': 2, 'benchmark': 3, 'gold standard': 4,
    'evaluation': 2, 'validation': 2, 'negative': 2, 'baseline': 3,
    'bias': 4, 'confound': 5, 'incomplete': 4, 'pitfall': 5,
    'degree': 2, 'hub': 2, 'study bias': 6, 'literature bias': 5,
    'reproducib': 4, 'critical assessment': 5, 'overestimat': 5,
    # tissue
    'tissue-specific': 5, 'tissue specific': 5, 'expression': 1,
    # clinical outcome
    'clinical trial': 4, 'repodb': 6, 'approved': 2, 'failed': 3,
    'phase ii': 3, 'phase 2': 2, 'success': 2, 'attrition': 4,
    # methods
    'separation': 3, 'random walk': 2, 'diffusion': 2, 'proximity measure': 5,
    'guilt-by-association': 4, 'knowledge graph': 2,
}
NEG_TERMS = ['traditional chinese medicine', 'tcm', 'herbal', 'formula',
             'molecular docking', 'molecular dynamics simulation', 'in vitro',
             'rat', 'mice model', 'plant extract', 'decoction']

def score(p):
    title = (p.get('title') or '').lower()
    abst = (p.get('abstract') or '').lower()
    text = title + ' ' + abst
    s = 0.0
    hits = []
    for t, w in HYP_TERMS.items():
        if t in title:
            s += w * 2; hits.append(t + '(T)')
        elif t in abst:
            s += w; hits.append(t)
    # penalty for off-topic network pharmacology (TCM-style papers dominate this literature)
    for t in NEG_TERMS:
        if t in text: s -= 6
    # citation boost (log-ish)
    c = p.get('citations') or 0
    try: c = int(c)
    except: c = 0
    s += min(c, 2000) ** 0.5 / 3.0
    # recency mild boost, seminal papers still rank via citations
    try:
        y = int(p.get('year') or 0)
        if y >= 2018: s += 2
        if y >= 2022: s += 1
    except: pass
    if not abst: s -= 3
    return s, hits

papers = []
seen = set()
for line in open('paper_search_results/raw_hits.jsonl'):
    try: p = json.loads(line)
    except: continue
    t = re.sub(r'[^a-z0-9]', '', (p.get('title') or '').lower())
    if not t: continue
    if t in seen:
        continue
    seen.add(t)
    papers.append(p)

for p in papers:
    p['_score'], p['_hits'] = score(p)

papers.sort(key=lambda x: -x['_score'])
with open('paper_search_results/ranked.jsonl', 'w') as f:
    for p in papers: f.write(json.dumps(p) + '\n')

print(f'unique papers: {len(papers)}\n')
for i, p in enumerate(papers[:60]):
    print(f"{i+1:3d}. [{p['_score']:6.1f}] ({p.get('year')}) c={p.get('citations')} {p.get('source')}")
    print(f"     {(p.get('title') or '')[:150]}")
    print(f"     doi={p.get('doi')} pmcid={p.get('pmcid')} oa={p.get('isOpenAccess') or p.get('oa_pdf')}")
