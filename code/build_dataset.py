#!/usr/bin/env python
"""
Assemble the benchmark for: does network proximity predict CLINICAL repurposing success?

Outputs (datasets/processed/):
  ppi_<name>.tsv          edge lists (Entrez gene id pairs), largest connected component
  drug_targets.tsv        drug_id (DrugBank) -> entrez targets
  disease_genes_<src>.tsv disease cui -> entrez genes
  pairs.tsv               repoDB drug-indication pairs w/ label + coverage flags
"""
import pandas as pd, numpy as np, gzip, os, re, json, sys, collections
import networkx as nx

OUT = 'datasets/processed'; os.makedirs(OUT, exist_ok=True)
def log(*a): print(*a, flush=True)

# ---------------------------------------------------------------- ID mappings
log('== gene id mappings ==')
gi = pd.read_csv('datasets/mappings/Homo_sapiens.gene_info.gz', sep='\t', low_memory=False,
                 usecols=['#tax_id','GeneID','Symbol','Synonyms','type_of_gene'])
gi = gi[gi['#tax_id'] == 9606]
sym2entrez = dict(zip(gi.Symbol, gi.GeneID))
for _, r in gi.iterrows():
    if r.Synonyms and r.Synonyms != '-':
        for s in str(r.Synonyms).split('|'):
            sym2entrez.setdefault(s, r.GeneID)
valid_entrez = set(gi.GeneID)
log(f'  symbols->entrez: {len(sym2entrez)}, human genes: {len(valid_entrez)}')

# Ensembl gene -> Entrez (for HuRI)
ens2entrez = {}
with gzip.open('datasets/mappings/gene2ensembl.gz', 'rt') as f:
    hdr = f.readline().rstrip('\n').lstrip('#').split('\t')
    ix = {c: i for i, c in enumerate(hdr)}
    for line in f:
        p = line.rstrip('\n').split('\t')
        if p[ix['tax_id']] != '9606': continue
        ens2entrez.setdefault(p[ix['Ensembl_gene_identifier']], int(p[ix['GeneID']]))
log(f'  ensembl->entrez: {len(ens2entrez)}')

# STRING protein id -> gene symbol
str2sym = {}
with gzip.open('datasets/ppi/9606.protein.info.v12.0.txt.gz', 'rt') as f:
    next(f)
    for line in f:
        p = line.rstrip('\n').split('\t')
        str2sym[p[0]] = p[1]
log(f'  string->symbol: {len(str2sym)}')

# ---------------------------------------------------------------- PPI networks
def save_net(name, edges):
    """edges: iterable of (entrez,entrez). Saves LCC edge list."""
    G = nx.Graph()
    G.add_edges_from(edges)
    G.remove_edges_from(nx.selfloop_edges(G))
    if G.number_of_nodes() == 0:
        log(f'  {name}: EMPTY'); return None
    lcc = max(nx.connected_components(G), key=len)
    G = G.subgraph(lcc).copy()
    df = pd.DataFrame(sorted(G.edges()), columns=['gene_a', 'gene_b'])
    df.to_csv(f'{OUT}/ppi_{name}.tsv', sep='\t', index=False)
    log(f'  {name}: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges (LCC)')
    return G

log('== PPI networks ==')
nets = {}

# 1) BioGRID human physical (literature-curated; carries study bias)
bg = pd.read_csv('datasets/ppi/biogrid_human_physical.tsv.gz', sep='\t', low_memory=False)
bg = bg[['Entrez Gene Interactor A', 'Entrez Gene Interactor B', 'Throughput']].dropna()
bg = bg[bg['Entrez Gene Interactor A'].astype(str).str.isdigit() &
        bg['Entrez Gene Interactor B'].astype(str).str.isdigit()]
bg['a'] = bg['Entrez Gene Interactor A'].astype(int); bg['b'] = bg['Entrez Gene Interactor B'].astype(int)
bg = bg[bg.a.isin(valid_entrez) & bg.b.isin(valid_entrez)]
nets['biogrid'] = save_net('biogrid', zip(bg.a, bg.b))
# low-throughput only == the most literature-biased subset
lt = bg[bg.Throughput.astype(str).str.contains('Low Throughput', na=False)]
nets['biogrid_lt'] = save_net('biogrid_lt', zip(lt.a, lt.b))

# 2) HuRI: systematic Y2H, study-bias-free by construction
for nm, fn in [('huri', 'HuRI.tsv'), ('hi_union', 'HI-union.tsv'), ('lit_bm', 'Lit-BM.tsv')]:
    h = pd.read_csv(f'datasets/ppi/{fn}', sep='\t', header=None, names=['a', 'b'])
    e = [(ens2entrez[a], ens2entrez[b]) for a, b in zip(h.a, h.b)
         if a in ens2entrez and b in ens2entrez]
    nets[nm] = save_net(nm, e)

# 3) STRING high-confidence (>=700) physical+functional
rows = []
with gzip.open('datasets/ppi/9606.protein.links.v12.0.txt.gz', 'rt') as f:
    next(f)
    for line in f:
        a, b, s = line.split()
        if int(s) < 700: continue
        sa, sb = str2sym.get(a), str2sym.get(b)
        if sa in sym2entrez and sb in sym2entrez:
            rows.append((sym2entrez[sa], sym2entrez[sb]))
nets['string700'] = save_net('string700', rows)

# 4) Hetionet Gene-interacts-Gene (the PPI used by Rephetio)
he = pd.read_csv('datasets/hetionet/hetionet-v1.0-edges.sif.gz', sep='\t')
gig = he[he.metaedge == 'GiG']
e = [(int(a.split('::')[1]), int(b.split('::')[1])) for a, b in zip(gig.source, gig.target)]
nets['hetionet_gig'] = save_net('hetionet_gig', e)

# ---------------------------------------------------------------- drug targets
log('== drug targets ==')
cbg = he[he.metaedge == 'CbG']
dt = pd.DataFrame({'drug_id': cbg.source.str.replace('Compound::', '', regex=False),
                   'entrez': [int(x.split('::')[1]) for x in cbg.target]})
dt = dt.drop_duplicates()
dt.to_csv(f'{OUT}/drug_targets_hetionet.tsv', sep='\t', index=False)
log(f'  hetionet CbG: {len(dt)} edges, {dt.drug_id.nunique()} drugs, {dt.entrez.nunique()} targets')

# DrugCentral (alternative drug-target source, name-matched)
dc = pd.read_csv('datasets/drugcentral/drug.target.interaction.tsv.gz', sep='\t', low_memory=False)
dc = dc[dc.ORGANISM == 'Homo sapiens'][['DRUG_NAME', 'GENE']].dropna()
rows = []
for nm, g in zip(dc.DRUG_NAME, dc.GENE):
    for sym in str(g).split('|'):
        if sym in sym2entrez: rows.append((str(nm).lower(), sym2entrez[sym]))
dcdf = pd.DataFrame(rows, columns=['drug_name_lc', 'entrez']).drop_duplicates()
dcdf.to_csv(f'{OUT}/drug_targets_drugcentral.tsv', sep='\t', index=False)
log(f'  drugcentral: {len(dcdf)} edges, {dcdf.drug_name_lc.nunique()} drugs')

# ---------------------------------------------------------------- disease genes
log('== disease genes ==')
D = 'datasets/disgenet/dhimmel-disgenet-a636cc0/download'
for name in ['curated', 'all']:
    dg = pd.read_csv(f'{D}/{name}_gene_disease_associations.txt.gz', sep='\t')
    dg['cui'] = dg.diseaseId.str.replace('umls:', '', regex=False)
    dg = dg[dg.geneId.isin(valid_entrez)][['cui', 'geneId', 'score']].rename(columns={'geneId': 'entrez'})
    dg.drop_duplicates().to_csv(f'{OUT}/disease_genes_disgenet_{name}.tsv', sep='\t', index=False)
    log(f'  disgenet_{name}: {len(dg)} assoc, {dg.cui.nunique()} diseases')

# ---------------------------------------------------------------- repoDB pairs
log('== repoDB pairs ==')
repo = pd.read_csv('datasets/repodb/repodb_full.csv')
repo = repo[['drug_name', 'drug_id', 'ind_name', 'ind_id', 'sem_type', 'status', 'phase', 'DetailedStatus']]
repo['label'] = (repo.status == 'Approved').astype(int)
repo['drug_name_lc'] = repo.drug_name.str.lower()
repo.to_csv(f'{OUT}/pairs.tsv', sep='\t', index=False)
log(f'  pairs: {len(repo)} ({repo.label.sum()} approved / {(1-repo.label).sum()} failed)')

# ---------------------------------------------------------------- coverage report
log('== coverage (drug has targets AND disease has genes AND both in network) ==')
cov = {}
for net_name, G in nets.items():
    if G is None: continue
    nodes = set(G.nodes())
    d_ok = dt[dt.entrez.isin(nodes)].groupby('drug_id').size()
    d_ok = set(d_ok.index)
    for dgname in ['curated', 'all']:
        dg = pd.read_csv(f'{OUT}/disease_genes_disgenet_{dgname}.tsv', sep='\t')
        s_ok = set(dg[dg.entrez.isin(nodes)].groupby('cui').size().index)
        sub = repo[repo.drug_id.isin(d_ok) & repo.ind_id.isin(s_ok)]
        cov[(net_name, dgname)] = dict(pairs=len(sub), pos=int(sub.label.sum()),
                                       neg=int((1 - sub.label).sum()),
                                       drugs=sub.drug_id.nunique(), diseases=sub.ind_id.nunique(),
                                       nodes=G.number_of_nodes(), edges=G.number_of_edges())
        log(f'  {net_name:14s} x {dgname:8s}: {len(sub):5d} pairs '
            f'({int(sub.label.sum())}+/{int((1-sub.label).sum())}-) '
            f'{sub.drug_id.nunique()} drugs, {sub.ind_id.nunique()} diseases')
json.dump({f'{k[0]}|{k[1]}': v for k, v in cov.items()}, open(f'{OUT}/coverage.json', 'w'), indent=1)
log('\nDONE')
