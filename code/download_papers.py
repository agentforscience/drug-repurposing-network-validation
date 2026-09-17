#!/usr/bin/env python
"""Download PDFs for top-ranked papers.
Resolution order: Europe PMC OA render -> S2 openAccessPdf -> OpenAlex OA location -> arXiv/bioRxiv.
Only free/open-access routes are used."""
import json, os, re, sys, time
import httpx

os.makedirs('papers', exist_ok=True)
UA = {'User-Agent': 'Mozilla/5.0 (compatible; academic-literature-review)'}
N = int(sys.argv[1]) if len(sys.argv) > 1 else 60

def slug(t, n=70):
    return re.sub(r'[^A-Za-z0-9]+', '_', (t or 'untitled')).strip('_')[:n]

def try_get(url, path, client):
    try:
        r = client.get(url, headers=UA, follow_redirects=True, timeout=180)
        if r.status_code != 200: return False, f'http{r.status_code}'
        c = r.content
        if not c.startswith(b'%PDF'): return False, 'not-pdf'
        if len(c) < 20000: return False, f'small({len(c)})'
        open(path, 'wb').write(c)
        return True, f'{len(c)//1024}KB'
    except Exception as e:
        return False, type(e).__name__

def openalex_pdf(doi, client):
    if not doi: return []
    try:
        r = client.get(f'https://api.openalex.org/works/https://doi.org/{doi}',
                       headers=UA, timeout=60, follow_redirects=True)
        if r.status_code != 200: return []
        d = r.json()
        urls = []
        for loc in ([d.get('best_oa_location')] + (d.get('locations') or [])):
            if not loc: continue
            u = loc.get('pdf_url')
            if u: urls.append(u)
        pmcid = (d.get('ids') or {}).get('pmid')
        return list(dict.fromkeys(urls))
    except Exception:
        return []

papers = [json.loads(l) for l in open('paper_search_results/ranked.jsonl')][:N]
manifest = []
with httpx.Client() as client:
    for i, p in enumerate(papers):
        title = p.get('title') or ''
        fn = f"papers/{i+1:03d}_{slug(title)}.pdf"
        if os.path.exists(fn) and os.path.getsize(fn) > 20000:
            print(f'{i+1:3d} SKIP', flush=True); manifest.append((p, fn, 'cached')); continue
        doi = p.get('doi') or ''
        cands = []
        if p.get('pmcid'):
            cands.append(f"https://europepmc.org/articles/{p['pmcid']}?pdf=render")
        if p.get('oa_pdf'): cands.append(p['oa_pdf'])
        if p.get('arxiv'): cands.append(f"https://arxiv.org/pdf/{p['arxiv']}.pdf")
        cands += openalex_pdf(doi, client)
        if doi.startswith('10.1101'):
            cands.append(f"https://www.biorxiv.org/content/{doi}v1.full.pdf")
            cands.append(f"https://www.medrxiv.org/content/{doi}v1.full.pdf")
        ok, why = False, 'no-candidates'
        for u in dict.fromkeys(cands):
            ok, why = try_get(u, fn, client)
            if ok: break
            time.sleep(0.4)
        print(f'{i+1:3d} {"OK  " if ok else "FAIL"} {why:20s} {title[:75]}', flush=True)
        manifest.append((p, fn if ok else None, why))
        time.sleep(0.4)

json.dump([{'rank': i+1, 'title': p.get('title'), 'file': fn, 'status': why,
            'doi': p.get('doi'), 'pmcid': p.get('pmcid'), 'pmid': p.get('pmid'),
            'year': p.get('year'), 'authors': p.get('authors'), 'journal': p.get('journal'),
            'citations': p.get('citations'), 'abstract': p.get('abstract'), 'score': p.get('_score')}
           for i, (p, fn, why) in enumerate(manifest)],
          open('paper_search_results/download_manifest.json', 'w'), indent=1)
print(f'\nAvailable: {sum(1 for _, fn, _ in manifest if fn)}/{len(manifest)}')
