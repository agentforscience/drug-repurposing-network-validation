#!/usr/bin/env python
"""Retry failed paper downloads with backoff and extra OA routes."""
import json, os, re, time
import httpx
UA = {'User-Agent': 'Mozilla/5.0 (compatible; academic-literature-review)'}
def slug(t, n=70): return re.sub(r'[^A-Za-z0-9]+','_',(t or 'x')).strip('_')[:n]
man = json.load(open('paper_search_results/download_manifest.json'))

def try_get(url, path, client, sleep=0):
    if sleep: time.sleep(sleep)
    try:
        r = client.get(url, headers=UA, follow_redirects=True, timeout=180)
        if r.status_code != 200: return False, f'http{r.status_code}'
        c = r.content
        if not c.startswith(b'%PDF'): return False, 'not-pdf'
        if len(c) < 20000: return False, f'small({len(c)})'
        open(path,'wb').write(c); return True, f'{len(c)//1024}KB'
    except Exception as e: return False, type(e).__name__

def routes(p, client):
    out = []
    pm = p.get('pmcid'); doi = p.get('doi') or ''
    if pm:
        out += [f"https://europepmc.org/articles/{pm}?pdf=render",
                f"https://europepmc.org/api/fulltextRepo?pprId={pm}&type=FILE&fileName=EMS.pdf",
                f"https://pmc.ncbi.nlm.nih.gov/articles/{pm}/pdf/{pm}.pdf"]
    if doi:
        try:
            r = client.get(f'https://api.openalex.org/works/https://doi.org/{doi}', headers=UA, timeout=60, follow_redirects=True)
            if r.status_code == 200:
                d = r.json()
                for loc in ([d.get('best_oa_location')] + (d.get('locations') or [])):
                    if loc and loc.get('pdf_url'): out.append(loc['pdf_url'])
                    if loc and loc.get('landing_page_url','').startswith('https://www.ncbi.nlm.nih.gov/pmc/articles/PMC'):
                        p2 = loc['landing_page_url'].rstrip('/').split('/')[-1]
                        out.append(f"https://europepmc.org/articles/{p2}?pdf=render")
        except Exception: pass
    if doi.startswith('10.1101'):
        out += [f"https://www.biorxiv.org/content/{doi}v1.full.pdf",
                f"https://www.medrxiv.org/content/{doi}v1.full.pdf"]
    return list(dict.fromkeys(out))

fixed = 0
with httpx.Client() as client:
    for e in man:
        if e.get('file'): continue
        fn = f"papers/{e['rank']:03d}_{slug(e['title'])}.pdf"
        if os.path.exists(fn) and os.path.getsize(fn) > 20000:
            e['file'] = fn; e['status'] = 'ok'; continue
        ok, why = False, e.get('status')
        for u in routes(e, client):
            ok, why = try_get(u, fn, client, sleep=2.5)
            if ok: break
        if ok:
            e['file'] = fn; e['status'] = why; fixed += 1
        else:
            e['status'] = why
        print(f"{e['rank']:3d} {'OK  ' if ok else 'FAIL'} {str(why)[:18]:20s} {e['title'][:70]}", flush=True)
json.dump(man, open('paper_search_results/download_manifest.json','w'), indent=1)
print(f'\nnewly fixed: {fixed}; total available: {sum(1 for e in man if e.get("file"))}/{len(man)}')
