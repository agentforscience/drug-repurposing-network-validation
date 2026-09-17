#!/usr/bin/env python
"""Literature search across Europe PMC, Semantic Scholar, arXiv. Incremental + resumable."""
import json, time, sys, os, re, urllib.parse, traceback
import httpx

OUT = "paper_search_results"
os.makedirs(OUT, exist_ok=True)
RAW = open(f"{OUT}/raw_hits.jsonl", "a", buffering=1)

QUERIES = [
    "network proximity drug target disease module drug repurposing",
    "network-based drug repurposing protein-protein interaction proximity",
    "drug-disease network proximity clinical trial outcome validation",
    "network medicine uncovering disease-disease relationships interactome",
    "protein interaction network incompleteness bias evaluation network medicine",
    "tissue-specific protein interaction network drug target disease module",
    "computational drug repurposing benchmark evaluation AUROC negative samples",
    "drug repurposing prediction machine learning drug-disease association benchmark",
    "network proximity COVID-19 drug repurposing interactome",
    "degree bias hubs network proximity confound drug target",
    "repoDB gold standard drug repositioning approved failed indications",
    "drug target overlap disease genes simple baseline repurposing prediction",
    "critical assessment evaluation pitfalls drug repurposing prediction models",
    "separation measure disease module overlap interactome Menche",
]

def log(*a):
    print(*a, flush=True)

def emit(recs):
    for r in recs:
        RAW.write(json.dumps(r) + "\n")

def epmc(query, page_size=50):
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    params = {"query": query, "format": "json", "pageSize": page_size,
              "resultType": "core", "sort": "CITED desc"}
    r = httpx.get(url, params=params, timeout=120)
    r.raise_for_status()
    res = r.json().get("resultList", {}).get("result", [])
    out = []
    for p in res:
        out.append({
            "source": "europepmc",
            "id": p.get("id"), "pmid": p.get("pmid"), "pmcid": p.get("pmcid"),
            "doi": p.get("doi"), "title": p.get("title"),
            "authors": p.get("authorString"), "year": p.get("pubYear"),
            "journal": (p.get("journalInfo") or {}).get("journal", {}).get("title"),
            "citations": p.get("citedByCount"),
            "abstract": p.get("abstractText"),
            "isOpenAccess": p.get("isOpenAccess"), "hasPDF": p.get("hasPDF"),
            "query": query,
        })
    return out

def s2(query, limit=25, retries=4):
    url = "https://api.semanticscholar.org/graph/v1/paper/search"
    fields = "title,abstract,year,citationCount,authors,externalIds,openAccessPdf,venue"
    for i in range(retries):
        try:
            r = httpx.get(url, params={"query": query, "limit": limit, "fields": fields}, timeout=60)
            if r.status_code in (429, 503):
                time.sleep(8 * (i + 1)); continue
            r.raise_for_status()
            data = r.json().get("data", [])
            out = []
            for p in data:
                ext = p.get("externalIds") or {}
                out.append({
                    "source": "semanticscholar", "id": p.get("paperId"),
                    "doi": ext.get("DOI"), "pmid": ext.get("PubMed"), "arxiv": ext.get("ArXiv"),
                    "title": p.get("title"),
                    "authors": ", ".join(a.get("name","") for a in (p.get("authors") or [])[:8]),
                    "year": p.get("year"), "journal": p.get("venue"),
                    "citations": p.get("citationCount"), "abstract": p.get("abstract"),
                    "oa_pdf": (p.get("openAccessPdf") or {}).get("url"), "query": query,
                })
            return out
        except Exception:
            time.sleep(8 * (i + 1))
    return []

def arxiv(query, limit=12):
    q = "all:" + urllib.parse.quote(query)
    try:
        r = httpx.get(f"http://export.arxiv.org/api/query?search_query={q}&start=0&max_results={limit}&sortBy=relevance", timeout=60)
        r.raise_for_status()
    except Exception:
        return []
    entries = re.findall(r"<entry>(.*?)</entry>", r.text, re.S)
    out = []
    for e in entries:
        def g(tag):
            m = re.search(rf"<{tag}>(.*?)</{tag}>", e, re.S)
            return re.sub(r"\s+", " ", m.group(1)).strip() if m else None
        aid = (g("id") or "").rsplit("/", 1)[-1]
        out.append({"source": "arxiv", "id": aid, "arxiv": aid, "title": g("title"),
                    "abstract": g("summary"), "year": (g("published") or "")[:4],
                    "authors": ", ".join(re.findall(r"<name>(.*?)</name>", e)[:8]),
                    "query": query})
    return out

for qi, q in enumerate(QUERIES):
    log(f"[{qi+1}/{len(QUERIES)}] {q}")
    for name, fn in (("epmc", epmc), ("s2", s2), ("arxiv", arxiv)):
        try:
            recs = fn(q)
            emit(recs); log(f"   {name}: {len(recs)}")
        except Exception as ex:
            log(f"   {name} FAIL: {ex!r}")
        time.sleep(1.5)

RAW.close()

# dedupe
seen = {}
with open(f"{OUT}/raw_hits.jsonl") as f:
    for line in f:
        try: p = json.loads(line)
        except: continue
        t = re.sub(r"[^a-z0-9]", "", (p.get("title") or "").lower())
        if not t: continue
        if t in seen:
            old = seen[t]
            for k, v in p.items():
                if v and not old.get(k): old[k] = v
            old.setdefault("queries", set())
        else:
            seen[t] = p
uniq = list(seen.values())
with open(f"{OUT}/all_candidates.jsonl", "w") as f:
    for p in uniq: f.write(json.dumps(p) + "\n")
log(f"DONE unique={len(uniq)}")
