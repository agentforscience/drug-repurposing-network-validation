#!/usr/bin/env python
import json, time, re, urllib.parse
import httpx
OUT="paper_search_results"; RAW=open(f"{OUT}/raw_hits.jsonl","a",buffering=1)
QUERIES = [
 "limitations network-based drug repurposing criticism proximity",
 "systematic integration biomedical knowledge prioritizes drugs repurposing Rephetio",
 "topological bias network prediction protein interaction study bias",
 "degree-preserving null model network proximity significance z-score",
 "evaluation drug repositioning methods cross-validation bias gold standard negatives",
 "does network topology predict drug efficacy negative result",
 "interactome incompleteness impact disease gene prediction performance",
 "drug repurposing prediction simple baseline outperforms complex model",
 "network proximity separation drug disease Guney efficacy screening",
 "hub proteins bias drug targets promiscuity network analysis",
 "tissue-specific interactome disease module network proximity improvement",
 "benchmark comparison network based drug disease association prediction methods",
 "clinical trial failure prediction efficacy drug indication computational",
 "data leakage evaluation biomedical link prediction knowledge graph",
 "drug repurposing hub clinical phase drug indication dataset",
 "shortest path closest measure drug targets disease proteins AUC",
]
def epmc(q,ps=40):
    r=httpx.get("https://www.ebi.ac.uk/europepmc/webservices/rest/search",
        params={"query":q,"format":"json","pageSize":ps,"resultType":"core","sort":"CITED desc"},timeout=120)
    r.raise_for_status()
    out=[]
    for p in r.json().get("resultList",{}).get("result",[]):
        out.append({"source":"europepmc","id":p.get("id"),"pmid":p.get("pmid"),"pmcid":p.get("pmcid"),
            "doi":p.get("doi"),"title":p.get("title"),"authors":p.get("authorString"),"year":p.get("pubYear"),
            "journal":(p.get("journalInfo") or {}).get("journal",{}).get("title"),"citations":p.get("citedByCount"),
            "abstract":p.get("abstractText"),"isOpenAccess":p.get("isOpenAccess"),"hasPDF":p.get("hasPDF"),"query":q})
    return out
def s2(q,limit=25,retries=5):
    for i in range(retries):
        try:
            r=httpx.get("https://api.semanticscholar.org/graph/v1/paper/search",
                params={"query":q,"limit":limit,"fields":"title,abstract,year,citationCount,authors,externalIds,openAccessPdf,venue"},timeout=60)
            if r.status_code in (429,503): time.sleep(10*(i+1)); continue
            r.raise_for_status()
            out=[]
            for p in r.json().get("data",[]):
                ext=p.get("externalIds") or {}
                out.append({"source":"semanticscholar","id":p.get("paperId"),"doi":ext.get("DOI"),
                    "pmid":ext.get("PubMed"),"arxiv":ext.get("ArXiv"),"title":p.get("title"),
                    "authors":", ".join(a.get("name","") for a in (p.get("authors") or [])[:8]),
                    "year":p.get("year"),"journal":p.get("venue"),"citations":p.get("citationCount"),
                    "abstract":p.get("abstract"),"oa_pdf":(p.get("openAccessPdf") or {}).get("url"),"query":q})
            return out
        except Exception: time.sleep(10*(i+1))
    return []
for i,q in enumerate(QUERIES):
    print(f"[{i+1}/{len(QUERIES)}] {q}",flush=True)
    for nm,fn in (("epmc",epmc),("s2",s2)):
        try:
            recs=fn(q)
            for r in recs: RAW.write(json.dumps(r)+"\n")
            print(f"  {nm}: {len(recs)}",flush=True)
        except Exception as e: print(f"  {nm} FAIL {e!r}",flush=True)
        time.sleep(2)
print("DONE2",flush=True)
