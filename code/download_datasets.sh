#!/bin/bash
# Downloads core datasets. Idempotent-ish: skips non-empty existing files.
set -u
cd "$(dirname "$0")/.."
UA="Mozilla/5.0 (research; academic use)"
get(){ # url outpath
  if [ -s "$2" ]; then echo "SKIP $2"; return; fi
  mkdir -p "$(dirname "$2")"
  echo "GET $1 -> $2"
  curl -sL --retry 4 --retry-delay 5 --max-time 1800 -A "$UA" "$1" -o "$2" \
    && echo "  ok $(du -h "$2" | cut -f1)" || echo "  FAIL $1"
}

# --- Hetionet v1.0 (DrugBank-keyed drug-target, DOID-keyed disease-gene, PPI, CtD gold standard)
get "https://github.com/hetio/hetionet/raw/master/hetnet/tsv/hetionet-v1.0-edges.sif.gz" datasets/hetionet/hetionet-v1.0-edges.sif.gz
get "https://github.com/hetio/hetionet/raw/master/hetnet/tsv/hetionet-v1.0-nodes.tsv" datasets/hetionet/hetionet-v1.0-nodes.tsv

# --- Disease Ontology (DOID <-> UMLS CUI xrefs, to link repoDB CUIs to DOID resources)
get "https://github.com/DiseaseOntology/HumanDiseaseOntology/raw/main/src/ontology/doid.obo" datasets/ontology/doid.obo

# --- Jensen lab DISEASES (curated knowledge channel, DOID-keyed disease-gene)
get "https://download.jensenlab.org/human_disease_knowledge_filtered.tsv" datasets/diseases_jensen/human_disease_knowledge_filtered.tsv
get "https://download.jensenlab.org/human_disease_experiments_filtered.tsv" datasets/diseases_jensen/human_disease_experiments_filtered.tsv

# --- HuRI: systematic, study-bias-free human binary interactome (Luck et al. 2020)
get "https://www.interactome-atlas.org/data/HuRI.tsv" datasets/ppi/HuRI.tsv
get "https://www.interactome-atlas.org/data/HI-union.tsv" datasets/ppi/HI-union.tsv
get "https://www.interactome-atlas.org/data/Lit-BM.tsv" datasets/ppi/Lit-BM.tsv

# --- STRING v12 human (weighted, multi-evidence) + protein aliases for ID mapping
get "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz" datasets/ppi/9606.protein.links.v12.0.txt.gz
get "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz" datasets/ppi/9606.protein.info.v12.0.txt.gz

# --- BioGRID human physical interactions (literature-curated; study bias present)
get "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.tab3.zip" datasets/ppi/BIOGRID-ORGANISM-LATEST.tab3.zip

# --- DrugCentral (drug-target bioactivity + identifier crossrefs incl. DrugBank)
get "https://unmtid-shinyapps.net/download/DrugCentral/2021_09_01/drug.target.interaction.tsv.gz" datasets/drugcentral/drug.target.interaction.tsv.gz
get "https://unmtid-shinyapps.net/download/DrugCentral/2021_09_01/drugcentral.dump.010_05_2021.sql.gz" datasets/drugcentral/drugcentral_dump.sql.gz

# --- NCBI gene info (Entrez <-> symbol mapping, human)
get "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz" datasets/mappings/Homo_sapiens.gene_info.gz
echo "ALL DONE"
