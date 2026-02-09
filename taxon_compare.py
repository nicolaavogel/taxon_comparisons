#!/usr/bin/env python3

#######################################################################################
### importing libraries
import argparse
from pygbif import species
from Bio import Entrez
from dataclasses import dataclass
from pathlib import Path
from functools import lru_cache
from typing import Dict, Optional, Set, Tuple, List
import json
import gzip

#######################################################################################
### Functions for calling the taxonomy from GBIF 

def gbif_find_taxon_key(name: str, rank: str, prefer_kingdom: str = "Animalia"):
    """Find best GBIF taxon key for a name at a given rank."""
    rank = rank.upper()
    suggestions = species.name_suggest(q=name, rank=rank, limit=20)

    hits = [s for s in suggestions if s.get("rank") == rank and s.get("canonicalName", "").lower() == name.lower()]
    if not hits:
        # fallback: take any hit at that rank
        hits = [s for s in suggestions if s.get("rank") == rank]

    if not hits:
        raise ValueError(f"No GBIF match found for {rank}: {name}")

    # prefer a particular kingdom if possible
    preferred = [h for h in hits if h.get("kingdom") == prefer_kingdom]
    best = preferred[0] if preferred else hits[0]

    return best["key"], best


def gbif_list_species_under_taxon(
    taxon_name: str,
    taxon_rank: str,
    accepted_only: bool = True,
    limit_per_page: int = 300,
    return_records: bool = False,
    prefer_kingdom: str = "Animalia",
    include_subspecies: bool = False,
):
    taxon_key, match = gbif_find_taxon_key(taxon_name, taxon_rank, prefer_kingdom=prefer_kingdom)

    offset = 0
    all_results = []

    while True:
        params = {
            "higherTaxonKey": taxon_key,
            "rank": "SPECIES" if not include_subspecies else "SUBSPECIES",
            "limit": limit_per_page,
            "offset": offset
        }
        if accepted_only:
            params["status"] = "ACCEPTED"

        res = species.name_lookup(**params)
        results = res.get("results", [])
        all_results.extend(results)

        # Stop when we've retrieved all results
        retrieved = offset + len(results)
        if retrieved >= res.get("count", 0) or len(results) == 0:
            break

        offset += limit_per_page

    if return_records:
        return all_results, match

    # Return just species/subspecies names
    names = []
    for r in all_results:
        names.append(r.get("canonicalName") or r.get("scientificName"))
    return sorted(set(names)), match

    
#######################################################################################    
### NCBI call functions (references available for a taxon in all of NCBI nuccore)
def ncbi_setup(email: str, api_key: Optional[str] = None) -> None:
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

def ncbi_get_taxid_by_rank(name: str, rank: str) -> int:
    rank = rank.lower()
    term = f'{name}[Scientific Name] AND {rank}[Rank]'
    h = Entrez.esearch(db="taxonomy", term=term, retmode="xml")
    r = Entrez.read(h)
    ids = r.get("IdList", [])
    if not ids:
        raise ValueError(f"No NCBI taxid for {rank}: {name}")
    return int(ids[0])


def ncbi_count_species_under_taxid(taxid: int) -> int:
    term = f"txid{taxid}[Subtree] AND species[Rank]"
    h = Entrez.esearch(db="taxonomy", term=term, retmax=0, retmode="xml")
    r = Entrez.read(h)
    return int(r["Count"])

def ncbi_count_nuccore_records_for_taxid(taxid: int) -> int:
    term = f"txid{taxid}[Organism:exp]"
    h = Entrez.esearch(db="nuccore", term=term, retmax=0, retmode="xml")
    r = Entrez.read(h)
    return int(r["Count"])


def ncbi_list_species_under_taxid(taxid: int, batch_size: int = 500, include_subspecies: bool = False) -> Dict[int, str]:
    rank = "species" if not include_subspecies else "subspecies"
    term = f"txid{taxid}[Subtree] AND {rank}[Rank]"

    # First get count
    h0 = Entrez.esearch(db="taxonomy", term=term, retmax=0, retmode="xml")
    r0 = Entrez.read(h0)
    total = int(r0["Count"])

    out: Dict[int, str] = {}

    # Paginate through IDs
    for start in range(0, total, batch_size):
        h = Entrez.esearch(db="taxonomy", term=term, retstart=start, retmax=batch_size, retmode="xml")
        r = Entrez.read(h)
        ids = r.get("IdList", [])
        if not ids:
            continue

        # Fetch records for these taxonomy IDs
        fh = Entrez.efetch(db="taxonomy", id=",".join(ids), retmode="xml")
        fr = Entrez.read(fh)

        for rec in fr:
            try:
                sid = int(rec["TaxId"])
                nm = rec.get("ScientificName", "")
                if nm:
                    out[sid] = nm
            except Exception:
                continue

    return out


def ncbi_taxon_nuccore_summary(
    taxon_name: str,
    taxon_rank: str,
    include_species_list: bool = False,
) -> Dict[str, object]:
    taxid = ncbi_get_taxid_by_rank(taxon_name, taxon_rank)

    d: Dict[str, object] = {
        f"{taxon_rank.lower()}_taxid": taxid,
        "n_species_ncbi_taxonomy": ncbi_count_species_under_taxid(taxid),
        "n_nuccore_records_ncbi": ncbi_count_nuccore_records_for_taxid(taxid),
    }

    if include_species_list:
        # species only
        sp_map = ncbi_list_species_under_taxid(taxid, include_subspecies=False)
        # subspecies only
        ssp_map = ncbi_list_species_under_taxid(taxid, include_subspecies=True)

        d["ncbi_species_taxids"] = sorted(sp_map.keys())
        d["ncbi_species_names"] = sorted(set(sp_map.values()))
        d["ncbi_species_taxid_to_name"] = sp_map

        d["ncbi_subspecies_taxids"] = sorted(ssp_map.keys())
        d["ncbi_subspecies_names"] = sorted(set(ssp_map.values()))
        d["ncbi_subspecies_taxid_to_name"] = ssp_map

    return d




#######################################################################################
# class / function for call to internal DB acc2tax file (how many references are represented in our DB)

@dataclass
class TaxNode:
    parent: int
    rank: str

def load_ncbi_taxdump(nodes_dmp_path: str, names_dmp_path: Optional[str] = None) -> Tuple[Dict[int, TaxNode], Dict[int, str]]:
    nodes: Dict[int, TaxNode] = {}
    with open(nodes_dmp_path, "r", encoding="utf-8") as f:
        for line in f:
            parts = [p.strip() for p in line.split("|")]
            taxid = int(parts[0])
            parent = int(parts[1])
            rank = parts[2]
            nodes[taxid] = TaxNode(parent=parent, rank=rank)

    names: Dict[int, str] = {}
    if names_dmp_path:
        with open(names_dmp_path, "r", encoding="utf-8") as f:
            for line in f:
                parts = [p.strip() for p in line.split("|")]
                taxid = int(parts[0])
                name_txt = parts[1]
                name_class = parts[3]
                if name_class == "scientific name":
                    names[taxid] = name_txt

    return nodes, names

def lineage_to_root(taxid: int, nodes: Dict[int, TaxNode], max_steps: int = 10_000) -> List[int]:
    lin = []
    cur = taxid
    steps = 0
    while True:
        lin.append(cur)
        if cur == 1:
            break
        node = nodes.get(cur)
        if node is None:
            break
        nxt = node.parent
        if nxt == cur:
            break
        cur = nxt
        steps += 1
        if steps > max_steps:
            raise RuntimeError(f"Lineage too deep / loop for taxid={taxid}")
    return lin

def make_ancestor_fns(nodes: Dict[int, TaxNode]):
    @lru_cache(maxsize=2_000_000)
    def ancestor_cached(taxid: int, rank: str) -> Optional[int]:
        rank = rank.lower()
        cur = taxid
        steps = 0
        while True:
            n = nodes.get(cur)
            if n is None:
                return None
            if n.rank.lower() == rank:
                return cur
            if cur == 1 or n.parent == cur:
                return None
            cur = n.parent
            steps += 1
            if steps > 10_000:
                return None

    return ancestor_cached

def ancestor_at_rank(taxid: int, nodes: Dict[int, TaxNode], target_rank: str) -> Optional[int]:
    target_rank = target_rank.lower()
    for t in lineage_to_root(taxid, nodes):
        n = nodes.get(t)
        if n and n.rank.lower() == target_rank:
            return t
    return None

def iter_acc2tax_4col(acc2tax_path: str, default_db_id: str = "UNKNOWN"):
    """
    Yields (acc_unver, acc_ver, taxid, db_id).

    Expects whitespace-separated columns with >= 3 cols.
    If a 4th column exists, it is used as db_id.
    """
    opener = gzip.open if acc2tax_path.endswith(".gz") else open
    with opener(acc2tax_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            # skip header-ish rows
            if parts[2].lower() == "taxid" or not parts[2].isdigit():
                continue

            acc_unver = parts[0]
            acc_ver = parts[1]
            taxid = int(parts[2])
            db_id = parts[3] if len(parts) >= 4 else default_db_id
            yield acc_unver, acc_ver, taxid, db_id


def local_db_details_for_taxon_by_dbid(
    taxon_taxid: int,
    taxon_rank: str,
    acc2tax_rows,
    nodes: Dict[int, TaxNode],
    names: Optional[Dict[int, str]] = None,
    count_by: str = "versioned",
) -> Dict[str, object]:
    """
    Returns:
      - overall union across all db_ids
      - per-db_id record/species/subspecies sets
    """
    taxon_rank = taxon_rank.lower()
    if count_by not in ("versioned", "unversioned"):
        raise ValueError("count_by must be 'versioned' or 'unversioned'")

    ancestor_cached = make_ancestor_fns(nodes)

    # overall
    all_record_ids: Set[str] = set()
    all_species_ids: Set[int] = set()
    all_subspecies_ids: Set[int] = set()

    # per db_id
    per_db: Dict[str, Dict[str, Set]] = {}

    for acc_unver, acc_ver, taxid, db_id in acc2tax_rows:
        anc = ancestor_cached(taxid, taxon_rank)
        if anc != taxon_taxid:
            continue

        rec = acc_ver if count_by == "versioned" else acc_unver

        if db_id not in per_db:
            per_db[db_id] = {
                "record_ids": set(),
                "species_ids": set(),
                "subspecies_ids": set(),
            }

        per_db[db_id]["record_ids"].add(rec)
        all_record_ids.add(rec)

        sp = ancestor_cached(taxid, "species")
        if sp is not None:
            per_db[db_id]["species_ids"].add(sp)
            all_species_ids.add(sp)

        ssp = ancestor_cached(taxid, "subspecies")
        if ssp is not None:
            per_db[db_id]["subspecies_ids"].add(ssp)
            all_subspecies_ids.add(ssp)

    def ids_to_names(ids: List[int]) -> List[str]:
        if not names:
            return [f"taxid:{t}" for t in ids]
        return [names.get(t, f"taxid:{t}") for t in ids]

    # finalize per-db (sorted + names)
    per_db_out: Dict[str, object] = {}
    for db_id, d in per_db.items():
        sp_ids = sorted(d["species_ids"])
        ssp_ids = sorted(d["subspecies_ids"])

        per_db_out[db_id] = {
            "n_records_in_taxon": len(d["record_ids"]),
            "n_species_in_taxon": len(d["species_ids"]),
            "n_subspecies_in_taxon": len(d["subspecies_ids"]),
            "species_taxids": sp_ids,
            "subspecies_taxids": ssp_ids,
            "species_names": sorted(set(ids_to_names(sp_ids))),
            "subspecies_names": sorted(set(ids_to_names(ssp_ids))),
        }

    # overall union
    all_species_taxids = sorted(all_species_ids)
    all_subspecies_taxids = sorted(all_subspecies_ids)

    overall = {
        "n_records_in_taxon": len(all_record_ids),
        "n_species_in_taxon": len(all_species_ids),
        "n_subspecies_in_taxon": len(all_subspecies_ids),
        "species_taxids": all_species_taxids,
        "subspecies_taxids": all_subspecies_taxids,
        "species_names": sorted(set(ids_to_names(all_species_taxids))) if names else None,
        "subspecies_names": sorted(set(ids_to_names(all_subspecies_taxids))) if names else None,
    }

    return {"overall": overall, "per_db": per_db_out}



    
        
###################
# COMPARISON SECTION 

def normalize_name(s: str) -> str:
    return " ".join(s.strip().split()).lower()

def compare_species_sets(
    gbif_names: List[str],
    ncbi_names: Optional[List[str]],
    local_names: Optional[List[str]],
) -> Dict[str, object]:
    gbif = {normalize_name(x) for x in (gbif_names or [])}
    ncbi = {normalize_name(x) for x in (ncbi_names or [])} if ncbi_names else set()
    local = {normalize_name(x) for x in (local_names or [])} if local_names else set()

    return {
        "overlap_all_3": sorted(gbif & ncbi & local),
        "gbif_only": sorted(gbif - ncbi - local),
        "ncbi_only": sorted(ncbi - gbif - local),
        "local_only": sorted(local - gbif - ncbi),
        "gbif_ncbi_overlap": sorted(gbif & ncbi),
        "gbif_local_overlap": sorted(gbif & local),
        "ncbi_local_overlap": sorted(ncbi & local),
        "n_overlap_all_3": len(gbif & ncbi & local),
        "n_gbif_only": len(gbif - ncbi - local),
        "n_ncbi_only": len(ncbi - gbif - local),
        "n_local_only": len(local - gbif - ncbi),
    }
    


#######################################################################################
### CALL ALL FUNCTIONS FOR COMAPRISON
def compare_taxon_gbif_ncbi_local(
    taxon_name: str,
    taxon_rank: str,
    gbif_species_list: List[str],
    gbif_subspecies_list: List[str],
    acc2tax_path: str,
    nodes_dmp_path: str,
    names_dmp_path: Optional[str] = None,
    email_for_ncbi: str = "your.email@domain.com",
    api_key: Optional[str] = None,
    count_by: str = "versioned",
    include_subspecies: bool = False,
):
    taxon_rank_upper = taxon_rank.upper()

    # NCBI online
    ncbi_setup(email=email_for_ncbi, api_key=api_key)
    ncbi_summary = ncbi_taxon_nuccore_summary(
        taxon_name,
        taxon_rank_upper,
        include_species_list=True,
    )
    

    # local DB
    nodes, names = load_ncbi_taxdump(nodes_dmp_path, names_dmp_path)
    rows = iter_acc2tax_4col(acc2tax_path)

    taxid_key = f"{taxon_rank_upper.lower()}_taxid"
    local = local_db_details_for_taxon_by_dbid(
        taxon_taxid=ncbi_summary[taxid_key],
        taxon_rank=taxon_rank_upper,
        acc2tax_rows=rows,
        nodes=nodes,
        names=names,
        count_by=count_by,
    )

    # Build a unified output model with explicit species/subspecies lists per source
    # GBIF
    gbif_species_names = gbif_species_list or []
    gbif_subspecies_names = gbif_subspecies_list or []

    # NCBI: ncbi_list_species_under_taxid returns (species OR subspecies) depending on include_subspecies flag,
    # so to do BOTH, call twice.
    taxid = ncbi_summary[taxid_key]
    ncbi_species_map = ncbi_list_species_under_taxid(taxid, include_subspecies=False)
    ncbi_subspecies_map = ncbi_list_species_under_taxid(taxid, include_subspecies=True)

    ncbi_species_names = sorted(set(ncbi_species_map.values()))
    ncbi_subspecies_names = sorted(set(ncbi_subspecies_map.values()))

    # Local overall
    local_overall = local["overall"]
    local_species_names = local_overall.get("species_names") or []
    local_subspecies_names = local_overall.get("subspecies_names") or []

    # Overall comparison (combined sets)
    overall_comparison = compare_species_sets(
        gbif_names=gbif_species_names + gbif_subspecies_names,
        ncbi_names=ncbi_species_names + ncbi_subspecies_names,
        local_names=local_species_names + local_subspecies_names,
    )

    # Per-db comparisons
    per_db_comparisons: Dict[str, object] = {}
    for db_id, dbd in local["per_db"].items():
        per_db_comparisons[db_id] = compare_species_sets(
            gbif_names=gbif_species_names + gbif_subspecies_names,
            ncbi_names=ncbi_species_names + ncbi_subspecies_names,
            local_names=(dbd.get("species_names") or []) + (dbd.get("subspecies_names") or []),
        )

    return {
        "taxon": taxon_name,
        "taxon_rank": taxon_rank_upper,

        # GBIF
        "gbif_species_names": gbif_species_names,
        "gbif_subspecies_names": gbif_subspecies_names,
        "gbif_n_species": len(gbif_species_names),
        "gbif_n_subspecies": len(gbif_subspecies_names),

        # NCBI summary counts
        **ncbi_summary,

        # NCBI explicit lists (both ranks)
        "ncbi_species_names": ncbi_species_names,
        "ncbi_subspecies_names": ncbi_subspecies_names,
        "ncbi_n_species": len(ncbi_species_names),
        "ncbi_n_subspecies": len(ncbi_subspecies_names),

        # Local overall + per-db
        "local_count_by": count_by,
        "local_overall": local_overall,
        "local_per_db": local["per_db"],

        # Comparisons
        "comparison_overall_combined": overall_comparison,
        "comparison_per_db_combined": per_db_comparisons,
    }

def build_species_buckets(res: dict) -> dict:
    """
    Returns disjoint buckets (normalized names) for GBIF/NCBI/LOCAL
    using combined species+subspecies for each source.
    """
    gbif_raw = (res.get("gbif_species_names") or []) + (res.get("gbif_subspecies_names") or [])
    ncbi_raw = (res.get("ncbi_species_names") or []) + (res.get("ncbi_subspecies_names") or [])

    local_overall = res.get("local_overall") or {}
    local_raw = (local_overall.get("species_names") or []) + (local_overall.get("subspecies_names") or [])

    gbif = {normalize_name(x) for x in gbif_raw if str(x).strip()}
    ncbi = {normalize_name(x) for x in ncbi_raw if str(x).strip()}
    local = {normalize_name(x) for x in local_raw if str(x).strip()}

    return {
        "ALL_3": gbif & ncbi & local,
        "GBIF_NCBI": (gbif & ncbi) - local,
        "GBIF_LOCAL": (gbif & local) - ncbi,
        "NCBI_LOCAL": (ncbi & local) - gbif,
        "GBIF_ONLY": gbif - ncbi - local,
        "NCBI_ONLY": ncbi - gbif - local,
        "LOCAL_ONLY": local - gbif - ncbi,
    }


def summarize_buckets(buckets: dict) -> dict:
    return {k: len(v) for k, v in buckets.items()}

def write_bucket_lists(
    buckets: dict,
    out_dir: str = ".",
    prefix: str = "species",
) -> None:
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    for k, v in buckets.items():
        p = Path(out_dir) / f"{prefix}_{k}.txt"
        with open(p, "w", encoding="utf-8") as f:
            f.write("\n".join(sorted(v)))
            
def build_buckets_for_sets(gbif_list: List[str], ncbi_list: List[str], local_list: List[str]) -> dict:
    gbif = {normalize_name(x) for x in (gbif_list or []) if str(x).strip()}
    ncbi = {normalize_name(x) for x in (ncbi_list or []) if str(x).strip()}
    local = {normalize_name(x) for x in (local_list or []) if str(x).strip()}

    return {
        "ALL_3": gbif & ncbi & local,
        "GBIF_NCBI": (gbif & ncbi) - local,
        "GBIF_LOCAL": (gbif & local) - ncbi,
        "NCBI_LOCAL": (ncbi & local) - gbif,
        "GBIF_ONLY": gbif - ncbi - local,
        "NCBI_ONLY": ncbi - gbif - local,
        "LOCAL_ONLY": local - gbif - ncbi,
    }

def write_bucket_lists_per_db(res: dict, out_dir: str = ".", prefix: str = "species") -> None:
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    gbif_all = (res.get("gbif_species_names", []) or []) + (res.get("gbif_subspecies_names", []) or [])
    ncbi_all = (res.get("ncbi_species_names", []) or []) + (res.get("ncbi_subspecies_names", []) or [])

    # overall local
    local_overall = res.get("local_overall", {})
    local_all = (local_overall.get("species_names") or []) + (local_overall.get("subspecies_names") or [])
    buckets_overall = build_buckets_for_sets(gbif_all, ncbi_all, local_all)
    write_bucket_lists(buckets_overall, out_dir=out_dir, prefix=f"{prefix}_LOCAL-ALL")

    # per db_id
    for db_id, d in (res.get("local_per_db") or {}).items():
        local_db_all = (d.get("species_names") or []) + (d.get("subspecies_names") or [])
        b = build_buckets_for_sets(gbif_all, ncbi_all, local_db_all)
        safe_db = "".join(c if c.isalnum() or c in "-_." else "_" for c in str(db_id))
        write_bucket_lists(b, out_dir=out_dir, prefix=f"{prefix}_DB-{safe_db}")

def write_species_membership_html(res: dict, out_html: str = "species_membership.html") -> None:
    gbif_raw = (res.get("gbif_species_names", []) or []) + (res.get("gbif_subspecies_names", []) or [])
    ncbi_raw = (res.get("ncbi_species_names", []) or []) + (res.get("ncbi_subspecies_names", []) or [])

    local_per_db = res.get("local_per_db") or {}
    db_ids = sorted(local_per_db.keys())

    def to_map(lst):
        m = {}
        for x in lst:
            nx = normalize_name(x)
            if nx and nx not in m:
                m[nx] = str(x).strip()
        return m

    gbif_map = to_map(gbif_raw)
    ncbi_map = to_map(ncbi_raw)

    local_maps = {}
    for db_id in db_ids:
        d = local_per_db[db_id]
        local_raw = (d.get("species_names") or []) + (d.get("subspecies_names") or [])
        local_maps[db_id] = to_map(local_raw)

    all_keys = sorted(set(gbif_map) | set(ncbi_map) | set().union(*(set(m) for m in local_maps.values())))

    rows = []
    for k in all_keys:
        pretty = gbif_map.get(k) or ncbi_map.get(k) or next((local_maps[d].get(k) for d in db_ids if k in local_maps[d]), k) or k
        row = {
            "name": pretty,
            "key": k,
            "GBIF": k in gbif_map,
            "NCBI": k in ncbi_map,
        }
        present_dbs = []
        for db_id in db_ids:
            in_db = k in local_maps[db_id]
            row[f"LOCAL::{db_id}"] = in_db
            if in_db:
                present_dbs.append(str(db_id))
        row["LOCAL_DBS"] = ", ".join(present_dbs)
        rows.append(row)

    # Build HTML
    taxon_title = f'{res.get("taxon","Taxon")} {res.get("taxon_rank","")}'.strip()

    # header with dynamic columns
    local_th = "\n".join([f'        <th onclick="sortTable({3+i})">{db}</th>' for i, db in enumerate(db_ids)])
    # columns: Species, GBIF, NCBI, (each db), Local DBs, Bucket(Overall Local)
    # Bucket here is computed vs "any local db" (overall local union)
    html = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8"/>
  <title>{taxon_title} membership</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 18px; }}
    input {{ padding: 8px; width: 520px; }}
    table {{ border-collapse: collapse; margin-top: 12px; width: 100%; }}
    th, td {{ border: 1px solid #ddd; padding: 8px; }}
    th {{ cursor: pointer; background: #f6f6f6; position: sticky; top: 0; }}
    tr:nth-child(even):not(.badge-gbif-only):not(.badge-strong):not(.badge-all):not(.badge-warn):not(.badge-only) {{
      background: #fafafa;
    }}
    
    .mono {{ font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace; }}

    /* row highlights */
    .badge-gbif-only {{ background: #ffecec !important; }}
    .badge-strong    {{ background: #dff6df !important; }}
    .badge-all       {{ background: #eefbee !important; }}
    .badge-warn      {{ background: #fff7cc !important; }}
    .badge-only      {{ background: #f2f2f2 !important; }}
  </style>
</head>
<body>
  <h2>{taxon_title} — membership (GBIF / NCBI / Local DBs)</h2>
  <div>
    <input id="q" placeholder="Search..." onkeyup="filterRows()"/>
    <span id="count"></span>
  </div>

  <table id="tbl">
    <thead>
      <tr>
        <th onclick="sortTable(0)">Name</th>
        <th onclick="sortTable(1)">GBIF</th>
        <th onclick="sortTable(2)">NCBI</th>
{local_th}
        <th onclick="sortTable({3+len(db_ids)})">Local DBs</th>
        <th onclick="sortTable({4+len(db_ids)})">Bucket (GBIF/NCBI/Any Local)</th>
      </tr>
    </thead>
    <tbody>
"""

    # overall local union for bucket calc
    overall_local = res.get("local_overall") or {}
    overall_local_raw = (overall_local.get("species_names") or []) + (overall_local.get("subspecies_names") or [])
    overall_local_set = {normalize_name(x) for x in overall_local_raw if str(x).strip()}

    def bucket_of(k: str, g: bool, n: bool) -> str:
        l = k in overall_local_set
        if g and n and l: return "ALL_3"
        if g and n and not l: return "GBIF_NCBI"
        if g and l and not n: return "GBIF_LOCAL"
        if n and l and not g: return "NCBI_LOCAL"
        if g and not n and not l: return "GBIF_ONLY"
        if n and not g and not l: return "NCBI_ONLY"
        if l and not g and not n: return "LOCAL_ONLY"
        return "?"

    for r in rows:
        k = r["key"]
        bucket = bucket_of(k, r["GBIF"], r["NCBI"])

        # Which local DBs contain this name?
        present_dbs = {db_id for db_id in db_ids if r.get(f"LOCAL::{db_id}")}

        # Helper flags for your named DBs (adjust matching if your IDs differ)
        has_core_nt = any("core_nt" in str(d).lower() for d in present_dbs)
        has_wgs     = any("wgs"     in str(d).lower() for d in present_dbs)
        has_refseq  = any("refseq"  in str(d).lower() for d in present_dbs)

        # Strong green rule:
        # present in GBIF+NCBI+Local, AND (in all local DBs OR at least in wgs/refseq)
        strong_green = (
            bucket == "ALL_3"
            and (present_dbs == set(db_ids) or has_wgs or has_refseq)
        )

        # Light yellow rule:
        # (A) present in ALL_3 but only local in core_nt
        # OR
        # (B) present in NCBI+Local but only local in core_nt
        only_core_nt_local = (len(present_dbs) > 0) and all("core_nt" in str(d).lower() for d in present_dbs)
        light_yellow = (
            (bucket == "ALL_3" and only_core_nt_local)
            or (bucket == "NCBI_LOCAL" and only_core_nt_local)
        )

        # GBIF-only stays red
        if bucket == "GBIF_ONLY":
            cls = "badge-gbif-only"
        elif strong_green:
            cls = "badge-strong"
        elif light_yellow:
            cls = "badge-warn"
        elif bucket == "ALL_3":
            cls = "badge-all"
        elif bucket.endswith("_ONLY"):
            cls = "badge-only"
        else:
            cls = ""
        

        html += f"""      <tr class="{cls}">
        <td>{r["name"]}</td>
        <td>{"✔" if r["GBIF"] else ""}</td>
        <td>{"✔" if r["NCBI"] else ""}</td>
"""
        for db_id in db_ids:
            html += f"""        <td>{"✔" if r.get(f"LOCAL::{db_id}") else ""}</td>
"""
        html += f"""        <td class="mono">{r["LOCAL_DBS"]}</td>
        <td>{bucket}</td>
      </tr>
"""

    html += """    </tbody>
  </table>

<script>
function filterRows() {
  const q = document.getElementById("q").value.toLowerCase();
  const tbl = document.getElementById("tbl").getElementsByTagName("tbody")[0];
  const rows = tbl.getElementsByTagName("tr");
  let shown = 0;
  for (let i = 0; i < rows.length; i++) {
    const txt = rows[i].innerText.toLowerCase();
    const ok = txt.indexOf(q) > -1;
    rows[i].style.display = ok ? "" : "none";
    if (ok) shown++;
  }
  document.getElementById("count").innerText = "Showing " + shown + " rows";
}
filterRows();

function sortTable(col) {
  const table = document.getElementById("tbl");
  const tbody = table.tBodies[0];
  const rows = Array.from(tbody.rows);
  const asc = table.getAttribute("data-sort-col") != col || table.getAttribute("data-sort-dir") != "asc";

  rows.sort((a,b) => {
    const A = a.cells[col].innerText;
    const B = b.cells[col].innerText;
    return asc ? A.localeCompare(B) : B.localeCompare(A);
  });

  rows.forEach(r => tbody.appendChild(r));
  table.setAttribute("data-sort-col", col);
  table.setAttribute("data-sort-dir", asc ? "asc" : "desc");
}
</script>
</body></html>
"""
    with open(out_html, "w", encoding="utf-8") as f:
        f.write(html)

def load_config(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        cfg = json.load(f)
    if not isinstance(cfg, dict):
        raise ValueError("Config file must contain a JSON object at the top level.")
    return cfg



def main(
    taxon: str = "Alca",
    taxon_rank: str = "GENUS",
    acc2tax: str = "/path/to/acc2tax.tsv",
    nodes: str = "/path/to/nodes.dmp",
    names: Optional[str] = None,
    email: str = "your.email@domain.com",
    api_key: Optional[str] = None,
    kingdom: str = "Plantae",
    count_by: str = "versioned",
    pretty: bool = True,
    include_subspecies: bool = False,
):
    gbif_species_list, gbif_match = gbif_list_species_under_taxon(
        taxon_name=taxon,
        taxon_rank=taxon_rank,
        accepted_only=True,
        prefer_kingdom=kingdom,
        include_subspecies=False,
    )

    gbif_subspecies_list, _ = gbif_list_species_under_taxon(
        taxon_name=taxon,
        taxon_rank=taxon_rank,
        accepted_only=True,
        prefer_kingdom=kingdom,
        include_subspecies=True,
    )
    

    res = compare_taxon_gbif_ncbi_local(
        taxon_name=taxon,
        taxon_rank=taxon_rank,
        gbif_species_list=gbif_species_list,
        gbif_subspecies_list=gbif_subspecies_list,
        acc2tax_path=acc2tax,
        nodes_dmp_path=nodes,
        names_dmp_path=names,
        email_for_ncbi=email,
        api_key=api_key,
        count_by=count_by,
        include_subspecies=True,
    )
    

    if pretty:
        print(json.dumps(res, indent=2, sort_keys=True))
    else:
        print(json.dumps(res))

    return res


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Taxon comparison: GBIF vs NCBI vs local DB")
    parser.add_argument("--config", default="run_config.json", help="Path to JSON config file")
    args = parser.parse_args()

    cfg = load_config(args.config)

    # pull optional outputs from config (and remove them before calling main)
    out_dir = cfg.pop("out_dir", ".")
    html_out = cfg.pop("html_out", None)

    # run main using config values
    res = main(**cfg)

    # post-processing (same as you had, but configurable)
    buckets = build_species_buckets(res)
    res["species_buckets_counts"] = summarize_buckets(buckets)

    taxon = res.get("taxon", cfg.get("taxon", "Taxon"))
    taxon_rank = res.get("taxon_rank", cfg.get("taxon_rank", "RANK"))

    write_bucket_lists_per_db(res, out_dir=out_dir, prefix=f"{taxon}_{taxon_rank}")

    if html_out is None:
        html_out = f"{taxon}_{taxon_rank}_membership_by_dbid.html"
    write_species_membership_html(res, out_html=html_out)

   

