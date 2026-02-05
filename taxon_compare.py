#!/usr/bin/env python3

#######################################################################################
### importing libraries
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
    prefer_kingdom: str = "Animalia"
):
    """
    List species under a taxon (genus/family/etc.) in GBIF.

    Parameters
    ----------
    taxon_name : str
        e.g., "Alca"
    taxon_rank : str
        e.g., "GENUS" or "FAMILY"
    accepted_only : bool
        If True, only accepted species.
    limit_per_page : int
        GBIF paging size. Max ~300 is usually safe.
    return_records : bool
        If True, return full GBIF result dicts rather than just names.
    """
    taxon_key, match = gbif_find_taxon_key(taxon_name, taxon_rank, prefer_kingdom=prefer_kingdom)

    offset = 0
    all_results = []

    while True:
        params = {
            "higherTaxonKey": taxon_key,
            "rank": "SPECIES",
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
    
    # Return just species names (canonical if available, else scientificName)
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


def ncbi_list_species_under_taxid(taxid: int, batch_size: int = 500) -> Dict[int, str]:
    """
    Returns dict {species_taxid: scientific_name} for all species under a subtree taxid.
    Uses NCBI taxonomy esearch + efetch.
    """
    term = f"txid{taxid}[Subtree] AND species[Rank]"

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
            # rec is a dict with fields like 'TaxId' and 'ScientificName'
            try:
                sid = int(rec["TaxId"])
                nm = rec.get("ScientificName", "")
                if nm:
                    out[sid] = nm
            except Exception:
                continue

    return out

def ncbi_taxon_nuccore_summary(taxon_name: str, taxon_rank: str, include_species_list: bool = False) -> Dict[str, object]:
    taxid = ncbi_get_taxid_by_rank(taxon_name, taxon_rank)

    d: Dict[str, object] = {
        f"{taxon_rank.lower()}_taxid": taxid,
        "n_species_ncbi_taxonomy": ncbi_count_species_under_taxid(taxid),
        "n_nuccore_records_ncbi": ncbi_count_nuccore_records_for_taxid(taxid),
    }

    if include_species_list:
        sp = ncbi_list_species_under_taxid(taxid)
        d["ncbi_species_taxids"] = sorted(sp.keys())
        d["ncbi_species_names"] = sorted(set(sp.values()))
        d["ncbi_species_taxid_to_name"] = sp

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

def iter_acc2tax_3col(acc2tax_path: str):
    opener = gzip.open if acc2tax_path.endswith(".gz") else open
    with opener(acc2tax_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            if parts[2].lower() == "taxid" or not parts[2].isdigit():
                continue
            yield parts[0], parts[1], int(parts[2])

def local_db_details_for_taxon(
    taxon_taxid: int,
    taxon_rank: str,
    acc2tax_rows: List[Tuple[str, str, int]],
    nodes: Dict[int, TaxNode],
    names: Optional[Dict[int, str]] = None,
    count_by: str = "versioned",
) -> Dict[str, object]:
    taxon_rank = taxon_rank.lower()

    if count_by not in ("versioned", "unversioned"):
        raise ValueError("count_by must be 'versioned' or 'unversioned'")

    record_ids: Set[str] = set()
    species_ids: Set[int] = set()
    ancestor_cached = make_ancestor_fns(nodes)

    for acc_unver, acc_ver, taxid in acc2tax_rows:
        anc = ancestor_cached(taxid, taxon_rank)
        if anc != taxon_taxid:
            continue

        rec = acc_ver if count_by == "versioned" else acc_unver
        record_ids.add(rec)

        sp = ancestor_cached(taxid, "species")
        if sp is not None:
            species_ids.add(sp)

    species_taxids = sorted(species_ids)
    species_names = None
    if names:
        species_names = sorted({names.get(t, f"taxid:{t}") for t in species_taxids})

    return {
        "n_records_in_taxon": len(record_ids),
        "n_species_in_taxon": len(species_ids),
        "species_taxids": species_taxids,
        "species_names": species_names,
    }


    
        
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
    taxon_rank: str,          # "GENUS" or "FAMILY" etc.
    gbif_species_list: List[str],
    acc2tax_path: str,
    nodes_dmp_path: str,
    names_dmp_path: Optional[str] = None,
    email_for_ncbi: str = "n.alexandra.vogel@gmail.com",
    api_key: Optional[str] = None,
    count_by: str = "versioned",
):
    taxon_rank_upper = taxon_rank.upper()

    # NCBI online
    ncbi_setup(email=email_for_ncbi, api_key=api_key)
    ncbi_summary = ncbi_taxon_nuccore_summary(taxon_name, taxon_rank_upper, include_species_list=True)

    # local DB
    nodes, names = load_ncbi_taxdump(nodes_dmp_path, names_dmp_path)
    rows = iter_acc2tax_3col(acc2tax_path)

    taxid_key = f"{taxon_rank_upper.lower()}_taxid"
    local_details = local_db_details_for_taxon(
        taxon_taxid=ncbi_summary[taxid_key],
        taxon_rank=taxon_rank_upper,
        acc2tax_rows=rows,
        nodes=nodes,
        names=names,
        count_by=count_by
    )

    comparison = compare_species_sets(
        gbif_names=gbif_species_list,
        ncbi_names=ncbi_summary.get("ncbi_species_names"),
        local_names=local_details.get("species_names"),
    )

    return {
        "taxon": taxon_name,
        "taxon_rank": taxon_rank_upper,

        # GBIF
        "gbif_accepted_species_count": len(gbif_species_list),
        "gbif_accepted_species_names": gbif_species_list,

        # NCBI
        **ncbi_summary,

        # Local
        "local_n_records_in_taxon": local_details["n_records_in_taxon"],
        "local_n_species_in_taxon": local_details["n_species_in_taxon"],
        "local_species_taxids": local_details["species_taxids"],
        "local_species_names": local_details["species_names"],
        "local_count_by": count_by,

        "species_set_comparison": comparison,
    }

def build_species_buckets(res: dict) -> dict:
    """
    Returns disjoint buckets (normalized names) for GBIF/NCBI/LOCAL.
    """
    gbif_raw = res.get("gbif_accepted_species_names", []) or []
    ncbi_raw = res.get("ncbi_species_names", []) or []
    local_raw = res.get("local_species_names", []) or []

    gbif = {normalize_name(x) for x in gbif_raw if str(x).strip()}
    ncbi = {normalize_name(x) for x in ncbi_raw if str(x).strip()}
    local = {normalize_name(x) for x in local_raw if str(x).strip()}

    buckets = {
        "ALL_3": gbif & ncbi & local,
        "GBIF_NCBI": (gbif & ncbi) - local,
        "GBIF_LOCAL": (gbif & local) - ncbi,
        "NCBI_LOCAL": (ncbi & local) - gbif,
        "GBIF_ONLY": gbif - ncbi - local,
        "NCBI_ONLY": ncbi - gbif - local,
        "LOCAL_ONLY": local - gbif - ncbi,
    }
    return buckets

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
def write_species_membership_html(res: dict, out_html: str = "species_membership.html") -> None:
    gbif_raw = res.get("gbif_accepted_species_names", []) or []
    ncbi_raw = res.get("ncbi_species_names", []) or []
    local_raw = res.get("local_species_names", []) or []

    # Keep a "pretty" name representative for each normalized key
    def to_map(lst):
        m = {}
        for x in lst:
            nx = normalize_name(x)
            if nx and nx not in m:
                m[nx] = str(x).strip()
        return m

    gbif_map = to_map(gbif_raw)
    ncbi_map = to_map(ncbi_raw)
    local_map = to_map(local_raw)

    all_keys = sorted(set(gbif_map) | set(ncbi_map) | set(local_map))

    rows = []
    for k in all_keys:
        pretty = gbif_map.get(k) or ncbi_map.get(k) or local_map.get(k) or k
        rows.append({
            "name": pretty,
            "key": k,
            "GBIF": k in gbif_map,
            "NCBI": k in ncbi_map,
            "LOCAL": k in local_map,
        })

    # Simple HTML table with search + sort
    html = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8"/>
  <title>{res.get("genus","Genus")} species membership</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 18px; }}
    input {{ padding: 8px; width: 420px; }}
    table {{ border-collapse: collapse; margin-top: 12px; width: 100%; }}
    th, td {{ border: 1px solid #ddd; padding: 8px; }}
    th {{ cursor: pointer; background: #f6f6f6; position: sticky; top: 0; }}
    tr:nth-child(even) {{ background: #fafafa; }}
    .yes {{ font-weight: bold; }}
    .badge-only {{ background: #ffecec; }}
    .badge-all {{ background: #eaffea; }}
  </style>
</head>
<body>
  <h2>{res.get("genus","Genus")} — species membership (GBIF / NCBI / Local)</h2>
  <div>
    <input id="q" placeholder="Search species..." onkeyup="filterRows()"/>
    <span id="count"></span>
  </div>

  <table id="tbl">
    <thead>
      <tr>
        <th onclick="sortTable(0)">Species</th>
        <th onclick="sortTable(1)">GBIF</th>
        <th onclick="sortTable(2)">NCBI</th>
        <th onclick="sortTable(3)">LOCAL</th>
        <th onclick="sortTable(4)">Bucket</th>
      </tr>
    </thead>
    <tbody>
"""
    def bucket_of(r):
        g, n, l = r["GBIF"], r["NCBI"], r["LOCAL"]
        if g and n and l: return "ALL_3"
        if g and n and not l: return "GBIF_NCBI"
        if g and l and not n: return "GBIF_LOCAL"
        if n and l and not g: return "NCBI_LOCAL"
        if g and not n and not l: return "GBIF_ONLY"
        if n and not g and not l: return "NCBI_ONLY"
        if l and not g and not n: return "LOCAL_ONLY"
        return "?"

    for r in rows:
        bucket = bucket_of(r)
        cls = "badge-all" if bucket == "ALL_3" else ("badge-only" if bucket.endswith("_ONLY") else "")
        html += f"""      <tr class="{cls}">
        <td>{r["name"]}</td>
        <td>{"✔" if r["GBIF"] else ""}</td>
        <td>{"✔" if r["NCBI"] else ""}</td>
        <td>{"✔" if r["LOCAL"] else ""}</td>
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



def main(
    taxon: str = "Alca",
    taxon_rank: str = "GENUS",
    acc2tax: str = "/path/to/acc2tax.tsv",
    nodes: str = "/path/to/nodes.dmp",
    names: Optional[str] = None,
    email: str = "your.email@domain.com",
    api_key: Optional[str] = None,
    kingdom: str = "Animalia",
    count_by: str = "versioned",
    pretty: bool = True,
):
    gbif_species_list, gbif_match = gbif_list_species_under_taxon(
        taxon_name=taxon,
        taxon_rank=taxon_rank,
        accepted_only=True,
        prefer_kingdom=kingdom
    )

    res = compare_taxon_gbif_ncbi_local(
        taxon_name=taxon,
        taxon_rank=taxon_rank,
        gbif_species_list=gbif_species_list,
        acc2tax_path=acc2tax,
        nodes_dmp_path=nodes,
        names_dmp_path=names,
        email_for_ncbi=email,
        api_key=api_key,
        count_by=count_by
    )

    if pretty:
        print(json.dumps(res, indent=2, sort_keys=True))
    else:
        print(json.dumps(res))

    return res


if __name__ == "__main__":
    pass

    res = main(
        taxon=taxon,
        taxon_rank=taxon_rank,
        acc2tax="all_acc2tax.with_dbid.gz",
        nodes="/datasets/caeg_dataset/references/ncbi/20250530/taxonomy/ncbi/taxdump/nodes.dmp",
        names="/datasets/caeg_dataset/references/ncbi/20250530/taxonomy/ncbi/taxdump/names.dmp",
        email="n.alexandra.vogel@gmail.com",
        count_by="versioned",
        pretty=True
    )

    buckets = build_species_buckets(res)
    res["species_buckets_counts"] = summarize_buckets(buckets)

    write_bucket_lists(buckets, out_dir=".", prefix=f"{taxon}_{taxon_rank}")
    write_species_membership_html(res, out_html=f"{taxon}_{taxon_rank}_species_membership.html")
    



