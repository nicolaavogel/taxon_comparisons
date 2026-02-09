# Taxonomy Comparison Tool

Compare species lists across **GBIF**, **NCBI**, and **local databases**

This tool compares species/subspecies membership for a given taxon (e.g. genus, family) across:

* **GBIF** (Global Biodiversity Information Facility)
* **NCBI** (taxonomy + nuccore)
* **Local database** (custom accession → taxonomy mapping with `db_id`)

It produces:

* a detailed JSON summary,
* bucket text files (overlaps / uniques),
* an interactive searchable HTML membership table.

## Quick Start

### 1) Quick install (repo only)

```bash
git clone https://github.com/nicolaavogel/taxon_comparisons.git
cd taxon_comparisons
```

*(Environment setup is described in Installation below.)*

### 2) Prepare required input files

You need:

* **Local acc2tax file** (tab-separated; can be `.gz`)

  * Columns: `accession`, `accession.version`, `taxid`, `db_id`
* **NCBI taxonomy dump**

  * `nodes.dmp`
  * `names.dmp` (recommended for readable output)
* Internet access for GBIF + NCBI requests

### 3) Create a JSON config file

Create `run_config.json`:

```json
{
  "taxon": "Poa",
  "taxon_rank": "GENUS",
  "acc2tax": "all_acc2tax.with_dbid.gz",
  "nodes": "/datasets/.../nodes.dmp",
  "names": "/datasets/.../names.dmp",
  "email": "your.email@example.com",
  "api_key": null,
  "kingdom": "Plantae",
  "count_by": "versioned",
  "include_subspecies": true,
  "pretty": true,
  "out_dir": ".",
  "html_out": "Poa_GENUS_membership_by_dbid.html"
}
```

### 4) Run

```bash
python taxon_compare.py --config run_config.json
```

Outputs will be written to the current directory (or `out_dir`).


## Installation

### Option A: Conda (recommended)

```bash
conda env create -f taxon_compare.yml
conda activate taxon_compare
```

### Option B: Existing Python environment

Ensure Python 3.9+ and install the packages listed in `taxon_compare.yml`.

## Configuration (JSON)

The tool is driven by a JSON config file. Keys map to `main()` arguments plus a couple of output helpers.

### Required keys

| Key          | Description                                       |
| ------------ | ------------------------------------------------- |
| `taxon`      | Taxon name (e.g. `"Poa"`, `"Alca"`)               |
| `taxon_rank` | Rank (e.g. `"GENUS"`, `"FAMILY"`)                 |
| `acc2tax`    | Path to acc2tax TSV/TSV.GZ (must include `db_id`) |
| `nodes`      | Path to NCBI `nodes.dmp`                          |
| `email`      | Email required by NCBI Entrez                     |
| `kingdom`    | GBIF kingdom (e.g. "Plantae", "Animalia")         |

### Optional keys

| Key                  |       Default | Description                                         |
| -------------------- | ------------: | --------------------------------------------------- |
| `names`              |        `null` | Path to NCBI `names.dmp` (recommended)              |
| `api_key`            |        `null` | NCBI API key (higher rate limits)                   |
| `count_by`           | `"versioned"` | `"versioned"` or `"unversioned"` accession counting |
| `include_subspecies` |       `false` | Include subspecies in comparisons and outputs       |
| `pretty`             |        `true` | Pretty JSON output to stdout                        |
| `out_dir`            |         `"."` | Where bucket text files are written                 |
| `html_out`           |          auto | HTML output filename                                |


## Input format details

### Local acc2tax file (`acc2tax.tsv` or `.gz`)

Must be **tab-separated** with 4 columns:

```
accession    accession.version    taxid    db_id
XM_032804418 XM_032804418.1       106734   core_nt
KR677567     KR677567.1           288806   core_nt
...
```

* `db_id` is used to split the local dataset into per-database membership columns in the HTML.

### NCBI taxonomy dump

* `nodes.dmp` and `names.dmp` from NCBI Taxonomy dump.

## Outputs

After a run, you get:

### 1) JSON summary (stdout)

Includes:

* counts and name lists from GBIF / NCBI / Local
* per-db local breakdown
* overlap/unique “bucket” counts

### 2) Bucket files (text)

Written under `out_dir`:

* one set for overall local union
* one set per local `db_id`

These correspond to buckets like:

* `ALL_3`
* `GBIF_ONLY`
* `NCBI_ONLY`
* `LOCAL_ONLY`
* `GBIF_NCBI`, `GBIF_LOCAL`, `NCBI_LOCAL`

### 3) HTML membership table

A searchable, sortable table showing:

* GBIF / NCBI presence
* per-local-DB presence columns
* bucket label per species


## Notice

GBIF and NCBI taxonomy are developed independently and do not correspond perfectly. Species and subspecies names may differ due to synonymy, accepted-name choices, rank treatment (species vs. subspecies), and ongoing taxonomic revisions. Differences observed in the comparison output may therefore reflect taxonomy decisions rather than true absence from a database.

This tool currently compares normalized name strings and does not attempt to resolve synonyms or reconcile differing taxonomic concepts across sources.

## Troubleshooting

* **NCBI rate limits:** add `api_key` in your config.
* **“No NCBI taxid…”:** check taxon spelling + rank.
* **Empty local columns:** verify `acc2tax` includes `db_id` and correct taxids.
* **Missing/incorrect paths:** update `nodes`, `names`, `acc2tax` in config.

## Contact

Nicola Vogel

---