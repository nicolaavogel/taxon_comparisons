```markdown
# Taxonomy Comparison Tool

**Compare species lists across GBIF, NCBI, and local databases**

This Python tool allows you to compare species lists for a given taxon (e.g., genus, family) across three sources:
- **GBIF** (Global Biodiversity Information Facility)
- **NCBI** (National Center for Biotechnology Information)
- **Local database** (custom accession-to-taxonomy mapping)

The tool generates detailed reports, including overlaps, unique entries, and interactive HTML summaries.

---

## Table of Contents
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Configuration](#configuration)
- [Outputs](#outputs)
- [Dependencies](#dependencies)


---

## Features
- Fetch species lists from GBIF and NCBI for any taxon/rank.
- Compare species names across all three sources.
- Generate HTML reports with searchable/sortable tables.
- Save species buckets (overlaps, unique entries) as text files.
- Support for both versioned and unversioned accession counts.

---

## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/your-username/your-repo-name.git
cd your-repo-name
```

### 2. Set Up Conda Environment
A `taxon_compare.yml` file is provided for easy setup:
```bash
conda env create -f taxon_compare.yml
conda activate taxon_compare
```

### 3. Install Dependencies
If not using conda please ensure python3 and the python packages needed in this script are installed.

---

## Usage

### Command Line
Run the script directly:
```bash
python taxon_compare.py
```
By default, it will use the parameters defined in the `main()` function.

### Custom Parameters
You can override defaults by passing arguments:
```python
res = main(
    taxon="Alca",
    taxon_rank="GENUS",
    acc2tax="path/to/acc2tax.tsv",
    nodes="path/to/nodes.dmp",
    names="path/to/names.dmp",
    email="your.email@example.com",
    api_key="your-ncbi-api-key",
    count_by="versioned",
    pretty=True
)
```

### Required Files
- `acc2tax.tsv`: Accession to taxonomy file for your local database. Must include 4 columns `accession, accession.version, taxid, db_id` and be tab seperated. 
	e.g. ```
	accession	accession.version	taxid	db_id
	XM_032804418	XM_032804418.1	106734	core_nt
	KR677567	KR677567.1	288806	core_nt
	KM651970	KM651970.1	215467	core_nt
	DQ848817	DQ848817.1	152344	core_nt
	MT408808	MT408808.1	4932	core_nt
	GU166420	GU166420.1	721828	core_nt
	XM_007924952	XM_007924952.1	383855	core_nt
	XM_009880046	XM_009880046.1	50402	core_nt
	MN403818	MN403818.1	2651588	core_nt
	```
- `nodes.dmp` and `names.dmp`: NCBI taxonomy dump files. 

---

## Configuration

### Environment Variables
- `NCBI_EMAIL`: Your email for NCBI API calls (required).
- `NCBI_API_KEY`: (Optional) Your NCBI API key for higher rate limits.

### Conda Environment
The `environment.yml` file includes all required dependencies:
```yaml
name: taxonomy-comparison
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.9
  - pygbif
  - biopython
  - pip
  - pip:
    - requests
```

---

## Outputs

### 1. JSON Report
A detailed JSON file with:
- Species counts and names from GBIF, NCBI, and local DB.
- Overlap statistics and unique entries.

### 2. Species Buckets
Text files for each bucket (e.g., `Alca_GENUS_ALL_3.txt`, `Alca_GENUS_GBIF_ONLY.txt`).

### 3. HTML Report
An interactive HTML table (`Alca_GENUS_species_membership.html`) with:
- Search and sort functionality.
- Visual indicators for overlaps.

---

## Dependencies
- Python 3.9+
- [pygbif](https://github.com/gbif/pygbif)
- [Biopython](https://biopython.org/)
- [Requests](https://docs.python-requests.org/)

---

## License
This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

## Example Workflow
1. **Fetch Data**: Run the script to fetch species lists from GBIF and NCBI, and compare with your local DB.
2. **Analyze Outputs**: Use the HTML report to explore overlaps and unique entries.
3. **Export Buckets**: Use the text files for further analysis or reporting.

---

## Troubleshooting
- **API Limits**: If you hit NCBI rate limits, add your API key.
- **Missing Files**: Ensure all required files (`nodes.dmp`, `names.dmp`, `acc2tax.tsv`) are in the correct paths.
- **Email Requirement**: NCBI requires a valid email for API access.

---

## Contributing
Pull requests are welcome! For major changes, please open an issue first.

---

## Contact
For questions or feedback, contact [Nicola Vogel](mailto:nicola.vogel@sund.ku.dk).
```
