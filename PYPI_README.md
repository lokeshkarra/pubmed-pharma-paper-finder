# PubMed Pharma Paper Finder

**PubMed Pharma Paper Finder** is a Python package that helps retrieve research papers from PubMed and filters results to include only those with authors affiliated with pharmaceutical or biotech companies.

## 🚀 Features

- 🔍 Fetches research papers from **PubMed** using the official API.
- 🏢 Filters results to **include only** papers with non-academic authors affiliated with pharmaceutical/biotech companies.
- 📄 Extracts key details, including:
  - **PubMed ID**
  - **Title**
  - **Publication Date**
  - **Company Affiliations**
  - **Corresponding Author Email** (when available)
- 📊 Outputs results as a CSV file or prints to the console.

---

## 📦 Installation

```bash
pip install biopython==1.85
pip install -i https://test.pypi.org/simple/ pubmed-pharma-paper-finder
```

‼️ Note: Requires Python 3.12+ and biopython>=1.85 ‼️

## 🛠️ Usage

Command-Line Interface

```bash
get-papers-list -e your.email@example.com "cancer AND therapy"
```

Options:

    -e, --email (Required): Email address for NCBI API access.
    -f, --file: Save results to a CSV file.
    -m, --max-results: Set a max number of papers to fetch (default: 100).
    -d, --debug: Enable debug mode.
    -h, --help: Show help.

## Example:

```bash
get-papers-list -e your.email@example.com -f results.csv "diabetes AND treatment"
```

🏗️ Using as a Library

```py
from pubmed_pharma_paper_finder import PubMedPaperFinder

finder = PubMedPaperFinder(email="your.email@example.com")
papers = finder.run_query("cancer AND therapy")

for paper in papers:
    print(f"Title: {paper['title']}")
    print(f"Company Affiliations: {paper['company_affiliations']}")

```

## 📜 License

This project is licensed under the Apache License 2.0.

For full documentation and source code, visit [GitHub](https://github.com/lokeshkarra/pubmed-pharma-paper-finder/tree/main?tab=Apache-2.0-1-ov-file).
