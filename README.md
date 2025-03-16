# PubMed Pharma Paper Finder

A Python command-line tool that fetches research papers from PubMed based on a user query and filters results to include only papers with at least one author affiliated with a pharmaceutical or biotech company.

## Features

- Fetches papers from PubMed using the official API
- Supports PubMed's full query syntax
- Filters papers to include only those with pharmaceutical/biotech company affiliations
- Extracts key details including:
  - PubMed ID
  - Title
  - Publication Date
  - Non-academic Authors (affiliated with companies)
  - Company Affiliations
  - Corresponding Author Email (when available)
- Outputs results as CSV (to file or console)

## Code Organization

The project follows a modular structure:

```
pubmed-pharma-paper-finder/
├── pubmed_pharma_paper_finder/
│   ├── __init__.py        # Package initialization
│   ├── core.py            # Core functionality module
│   └── cli.py             # Command-line interface
├── tests/                 # Test suite
│   └── test_mock.py       # Unit and integration tests
├── pyproject.toml         # Poetry configuration
├── README.md              # Documentation
└── LICENSE                # License information
```

## Installation

### Prerequisites

- Python 3.8 or higher
- [Poetry](https://python-poetry.org/docs/#installation) for dependency management

### Install from GitHub

```bash
# Clone the repository
git clone https://github.com/lokeshkarra/pubmed-pharma-paper-finder.git
cd pubmed-pharma-paper-finder

# Install with Poetry
poetry install
```

### Install from Test PyPI

```bash
pip install biopython==1.85
pip install -i https://test.pypi.org/simple/ pubmed-pharma-paper-finder
```

‼️ Note: Requires Python 3.12+ and biopython>=1.85 ‼️

## Usage

### Command-Line Interface

```bash
# Basic usage
poetry run get-papers-list -e your.email@example.com "cancer AND therapy"

# Save results to a file
poetry run get-papers-list -e your.email@example.com -f results.csv "cancer AND therapy"

# Enable debug mode
poetry run get-papers-list -e your.email@example.com -d "cancer AND therapy"

# Set maximum number of results
poetry run get-papers-list -e your.email@example.com -m 200 "cancer AND therapy"

# Show help
poetry run get-papers-list -h
```

### Options

- `query`: PubMed search query (supports full PubMed query syntax)
- `-h, --help`: Show usage instructions
- `-d, --debug`: Print debug information
- `-f, --file`: Specify a filename for saving the results (default: print to console)
- `-m, --max-results`: Maximum number of results to fetch (default: 100)
- `-e, --email`: Email address to use for NCBI API access (required)
- `-k, --api-key`: NCBI API key for higher rate limits (optional)

### Using as a Library

You can also use the package as a library in your Python code:

```python
from pubmed_pharma_paper_finder import PubMedPaperFinder

# Create a finder instance
finder = PubMedPaperFinder(email="your.email@example.com")

# Run a query
papers = finder.run_query("cancer AND therapy")

# Process the results
for paper in papers:
    print(f"Title: {paper['title']}")
    print(f"Company Affiliations: {paper['company_affiliations']}")
```

## External Dependencies

- [Biopython](https://biopython.org/): For interacting with NCBI APIs
- [Poetry](https://python-poetry.org/): For dependency management and packaging

## Testing

```bash
poetry run pytest
```

## License

This project is licensed under the Apache License - see the LICENSE file for details.

## Development

Contributions are welcome! Please feel free to submit a Pull Request.

### Development Setup

```bash
# Clone the repository
git clone https://github.com/yourusername/pubmed-pharma-paper-finder.git
cd pubmed-pharma-paper-finder

# Install development dependencies
poetry install


```
