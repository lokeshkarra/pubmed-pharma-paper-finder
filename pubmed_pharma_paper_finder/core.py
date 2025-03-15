"""
Core functionality for fetching and filtering PubMed papers.
"""

import re
import time
import logging
import csv
import io
import concurrent.futures
from typing import Dict, List, Optional, Any
from Bio import Entrez, Medline

# Configure logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('pubmed_pharma_paper_finder')

# Regex for identifying pharmaceutical/biotech company affiliations
PHARMA_BIOTECH_TERMS = [
    r'pfizer', r'novartis', r'roche', r'merck', r'johnson\s*&\s*johnson', r'j\s*&\s*j',
    r'sanofi', r'glaxosmithkline', r'gsk', r'abbvie', r'amgen', r'gilead', r'astrazeneca',
    r'eli\s*lilly', r'lilly', r'bristol\s*myers\s*squibb', r'bms', r'boehringer\s*ingelheim',
    r'genentech', r'biogen', r'regeneron', r'moderna', r'biontech', r'vertex', r'illumina',
    r'alexion', r'biomarin', r'seagen', r'incyte', r'alnylam', r'exact\s*sciences',
    r'inc\.', r'llc', r'ltd\.', r'corp\.', r'corporation', r'pharmaceuticals', r'pharma',
    r'therapeutics', r'biotech', r'biosciences', r'biopharmaceuticals', r'biotherapeutics',
    r'laboratories', r'labs'
]
COMBINED_PATTERN = re.compile('|'.join(PHARMA_BIOTECH_TERMS), re.IGNORECASE)

ACADEMIC_TERMS = [
    "university", "college", "school of", "faculty of", "research center",
    "institute of", "department of", "hospital", "medical center", "unit", "laboratory"
]

EMAIL_PATTERN = re.compile(
    r'[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}', re.IGNORECASE
)

class PubMedPaperFinder:
    def __init__(self, email: str, api_key: Optional[str] = None, debug: bool = False):
        self.email = email
        self.api_key = api_key
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        if debug:
            logger.setLevel(logging.DEBUG)
        logger.debug("PubMedPaperFinder initialized with email: %s", email)

    def search_papers(self, query: str, max_results: int = 100) -> List[str]:
        try:
            search_handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, usehistory="y")
            search_results = Entrez.read(search_handle)
            search_handle.close()
            id_list = search_results.get("IdList", [])
            logger.debug(f"Retrieved {len(id_list)} articles from PubMed")
            return id_list
        except Exception as e:
            logger.error(f"Error searching PubMed: {e}")
            return []

    def fetch_paper_details(self, id_list: List[str]) -> List[Dict[str, Any]]:
        if not id_list:
            return []
        
        def fetch_batch(batch, max_retries=5):
            attempt = 0
            while attempt < max_retries:
                try:
                    time.sleep(1 + (0.5 * attempt))  # Exponential backoff
                    fetch_handle = Entrez.efetch(db="pubmed", id=",".join(batch), rettype="medline", retmode="text")
                    papers = list(Medline.parse(fetch_handle))
                    fetch_handle.close()
                    return papers
                except Exception as e:
                    if "429" in str(e):  # Too many requests error
                        wait_time = 2 ** attempt
                        logger.warning(f"Rate limited (429). Retrying in {wait_time} seconds...")
                        time.sleep(wait_time)
                        attempt += 1
                    else:
                        logger.error(f"Error fetching batch details: {e}")
                        break
            return []
        
        all_papers = []
        batch_size = 100  # Fetch 100 articles per batch
        with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
            future_to_batch = {executor.submit(fetch_batch, id_list[i:i+batch_size]): i for i in range(0, len(id_list), batch_size)}
            for future in concurrent.futures.as_completed(future_to_batch):
                all_papers.extend(future.result())

        return all_papers


    def is_affiliated_with_company(self, affiliation: str) -> bool:
        return bool(COMBINED_PATTERN.search(affiliation)) and not self.is_academic_affiliation(affiliation)

    def is_academic_affiliation(self, affiliation: str) -> bool:
        return any(term in affiliation.lower() for term in ACADEMIC_TERMS)

    def extract_email(self, paper: Dict[str, Any]) -> Optional[str]:
        def search_email(data):
            if isinstance(data, str):
                match = EMAIL_PATTERN.search(data)
                return match.group(0) if match else None
            elif isinstance(data, list):
                for item in data:
                    result = search_email(item)
                    if result:
                        return result
            elif isinstance(data, dict):
                for value in data.values():
                    result = search_email(value)
                    if result:
                        return result
            return None
        return search_email(paper)


    def filter_company_affiliated_papers(self, papers: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        filtered_papers = []

        for paper in papers:
            pubmed_id = paper.get("PMID", "")
            title = paper.get("TI", "Unknown Title")
            publication_date = paper.get("DP", "Unknown Date")

            # Extract authors safely
            authors = paper.get("AU", paper.get("FAU", []))  # AU is preferred, FAU as fallback
            authors = authors if isinstance(authors, list) else []
            if not authors:
                logger.warning(f"Missing authors for PubMed ID {pubmed_id}")

            # Extract affiliations safely
            affiliations = paper.get("AD", [])
            affiliations = affiliations if isinstance(affiliations, list) else []

            # Log extracted authors and affiliations for debugging
            # logger.debug(f"PubMed ID {pubmed_id} - Extracted Authors: {authors}")
            # logger.debug(f"PubMed ID {pubmed_id} - Extracted Affiliations: {affiliations}")

            # Map authors to affiliations if possible
            author_affiliation_map = {}
            for i in range(min(len(authors), len(affiliations))):
                author_affiliation_map[authors[i]] = affiliations[i]

            # Identify company-affiliated authors
            company_authors = [
                author for author, aff in author_affiliation_map.items()
                if self.is_affiliated_with_company(aff)
            ]
            company_affiliations = list(set(aff for aff in affiliations if self.is_affiliated_with_company(aff)))

            # Fix: If authors exist but affiliations don't, still include authors
            if company_authors:
                non_academic_authors = ", ".join(company_authors)
            else:
                non_academic_authors = ", ".join(authors) if authors else "Unknown"  # Always include known authors

            # Store paper data
            if company_affiliations:
                filtered_papers.append({
                    "pubmed_id": pubmed_id,
                    "title": title,
                    "publication_date": publication_date,
                    "non_academic_authors": non_academic_authors,  # Fix: Always store authors
                    "company_affiliations": company_affiliations if company_affiliations else ["Unknown"],
                    "corresponding_email": self.extract_email(paper),
                })

        return filtered_papers

    def run_query(self, query: str, max_results: int = 100) -> List[Dict[str, Any]]:
        logger.info(f"Running query: {query}")
        id_list = self.search_papers(query, max_results)
        papers = self.fetch_paper_details(id_list)
        return self.filter_company_affiliated_papers(papers) if papers else []

    def save_to_csv(self, papers: List[Dict[str, Any]], filename: Optional[str] = None) -> Optional[str]:
        if not papers:
            return None
        logger.info(f"Saving {len(papers)} papers to {filename}")
        fieldnames = ['pubmed_id', 'title', 'publication_date', 'non_academic_authors', 'company_affiliations', 'corresponding_email']
        formatted_papers = [{k: '; '.join(v) if isinstance(v, list) else v for k, v in paper.items()} for paper in papers]

        if filename:
            with open(filename, 'w', newline='', encoding='utf-8') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(formatted_papers)
            logger.info(f"CSV file {filename} saved successfully")
        else:
            output = io.StringIO()
            writer = csv.DictWriter(output, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(formatted_papers)
            return output.getvalue()
