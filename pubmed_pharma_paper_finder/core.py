
"""
Core functionality for fetching and filtering PubMed papers.
"""

import re
import time
import logging
import csv
import io
from typing import Dict, List, Optional, Any

import requests
from Bio import Entrez, Medline

# Configure logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('pubmed_pharma_paper_finder')

# Common pharmaceutical/biotech companies and terms for matching
PHARMA_BIOTECH_TERMS = [
    r'pfizer', r'novartis', r'roche', r'merck', r'johnson\s*&\s*johnson', r'j\s*&\s*j',
    r'sanofi', r'glaxosmithkline', r'gsk', r'abbvie', r'amgen', r'gilead', r'astrazeneca',
    r'eli\s*lilly', r'lilly', r'bristol\s*myers\s*squibb', r'bms', r'boehringer\s*ingelheim',
    
    # Biotech companies
    r'genentech', r'biogen', r'regeneron', r'moderna', r'biontech', r'vertex', r'illumina',
    r'alexion', r'biomarin', r'seagen', r'incyte', r'alnylam', r'exact\s*sciences',
    
    # Generic company identifiers
    r'inc\.', r'llc', r'ltd\.', r'corp\.', r'corporation', r'pharmaceuticals', r'pharma',
    r'therapeutics', r'biotech', r'biosciences', r'biopharmaceuticals', r'biotherapeutics',
    r'laboratories', r'labs'
]
COMBINED_PATTERN = re.compile('|'.join(PHARMA_BIOTECH_TERMS), re.IGNORECASE)

# List of academic/research terms to exclude
ACADEMIC_TERMS = [
    "university", "college", "school of", "faculty of", "research center", 
    "institute of", "department of", "hospital", "medical center"
]

class PubMedPaperFinder:
    """
    Class to fetch and filter PubMed papers based on company affiliations.
    """
    
    def __init__(self, email: str, api_key: Optional[str] = None, debug: bool = False):
        """
        Initialize the PubMedPaperFinder.
        """
        self.email = email
        self.api_key = api_key
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        if debug:
            logger.setLevel(logging.DEBUG)
        logger.debug("PubMedPaperFinder initialized with email: %s", email)
    
    def search_papers(self, query: str, max_results: int = 100) -> List[str]:
        """
        Search PubMed for papers matching the query.
        """
        try:
            search_handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, usehistory="y")
            search_results = Entrez.read(search_handle)
            search_handle.close()
            return search_results.get("IdList", [])
        except Exception as e:
            logger.error(f"Error searching PubMed: {e}")
            return []
    
    def fetch_paper_details(self, id_list: List[str]) -> List[Dict[str, Any]]:
        """
        Fetch detailed information for a list of PubMed IDs.
        """
        if not id_list:
            return []
        
        try:
            batch_size = 100
            all_papers = []
            for i in range(0, len(id_list), batch_size):
                time.sleep(1)  # Avoid hitting API rate limits
                fetch_handle = Entrez.efetch(db="pubmed", id=",".join(id_list[i:i+batch_size]), rettype="medline", retmode="text")
                all_papers.extend(list(Medline.parse(fetch_handle)))
                fetch_handle.close()
            return all_papers
        except Exception as e:
            logger.error(f"Error fetching paper details: {e}")
            return []
    
   

    def is_affiliated_with_company(self, affiliation: str) -> bool:
        """
        Check if an affiliation string belongs to a pharmaceutical/biotech company.
        """
        return bool(COMBINED_PATTERN.search(affiliation)) and not self.is_academic_affiliation(affiliation)

    def is_academic_affiliation(self, affiliation: str) -> bool:
        """
        Check if an affiliation is an academic or research institution.
        """
        ACADEMIC_TERMS = [
            "university", "college", "school of", "faculty of", "research center", 
            "institute of", "department of", "hospital", "medical center", "unit", "laboratory"
        ]
        return any(term in affiliation.lower() for term in ACADEMIC_TERMS)


    

    def extract_company_names(self, affiliation: str) -> List[str]:
        """
        Extract company names from an affiliation string, ensuring correct formatting.
        """
        companies = []
        for part in re.split(r'[,;]', affiliation):
            part = part.strip()
            if self.is_affiliated_with_company(part):
                # Remove extra words like 'Inc.', 'LLC', etc. from the start
                cleaned_part = re.sub(r'^(Inc\.|LLC|Ltd\.|Corp\.)\s*', '', part, flags=re.IGNORECASE)
                companies.append(cleaned_part)

        return companies


    def extract_email_from_paper(self, paper: Dict[str, Any]) -> Optional[str]:
        """
        Extract the corresponding author's email from various fields.
        
        This improved version:
        1. Searches all available fields in the paper
        2. Uses a more comprehensive regex for email detection
        3. Handles both direct string values and nested lists/dictionaries
        """
        # More comprehensive email regex
        email_pattern = r'(?:[a-zA-Z0-9!#$%&\'*+/=?^_`{|}~-]+(?:\.[a-zA-Z0-9!#$%&\'*+/=?^_`{|}~-]+)*)@(?:[a-zA-Z0-9](?:[a-zA-Z0-9-]*[a-zA-Z0-9])?\.)+[a-zA-Z]{2,}'

        
        # Function to recursively search for emails in any data structure
        def search_for_email(data):
            if isinstance(data, str):
                # Search directly in string
                email_match = re.search(email_pattern, data)
                if email_match:
                    return email_match.group(0)
            elif isinstance(data, list):
                # Search in each list item
                for item in data:
                    result = search_for_email(item)
                    if result:
                        return result
            elif isinstance(data, dict):
                # Search in each dictionary value
                for value in data.values():
                    result = search_for_email(value)
                    if result:
                        return result
            return None
        
        # Search all fields in the paper
        for field, value in paper.items():
            email = search_for_email(value)
            if email:
                return email
        
        # If email is present in a specially formatted field like 'corresponding_email'
        if 'corresponding_email' in paper:
            return paper['corresponding_email'] if paper['corresponding_email'] else None
            
        return None
    

    def filter_company_affiliated_papers(self, papers: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Filter papers to include only those with at least one author affiliated with a company.
        """
        filtered_papers = []

        for paper in papers:
            authors = paper.get('AU', [])  # List of author names
            author_affiliations = paper.get('AD', [])  # List of affiliations
            company_authors = []
            company_affiliations = []

            if not authors or not author_affiliations:
                continue  # Skip if no authors or affiliations available

            # Match authors to affiliations using positional mapping
            author_affiliation_map = {}
            for i in range(min(len(authors), len(author_affiliations))):
                author_affiliation_map[authors[i]] = author_affiliations[i]

            # Identify authors affiliated with companies
            for author, affiliation in author_affiliation_map.items():
                if self.is_affiliated_with_company(affiliation):
                    company_authors.append(author)
                    company_affiliations.extend(self.extract_company_names(affiliation))

            # Extract corresponding author email
            email = self.extract_email_from_paper(paper)

            # Only save if company-affiliated authors exist
            if company_affiliations:
                filtered_papers.append({
                    'pubmed_id': paper.get('PMID', ''),
                    'title': paper.get('TI', ''),
                    'publication_date': paper.get('DP', ''),
                    'non_academic_authors': list(set(company_authors)),  # Remove duplicates
                    'company_affiliations': list(set(company_affiliations)),  # Remove duplicates
                    'corresponding_email': email
                })

        logger.info(f"Filtered to {len(filtered_papers)} papers with company affiliations")
        return filtered_papers


    def run_query(self, query: str, max_results: int = 100) -> List[Dict[str, Any]]:
        """
        Run a complete query pipeline: search, fetch details, and filter.
        """
        logger.info(f"Running query: {query}")
        id_list = self.search_papers(query, max_results)
        papers = self.fetch_paper_details(id_list)
        return self.filter_company_affiliated_papers(papers) if papers else []

    def save_to_csv(self, papers: List[Dict[str, Any]], filename: Optional[str] = None) -> Optional[str]:
        """
        Save the filtered papers to a CSV file or return as a string.
        """
        if not papers:
            logger.warning("No papers to save")
            return None

        fieldnames = ['pubmed_id', 'title', 'publication_date', 'non_academic_authors', 'company_affiliations', 'corresponding_email']
        formatted_papers = [{k: '; '.join(v) if isinstance(v, list) else v for k, v in paper.items()} for paper in papers]

        if filename:
            with open(filename, 'w', newline='', encoding='utf-8') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(formatted_papers)
            logger.info(f"Results saved to {filename}")
            return None
        else:
            output = io.StringIO()
            writer = csv.DictWriter(output, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(formatted_papers)
            return output.getvalue()


