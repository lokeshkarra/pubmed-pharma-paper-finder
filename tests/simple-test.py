"""
Tests for the core functionality of the PubMed Pharma Paper Finder.
"""

import unittest
from unittest.mock import patch, MagicMock

from pubmed_pharma_paper_finder.core import PubMedPaperFinder


class TestPubMedPaperFinder(unittest.TestCase):
    def setUp(self):
        self.finder = PubMedPaperFinder(email="test@example.com")
    
    def test_is_affiliated_with_company(self):
        # Test positive cases
        self.assertTrue(self.finder.is_affiliated_with_company("Pfizer Inc., New York, NY, USA"))
        self.assertTrue(self.finder.is_affiliated_with_company("Department of Research, Novartis Pharmaceuticals, Basel, Switzerland"))
        self.assertTrue(self.finder.is_affiliated_with_company("Johnson & Johnson Pharmaceutical R&D, San Diego, CA"))
        self.assertTrue(self.finder.is_affiliated_with_company("BioMarin Therapeutics LLC"))
        
        # Test negative cases
        self.assertFalse(self.finder.is_affiliated_with_company("Harvard Medical School, Boston, MA"))
        self.assertFalse(self.finder.is_affiliated_with_company("Department of Biology, Stanford University"))
        self.assertFalse(self.finder.is_affiliated_with_company("National Institutes of Health, Bethesda, MD"))
    
    def test_extract_company_names(self):
        # Test extraction from affiliation string
        affiliations = "John Smith, Pfizer Inc., New York, NY, USA"
        companies = self.finder.extract_company_names(affiliations)
        self.assertIn("Pfizer Inc.", companies[0])
        
        # Test with multiple companies
        affiliations = "Jane Doe, Novartis Pharmaceuticals, Basel, Switzerland; Previously at Genentech Inc."
        companies = self.finder.extract_company_names(affiliations)
        self.assertEqual(len(companies), 2)
    
    @patch('pubmed_pharma_paper_finder.core.Entrez')
    def test_search_papers(self, mock_entrez):
        # Mock the search response
        mock_search_handle = MagicMock()
        mock_search_results = {"Count": "10", "IdList": ["12345", "67890"], "WebEnv": "web_env", "QueryKey": "query_key"}
        mock_entrez.esearch.return_value = mock_search_handle
        mock_entrez.read.return_value = mock_search_results
        
        # Call the method
        result = self.finder.search_papers("cancer AND therapy")
        
        # Check the result
        self.assertEqual(result, ["12345", "67890"])
        
        # Verify Entrez was called correctly
        mock_entrez.esearch.assert_called()
    
    @patch('pubmed_pharma_paper_finder.core.Entrez')
    @patch('pubmed_pharma_paper_finder.core.Medline')
    def test_fetch_paper_details(self, mock_medline, mock_entrez):
        # Mock the fetch response
        mock_fetch_handle = MagicMock()
        mock_entrez.efetch.return_value = mock_fetch_handle
        
        # Mock the Medline parse result
        mock_papers = [
            {"PMID": "12345", "TI": "Test Paper 1", "DP": "2023"},
            {"PMID": "67890", "TI": "Test Paper 2", "DP": "2022"}
        ]
        mock_medline.parse.return_value = mock_papers
        
        # Call the method
        result = self.finder.fetch_paper_details(["12345", "67890"])
        
        # Check that Entrez and Medline were called correctly
        mock_entrez.efetch.assert_called_once()
        mock_medline.parse.assert_called_once()
    
    def test_filter_company_affiliated_papers(self):
        # Test papers with and without company affiliations
        papers = [
            {
                "PMID": "12345",
                "TI": "Test Paper with Company",
                "DP": "2023",
                "AU": ["Smith, J", "Doe, J"],
                "AD": "Smith, J, Harvard Medical School, Boston, MA\nDoe, J, Pfizer Inc., New York, NY, USA"
            },
            {
                "PMID": "67890",
                "TI": "Test Paper without Company",
                "DP": "2022",
                "AU": ["Johnson, A", "Williams, B"],
                "AD": "Johnson, A, Stanford University, CA\nWilliams, B, NIH, Bethesda, MD"
            }
        ]
        
        # Call the method
        result = self.finder.filter_company_affiliated_papers(papers)
        
        # Should only include the first paper
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]["pubmed_id"], "12345")
        self.assertIn("Doe, J", result[0]["non_academic_authors"])


if __name__ == "__main__":
    unittest.main()
