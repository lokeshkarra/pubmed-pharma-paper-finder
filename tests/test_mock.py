import pytest
from pytest_mock import mocker
from pubmed_pharma_paper_finder.core import PubMedPaperFinder  # Updated to correct module path
from Bio import Entrez
#import json

@pytest.fixture
def pubmed_finder():
    return PubMedPaperFinder(email="techhonesty11@gmail.com", api_key=None, debug=True)

def mock_read_response(id_list):
    """Return a valid XML response structure expected by Entrez.read()"""
    return f"""<?xml version="1.0" encoding="UTF-8"?>
    <eSearchResult>
        <IdList>
            {''.join(f'<Id>{id}</Id>' for id in id_list)}
        </IdList>
    </eSearchResult>"""


def test_search_papers_valid(mocker, pubmed_finder):
    mock_esearch = mocker.patch("Bio.Entrez.esearch")

    # Ensure the return value is a valid mock object
    mock_response = mocker.MagicMock()
    mock_esearch.return_value = mock_response

    # Mocking `Entrez.read` to return structured data
    mock_entrez_read = mocker.patch("Bio.Entrez.read", return_value={"IdList": ["12345", "67890"]})

    ids = pubmed_finder.search_papers("cancer", max_results=2)

    print("Received search response:", ids)  # Debugging statement
    assert ids == ["12345", "67890"]
    
    mock_esearch.assert_called_once()
    mock_entrez_read.assert_called_once_with(mock_response)



def test_search_papers_no_results(mocker, pubmed_finder):
    mock_esearch = mocker.patch("Bio.Entrez.esearch")
    
    # Properly mock empty response
    mock_response = mocker.MagicMock()
    mock_response.read.return_value = mock_read_response([])
    mock_esearch.return_value = mock_response

    ids = pubmed_finder.search_papers("nonexistentquery", max_results=2)
    assert ids == []

def test_fetch_paper_details_valid(mocker, pubmed_finder):
    mock_efetch = mocker.patch("Bio.Entrez.efetch")

    # Mock response to simulate a valid fetch
    mock_response = mocker.MagicMock()
    mock_response.read.return_value = """PMID- 12345\nTI  - Sample Title\nDP  - 2024\nAU  - John Doe\nAD  - Pfizer Inc., USA"""
    mock_efetch.return_value = mock_response

    mocker.patch("Bio.Medline.parse", return_value=[
        {"PMID": "12345", "TI": "Sample Title", "DP": "2024", "AU": ["John Doe"], "AD": ["Pfizer Inc., USA"]}
    ])

    papers = pubmed_finder.fetch_paper_details(["12345"])
    
    assert len(papers) == 1
    assert papers[0]["PMID"] == "12345"
    assert papers[0]["TI"] == "Sample Title"
    mock_efetch.assert_called_once()

def test_fetch_paper_details_network_error(mocker, pubmed_finder):
    mock_efetch = mocker.patch("Bio.Entrez.efetch", side_effect=Exception("Network Error"))
    papers = pubmed_finder.fetch_paper_details(["12345"])
    assert papers == []

def test_run_query(mocker, pubmed_finder):
    # Fixed incorrect module path
    mock_search = mocker.patch("pubmed_pharma_paper_finder.core.PubMedPaperFinder.search_papers", return_value=["12345"])
    mock_fetch = mocker.patch("pubmed_pharma_paper_finder.core.PubMedPaperFinder.fetch_paper_details", return_value=[
        {"PMID": "12345", "TI": "Sample Title", "DP": "2024", "AU": ["John Doe"], "AD": ["Pfizer Inc., USA"]}
    ])

    papers = pubmed_finder.run_query("cancer")
    
    assert len(papers) == 1
    assert papers[0]["pubmed_id"] == "12345"
    assert papers[0]["title"] == "Sample Title"

    mock_search.assert_called_once()
    mock_fetch.assert_called_once()

def test_save_to_csv(mocker, pubmed_finder):
    mock_open = mocker.patch("builtins.open", new_callable=mocker.mock_open)
    
    papers = [
        {"pubmed_id": "12345", "title": "Sample", "publication_date": "2024", "non_academic_authors": "John Doe", 
         "company_affiliations": ["Pfizer"], "corresponding_email": "jdoe@pfizer.com"}
    ]
    
    pubmed_finder.save_to_csv(papers, "output.csv")

    mock_open.assert_called_once_with("output.csv", "w", newline="", encoding="utf-8")
    handle = mock_open()
    handle.write.assert_called()
