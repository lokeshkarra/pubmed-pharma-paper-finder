"""
Command-line interface for the PubMed Pharma Paper Finder.
"""

import sys
import argparse
import logging

from .core import PubMedPaperFinder

def main():
    """
    Main entry point for the command-line interface.
    """
    parser = argparse.ArgumentParser(
        description="Fetch PubMed papers with authors affiliated to pharmaceutical or biotech companies."
    )
    
    parser.add_argument(
        "query",
        help="PubMed search query (supports full PubMed query syntax)"
    )
    
    parser.add_argument(
        "-f", "--file",
        help="Output file name for CSV results (default: print to console)",
        default=None
    )
    
    parser.add_argument(
        "-d", "--debug",
        help="Enable debug logging",
        action="store_true"
    )
    
    parser.add_argument(
        "-m", "--max-results",
        help="Maximum number of results to fetch (default: 100)",
        type=int,
        default=100
    )
    
    parser.add_argument(
        "-e", "--email",
        help="Email address to use for NCBI API access (required)",
        required=True
    )
    
    parser.add_argument(
        "-k", "--api-key",
        help="NCBI API key for higher rate limits (optional)",
        default=None
    )
    
    args = parser.parse_args()
    
    # Set up the PubMedPaperFinder
    finder = PubMedPaperFinder(
        email=args.email,
        api_key=args.api_key,
        debug=args.debug
    )
    
    # Run the query
    papers = finder.run_query(args.query, max_results=args.max_results)
    
    if not papers:
        print("No papers found matching the criteria.")
        return 1
    
    # Save or print results
    result = finder.save_to_csv(papers, args.file)
    
    if result and not args.file:
        print(result)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
