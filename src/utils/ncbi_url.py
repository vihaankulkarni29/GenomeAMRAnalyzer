"""
NCBI URL to Accessions Utility

Resolves NCBI nuccore URLs to accession lists using Bio.Entrez.
Supports complex queries with organism names, gene names, and filters.
"""

from __future__ import annotations
import re
import sys
from typing import List, Optional
from urllib.parse import unquote, urlparse

try:
    from Bio import Entrez
    ENTREZ_AVAILABLE = True
except ImportError:
    ENTREZ_AVAILABLE = False
    Entrez = None

# Valid NCBI hosts for nuccore URLs
NCBI_NUCCORE_HOSTS = {"www.ncbi.nlm.nih.gov", "ncbi.nlm.nih.gov"}


def _extract_query_from_url(url: str) -> str:
    """Extract the search query from an NCBI nuccore URL."""
    parsed = urlparse(url)
    
    # Handle different URL formats:
    # https://www.ncbi.nlm.nih.gov/nuccore/?term=Escherichia+coli...
    # https://www.ncbi.nlm.nih.gov/nuccore/(Escherichia%20coli...)
    
    query = ""
    if "term=" in parsed.query:
        # Extract from query parameter
        for part in parsed.query.split("&"):
            if part.startswith("term="):
                query = part[5:]  # Remove "term="
                break
    elif "/nuccore/" in parsed.path:
        # Extract from path
        query = parsed.path.split("/nuccore/")[-1]
    
    # URL decode and clean
    query = unquote(query or "").strip()
    
    # Strip surrounding parentheses if present
    if query.startswith("(") and query.endswith(")"):
        query = query[1:-1]
    
    return query


def url_to_accessions(
    url: str, 
    email: str, 
    api_key: Optional[str] = None, 
    max_ids: int = 100,
    db: str = "nuccore"
) -> List[str]:
    """
    Convert an NCBI URL to a list of accession numbers.
    
    Args:
        url: NCBI nuccore URL with search query
        email: Email for NCBI Entrez (required)
        api_key: Optional NCBI API key for higher rate limits
        max_ids: Maximum number of accessions to retrieve
        db: NCBI database to search (default: nuccore)
    
    Returns:
        List of accession numbers (accession.version format)
    
    Raises:
        ImportError: If Bio.Entrez is not available
        ValueError: If URL is invalid or query cannot be extracted
        RuntimeError: If NCBI query fails
    """
    if not ENTREZ_AVAILABLE:
        raise ImportError(
            "Bio.Entrez is required for URL processing. "
            "Install with: pip install biopython"
        )
    
    if not url or urlparse(url).hostname not in NCBI_NUCCORE_HOSTS:
        raise ValueError(
            f"Provide a valid NCBI nuccore URL from: {', '.join(NCBI_NUCCORE_HOSTS)}"
        )
    
    query = _extract_query_from_url(url)
    if not query:
        raise ValueError("Could not parse search query from URL")
    
    if not email or email.startswith("your.email"):
        raise ValueError("Valid email address is required for NCBI access")
    
    # Configure Entrez
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    
    try:
        # Search for matching records
        print(f"[INFO] Searching NCBI {db} with query: {query}")
        handle = Entrez.esearch(db=db, term=query, retmax=max_ids)
        search_result = Entrez.read(handle)
        handle.close()
        
        ids = search_result.get("IdList", [])
        if not ids:
            print(f"[WARNING] No results found for query: {query}")
            return []
        
        print(f"[INFO] Found {len(ids)} records, fetching accessions...")
        
        # Fetch accession information for each ID
        handle = Entrez.esummary(db=db, id=",".join(ids))
        summaries = Entrez.read(handle)
        handle.close()
        
        accessions = []
        for summary in summaries:
            # Try different fields for accession
            acc = (
                summary.get("AccessionVersion") or 
                summary.get("Caption") or 
                summary.get("Title", "").split()[0]
            )
            if acc and re.match(r"^[A-Z]{1,2}_?\d+\.\d+$", acc):
                accessions.append(acc)
        
        # Remove duplicates and sort
        unique_accessions = sorted(set(accessions))
        print(f"[INFO] Retrieved {len(unique_accessions)} unique accessions")
        
        return unique_accessions
        
    except Exception as e:
        raise RuntimeError(f"NCBI query failed: {e}")


def validate_accessions(accessions: List[str]) -> List[str]:
    """
    Validate and filter accession numbers.
    
    Args:
        accessions: List of potential accession numbers
    
    Returns:
        List of valid accession numbers
    """
    valid_accessions = []
    accession_pattern = re.compile(r"^[A-Z]{1,2}_?\d+(\.\d+)?$")
    
    for acc in accessions:
        acc = acc.strip()
        if accession_pattern.match(acc):
            valid_accessions.append(acc)
        else:
            print(f"[WARNING] Invalid accession format: {acc}")
    
    return valid_accessions


def main():
    """Command line interface for URL to accessions conversion."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Convert NCBI URL to accession list"
    )
    parser.add_argument("url", help="NCBI nuccore URL")
    parser.add_argument("--email", required=True, help="Email for NCBI access")
    parser.add_argument("--api-key", help="NCBI API key (optional)")
    parser.add_argument("--max-ids", type=int, default=100, help="Maximum results")
    parser.add_argument("--output", "-o", help="Output file (default: stdout)")
    
    args = parser.parse_args()
    
    try:
        accessions = url_to_accessions(
            args.url, 
            args.email, 
            args.api_key, 
            args.max_ids
        )
        
        if args.output:
            with open(args.output, "w") as f:
                for acc in accessions:
                    f.write(f"{acc}\n")
            print(f"[INFO] Wrote {len(accessions)} accessions to {args.output}")
        else:
            for acc in accessions:
                print(acc)
                
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()