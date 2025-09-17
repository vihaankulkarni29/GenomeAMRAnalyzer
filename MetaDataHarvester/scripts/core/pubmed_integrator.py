#!/usr/bin/env python3
"""
PubMedIntegrator - Clinical study linking and context integration
Links AMR mutations with clinical studies and resistance phenotypes

Author: MetaDataHarvester Pipeline
Version: 2.0 - PubMed Integration
"""

import os
import sys
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass
import pandas as pd
import numpy as np
from Bio import Entrez
import requests
import json
import time
import re
from datetime import datetime
from collections import defaultdict


@dataclass
class ClinicalStudy:
    """Clinical study information"""
    pmid: str
    title: str
    abstract: str
    authors: List[str]
    journal: str
    publication_date: str
    doi: Optional[str]
    keywords: List[str]
    study_type: str
    sample_size: Optional[int]
    location: Optional[str]
    resistance_mechanisms: List[str]
    clinical_outcomes: List[str]

    @property
    def relevance_score(self) -> float:
        """Calculate relevance score for AMR research"""
        score = 0.0

        # Keywords related to AMR
        amr_keywords = [
            'antimicrobial resistance', 'drug resistance', 'efflux pump',
            'mutation', 'MIC', 'susceptibility', 'resistance gene',
            'multidrug resistance', 'antibiotic resistance'
        ]

        for keyword in amr_keywords:
            if keyword.lower() in self.title.lower() or keyword.lower() in self.abstract.lower():
                score += 1.0

        # Study type bonus
        if self.study_type in ['Clinical Trial', 'Clinical Study', 'Observational Study']:
            score += 2.0

        # Sample size bonus
        if self.sample_size and self.sample_size > 50:
            score += 1.0

        return min(score, 10.0)  # Cap at 10


@dataclass
class MutationClinicalLink:
    """Link between mutation and clinical study"""
    mutation: str
    study: ClinicalStudy
    evidence_type: str  # 'direct', 'indirect', 'associated'
    confidence: float
    clinical_context: str
    resistance_phenotype: Optional[str]


class PubMedIntegrator:
    """
    Integrates PubMed clinical studies with AMR mutation data
    """

    def __init__(self, email: str, cache_dir: str = "pubmed_cache"):
        self.email = email
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)

        # Setup logging
        self.logger = logging.getLogger(__name__)
        self._setup_logging()

        # Initialize NCBI
        Entrez.email = email

        # Cache for API responses
        self.cache_file = self.cache_dir / "pubmed_cache.json"
        self.study_cache = self._load_cache()

    def _setup_logging(self):
        """Setup logging configuration"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )

    def _load_cache(self) -> Dict:
        """Load cached PubMed data"""
        if self.cache_file.exists():
            try:
                with open(self.cache_file, 'r') as f:
                    return json.load(f)
            except Exception as e:
                self.logger.warning(f"Failed to load cache: {e}")
        return {}

    def _save_cache(self):
        """Save PubMed data to cache"""
        try:
            with open(self.cache_file, 'w') as f:
                json.dump(self.study_cache, f, indent=2)
        except Exception as e:
            self.logger.warning(f"Failed to save cache: {e}")

    def search_clinical_studies(self, query_terms: List[str],
                              max_results: int = 100) -> List[ClinicalStudy]:
        """
        Search PubMed for clinical studies related to AMR

        Args:
            query_terms: List of search terms
            max_results: Maximum number of results to return

        Returns:
            List of clinical studies
        """
        self.logger.info(f"Searching PubMed for clinical studies: {query_terms}")

        # Build PubMed query
        base_terms = [
            'antimicrobial resistance',
            'bacterial infection',
            'clinical study',
            'resistance mechanism'
        ]

        all_terms = base_terms + query_terms
        query = ' AND '.join(f'"{term}"[Title/Abstract]' for term in all_terms)

        try:
            # Search PubMed
            handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
            search_results = Entrez.read(handle)
            handle.close()

            pmids = search_results.get('IdList', [])
            self.logger.info(f"Found {len(pmids)} potential studies")

            if not pmids:
                return []

            # Fetch detailed information
            studies = []
            for pmid in pmids:
                if pmid in self.study_cache:
                    study = self._load_study_from_cache(pmid)
                else:
                    study = self._fetch_study_details(pmid)
                    if study:
                        self.study_cache[pmid] = self._study_to_dict(study)

                if study:
                    studies.append(study)

            # Save cache
            self._save_cache()

            # Filter and rank by relevance
            relevant_studies = [s for s in studies if s.relevance_score > 2.0]
            relevant_studies.sort(key=lambda x: x.relevance_score, reverse=True)

            self.logger.info(f"Found {len(relevant_studies)} relevant clinical studies")
            return relevant_studies[:max_results]

        except Exception as e:
            self.logger.error(f"PubMed search failed: {e}")
            return []

    def _fetch_study_details(self, pmid: str) -> Optional[ClinicalStudy]:
        """Fetch detailed information for a single study"""
        try:
            # Fetch abstract and metadata
            handle = Entrez.efetch(db="pubmed", id=pmid, rettype="abstract", retmode="xml")
            records = Entrez.read(handle)
            handle.close()

            if not records or 'PubmedArticle' not in records[0]:
                return None

            article = records[0]['PubmedArticle']['MedlineCitation']['Article']

            # Extract basic information
            title = article.get('ArticleTitle', '')
            abstract_text = ''

            if 'Abstract' in article:
                abstract_parts = article['Abstract'].get('AbstractText', [])
                if isinstance(abstract_parts, list):
                    abstract_text = ' '.join(str(part) for part in abstract_parts)
                else:
                    abstract_text = str(abstract_parts)

            # Extract authors
            authors = []
            if 'AuthorList' in article:
                for author in article['AuthorList']:
                    if 'LastName' in author and 'ForeName' in author:
                        authors.append(f"{author['ForeName']} {author['LastName']}")

            # Extract journal
            journal = ''
            if 'Journal' in article:
                journal_info = article['Journal']
                if 'Title' in journal_info:
                    journal = journal_info['Title']

            # Extract publication date
            pub_date = ''
            if 'Journal' in article and 'JournalIssue' in article['Journal']:
                issue = article['Journal']['JournalIssue']
                if 'PubDate' in issue:
                    pub_date_info = issue['PubDate']
                    if 'Year' in pub_date_info:
                        pub_date = pub_date_info['Year']

            # Extract keywords
            keywords = []
            if 'KeywordList' in records[0]['PubmedArticle']['MedlineCitation']:
                for keyword_list in records[0]['PubmedArticle']['MedlineCitation']['KeywordList']:
                    if isinstance(keyword_list, list):
                        keywords.extend(keyword_list)

            # Determine study type from title/abstract
            study_type = self._classify_study_type(title, abstract_text)

            # Extract additional metadata
            sample_size = self._extract_sample_size(abstract_text)
            location = self._extract_location(abstract_text)
            resistance_mechanisms = self._extract_resistance_mechanisms(abstract_text)
            clinical_outcomes = self._extract_clinical_outcomes(abstract_text)

            return ClinicalStudy(
                pmid=pmid,
                title=title,
                abstract=abstract_text,
                authors=authors,
                journal=journal,
                publication_date=pub_date,
                doi=None,  # Would need additional API call
                keywords=keywords,
                study_type=study_type,
                sample_size=sample_size,
                location=location,
                resistance_mechanisms=resistance_mechanisms,
                clinical_outcomes=clinical_outcomes
            )

        except Exception as e:
            self.logger.warning(f"Failed to fetch study {pmid}: {e}")
            return None

    def _classify_study_type(self, title: str, abstract: str) -> str:
        """Classify the type of clinical study"""
        text = (title + ' ' + abstract).lower()

        if 'randomized controlled trial' in text or 'rct' in text:
            return 'Randomized Controlled Trial'
        elif 'clinical trial' in text:
            return 'Clinical Trial'
        elif 'cohort study' in text or 'cohort' in text:
            return 'Cohort Study'
        elif 'case control' in text:
            return 'Case-Control Study'
        elif 'observational' in text:
            return 'Observational Study'
        elif 'surveillance' in text or 'epidemiology' in text:
            return 'Surveillance Study'
        elif 'review' in text or 'meta-analysis' in text:
            return 'Review/Meta-analysis'
        else:
            return 'Clinical Study'

    def _extract_sample_size(self, text: str) -> Optional[int]:
        """Extract sample size from text"""
        patterns = [
            r'sample size of (\d+)',
            r'n\s*=\s*(\d+)',
            r'(\d+)\s+patients',
            r'(\d+)\s+isolates',
            r'(\d+)\s+strains'
        ]

        for pattern in patterns:
            matches = re.findall(pattern, text, re.IGNORECASE)
            if matches:
                try:
                    return int(matches[0])
                except ValueError:
                    continue
        return None

    def _extract_location(self, text: str) -> Optional[str]:
        """Extract study location from text"""
        # Common location patterns
        location_patterns = [
            r'in\s+([A-Z][a-z]+(?:\s+[A-Z][a-z]+)*)',
            r'([A-Z][a-z]+(?:\s+[A-Z][a-z]+)*)\s+hospital',
            r'([A-Z][a-z]+(?:\s+[A-Z][a-z]+)*)\s+university'
        ]

        for pattern in location_patterns:
            matches = re.findall(pattern, text)
            if matches:
                return matches[0]
        return None

    def _extract_resistance_mechanisms(self, text: str) -> List[str]:
        """Extract resistance mechanisms mentioned"""
        mechanisms = []

        resistance_terms = [
            'efflux pump', 'efflux', 'mutation', 'gene', 'plasmid',
            'transposon', 'integron', 'resistance gene', 'beta-lactamase',
            'carbapenemase', 'extended-spectrum beta-lactamase'
        ]

        for term in resistance_terms:
            if term.lower() in text.lower():
                mechanisms.append(term)

        return list(set(mechanisms))

    def _extract_clinical_outcomes(self, text: str) -> List[str]:
        """Extract clinical outcomes mentioned"""
        outcomes = []

        outcome_terms = [
            'mortality', 'survival', 'treatment failure', 'cure rate',
            'infection control', 'antibiotic therapy', 'clinical outcome'
        ]

        for term in outcome_terms:
            if term.lower() in text.lower():
                outcomes.append(term)

        return list(set(outcomes))

    def link_mutations_to_studies(self, mutations: List[str],
                                organism: str) -> List[MutationClinicalLink]:
        """
        Link mutations to relevant clinical studies

        Args:
            mutations: List of mutation identifiers
            organism: Bacterial organism name

        Returns:
            List of mutation-clinical study links
        """
        self.logger.info(f"Linking {len(mutations)} mutations to clinical studies for {organism}")

        links = []

        # Search for relevant studies
        search_terms = [organism, 'antimicrobial resistance', 'mutation']
        studies = self.search_clinical_studies(search_terms, max_results=50)

        # For each mutation, find relevant studies
        for mutation in mutations:
            mutation_links = self._find_mutation_studies(mutation, studies)
            links.extend(mutation_links)

        self.logger.info(f"Created {len(links)} mutation-study links")
        return links

    def _find_mutation_studies(self, mutation: str,
                             studies: List[ClinicalStudy]) -> List[MutationClinicalLink]:
        """Find studies relevant to a specific mutation"""
        links = []

        for study in studies:
            # Check if mutation or related mechanism is mentioned
            mutation_mentioned = mutation in study.title or mutation in study.abstract
            mechanism_related = any(mech in study.resistance_mechanisms
                                  for mech in ['efflux pump', 'mutation', 'resistance gene'])

            if mutation_mentioned or mechanism_related:
                evidence_type = 'direct' if mutation_mentioned else 'indirect'
                confidence = 0.8 if mutation_mentioned else 0.5

                # Extract resistance phenotype if available
                resistance_phenotype = self._extract_resistance_phenotype(study.abstract)

                link = MutationClinicalLink(
                    mutation=mutation,
                    study=study,
                    evidence_type=evidence_type,
                    confidence=confidence,
                    clinical_context=study.abstract[:200] + '...' if len(study.abstract) > 200 else study.abstract,
                    resistance_phenotype=resistance_phenotype
                )

                links.append(link)

        return links

    def _extract_resistance_phenotype(self, text: str) -> Optional[str]:
        """Extract resistance phenotype from text"""
        phenotypes = ['MDR', 'XDR', 'resistant', 'susceptible', 'intermediate']

        for phenotype in phenotypes:
            if phenotype.lower() in text.lower():
                return phenotype

        return None

    def _load_study_from_cache(self, pmid: str) -> Optional[ClinicalStudy]:
        """Load study from cache"""
        if pmid not in self.study_cache:
            return None

        data = self.study_cache[pmid]
        try:
            return ClinicalStudy(
                pmid=data['pmid'],
                title=data['title'],
                abstract=data['abstract'],
                authors=data.get('authors', []),
                journal=data.get('journal', ''),
                publication_date=data.get('publication_date', ''),
                doi=data.get('doi'),
                keywords=data.get('keywords', []),
                study_type=data.get('study_type', 'Clinical Study'),
                sample_size=data.get('sample_size'),
                location=data.get('location'),
                resistance_mechanisms=data.get('resistance_mechanisms', []),
                clinical_outcomes=data.get('clinical_outcomes', [])
            )
        except Exception as e:
            self.logger.warning(f"Failed to load cached study {pmid}: {e}")
            return None

    def _study_to_dict(self, study: ClinicalStudy) -> Dict:
        """Convert study to dictionary for caching"""
        return {
            'pmid': study.pmid,
            'title': study.title,
            'abstract': study.abstract,
            'authors': study.authors,
            'journal': study.journal,
            'publication_date': study.publication_date,
            'doi': study.doi,
            'keywords': study.keywords,
            'study_type': study.study_type,
            'sample_size': study.sample_size,
            'location': study.location,
            'resistance_mechanisms': study.resistance_mechanisms,
            'clinical_outcomes': study.clinical_outcomes
        }

    def generate_clinical_report(self, links: List[MutationClinicalLink]) -> None:
        """Generate clinical integration report"""
        report_file = self.cache_dir / "clinical_integration_report.txt"

        with open(report_file, 'w') as f:
            f.write("CLINICAL STUDY INTEGRATION REPORT\n")
            f.write("=" * 40 + "\n\n")

            f.write("SUMMARY:\n")
            f.write("-" * 10 + "\n")
            f.write(f"Total mutation-study links: {len(links)}\n")

            if links:
                mutations = set(link.mutation for link in links)
                studies = set(link.study.pmid for link in links)

                f.write(f"Unique mutations: {len(mutations)}\n")
                f.write(f"Unique studies: {len(studies)}\n\n")

                # Evidence types
                direct_links = [l for l in links if l.evidence_type == 'direct']
                indirect_links = [l for l in links if l.evidence_type == 'indirect']

                f.write("EVIDENCE TYPES:\n")
                f.write("-" * 15 + "\n")
                f.write(f"Direct evidence: {len(direct_links)}\n")
                f.write(f"Indirect evidence: {len(indirect_links)}\n\n")

                # Top mutations by evidence
                mutation_counts = defaultdict(int)
                for link in links:
                    mutation_counts[link.mutation] += 1

                f.write("TOP MUTATIONS BY CLINICAL EVIDENCE:\n")
                f.write("-" * 40 + "\n")
                sorted_mutations = sorted(mutation_counts.items(), key=lambda x: x[1], reverse=True)
                for mutation, count in sorted_mutations[:10]:
                    f.write(f"{mutation}: {count} studies\n")

                f.write("\nDETAILED LINKS:\n")
                f.write("-" * 15 + "\n")
                for i, link in enumerate(links[:20]):  # Show first 20
                    f.write(f"\n{i+1}. {link.mutation}\n")
                    f.write(f"   Study: {link.study.pmid} - {link.study.title[:80]}...\n")
                    f.write(f"   Evidence: {link.evidence_type} (confidence: {link.confidence:.2f})\n")
                    f.write(f"   Context: {link.clinical_context[:100]}...\n")

        self.logger.info(f"Clinical integration report saved: {report_file}")


def main():
    """Main entry point"""
    import argparse

    parser = argparse.ArgumentParser(
        description="PubMed clinical study integration system"
    )

    parser.add_argument(
        "--organism",
        required=True,
        help="Bacterial organism name"
    )

    parser.add_argument(
        "--mutation-file",
        help="File containing mutations to link"
    )

    parser.add_argument(
        "--email",
        required=True,
        help="NCBI email for API access"
    )

    parser.add_argument(
        "--output-dir",
        default="pubmed_data",
        help="Output directory for clinical data"
    )

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # Initialize PubMed integrator
    integrator = PubMedIntegrator(args.email, args.output_dir)

    # Load mutations if provided
    mutations = []
    if args.mutation_file and os.path.exists(args.mutation_file):
        with open(args.mutation_file, 'r') as f:
            mutations = [line.strip() for line in f if line.strip()]

    # Search for clinical studies
    studies = integrator.search_clinical_studies([args.organism])

    # Link mutations to studies if available
    if mutations:
        links = integrator.link_mutations_to_studies(mutations, args.organism)
        integrator.generate_clinical_report(links)
    else:
        # Just generate study report
        with open(Path(args.output_dir) / "clinical_studies.txt", 'w') as f:
            f.write(f"Found {len(studies)} clinical studies for {args.organism}\n\n")
            for study in studies[:10]:  # Show top 10
                f.write(f"PMID: {study.pmid}\n")
                f.write(f"Title: {study.title}\n")
                f.write(f"Relevance Score: {study.relevance_score:.2f}\n")
                f.write(f"Study Type: {study.study_type}\n")
                f.write("-" * 50 + "\n")

    print(f"\nClinical study integration complete!")
    print(f"Found {len(studies)} relevant studies")
    print(f"Results saved to: {args.output_dir}")


if __name__ == "__main__":
    main()