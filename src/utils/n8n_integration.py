"""
n8n Integration Manager for GenomeAMRAnalyzer
Handles workflow automation, database integration, and external notifications
"""

import json
import requests
import logging
from datetime import datetime
from typing import Dict, List, Optional, Any
from pathlib import Path
import uuid

class N8nIntegrationManager:
    """
    Manages integration between GenomeAMRAnalyzer pipeline and n8n workflows
    """
    
    def __init__(self, config_file: str = "config/n8n_config.yaml"):
        """Initialize n8n integration manager"""
        self.config = self.load_config(config_file)
        self.webhook_urls = self.config.get('webhook_urls', {})
        self.database_config = self.config.get('database', {})
        self.notification_config = self.config.get('notifications', {})
        
        # Setup logging
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)
    
    def load_config(self, config_file: str) -> Dict:
        """Load n8n integration configuration"""
        config_path = Path(config_file)
        if config_path.exists():
            import yaml
            with open(config_path, 'r') as f:
                return yaml.safe_load(f)
        else:
            # Return default configuration
            return {
                'webhook_urls': {
                    'analysis_complete': 'http://localhost:5678/webhook/amr-analysis-complete',
                    'high_priority_alert': 'http://localhost:5678/webhook/amr-high-priority',
                    'quality_check': 'http://localhost:5678/webhook/amr-quality-check'
                },
                'database': {
                    'enabled': True,
                    'store_results': True,
                    'store_metadata': True
                },
                'notifications': {
                    'email_enabled': True,
                    'slack_enabled': False,
                    'high_priority_threshold': 5
                }
            }
    
    def prepare_pipeline_results(self, pipeline_output: Dict) -> Dict:
        """
        Standardize pipeline results for n8n workflow consumption
        """
        analysis_id = str(uuid.uuid4())
        timestamp = datetime.now().isoformat()
        
        standardized_results = {
            'metadata': {
                'analysis_id': analysis_id,
                'timestamp': timestamp,
                'pipeline_version': pipeline_output.get('version', '1.0.0'),
                'user_email': pipeline_output.get('user_email', 'unknown@example.com'),
                'analysis_type': pipeline_output.get('analysis_type', 'standard')
            },
            'summary': {
                'total_genomes': pipeline_output.get('genome_count', 0),
                'total_resistance_genes': pipeline_output.get('resistance_gene_count', 0),
                'high_priority_mutations': pipeline_output.get('high_priority_count', 0),
                'processing_time_minutes': pipeline_output.get('processing_time', 0),
                'quality_score': pipeline_output.get('quality_score', 0.0)
            },
            'files': {
                'html_report': pipeline_output.get('html_report_path', ''),
                'csv_summary': pipeline_output.get('csv_summary_path', ''),
                'json_detailed': pipeline_output.get('json_results_path', ''),
                'cooccurrence_analysis': pipeline_output.get('cooccurrence_path', ''),
                'alignments_dir': pipeline_output.get('alignments_dir', '')
            },
            'findings': {
                'resistance_genes': pipeline_output.get('resistance_genes', []),
                'mutations': pipeline_output.get('mutations', []),
                'cooccurrence_patterns': pipeline_output.get('cooccurrence_patterns', [])
            },
            'quality_metrics': {
                'genome_completeness': pipeline_output.get('genome_completeness', {}),
                'alignment_quality': pipeline_output.get('alignment_quality', {}),
                'database_coverage': pipeline_output.get('database_coverage', {})
            },
            'parameters': {
                'genes_file': pipeline_output.get('genes_file', ''),
                'genome_source': pipeline_output.get('genome_source', ''),
                'analysis_parameters': pipeline_output.get('analysis_parameters', {})
            }
        }
        
        return standardized_results
    
    def trigger_analysis_complete_workflow(self, pipeline_results: Dict) -> bool:
        """
        Trigger n8n workflow when analysis is complete
        """
        try:
            webhook_url = self.webhook_urls.get('analysis_complete')
            if not webhook_url:
                self.logger.warning("No webhook URL configured for analysis_complete")
                return False
            
            standardized_results = self.prepare_pipeline_results(pipeline_results)
            
            response = requests.post(
                webhook_url,
                json=standardized_results,
                headers={'Content-Type': 'application/json'},
                timeout=30
            )
            
            if response.status_code == 200:
                self.logger.info(f"Successfully triggered n8n workflow for analysis {standardized_results['metadata']['analysis_id']}")
                return True
            else:
                self.logger.error(f"Failed to trigger n8n workflow: HTTP {response.status_code}")
                return False
                
        except Exception as e:
            self.logger.error(f"Error triggering n8n workflow: {str(e)}")
            return False
    
    def trigger_high_priority_alert(self, pipeline_results: Dict) -> bool:
        """
        Trigger high-priority alert workflow for significant findings
        """
        try:
            high_priority_count = pipeline_results.get('high_priority_count', 0)
            threshold = self.notification_config.get('high_priority_threshold', 5)
            
            if high_priority_count >= threshold:
                webhook_url = self.webhook_urls.get('high_priority_alert')
                if webhook_url:
                    alert_payload = {
                        'analysis_id': pipeline_results.get('analysis_id'),
                        'high_priority_findings': high_priority_count,
                        'user_email': pipeline_results.get('user_email'),
                        'timestamp': datetime.now().isoformat(),
                        'summary': f"Found {high_priority_count} high-priority resistance mutations",
                        'files': {
                            'html_report': pipeline_results.get('html_report_path'),
                            'csv_summary': pipeline_results.get('csv_summary_path')
                        }
                    }
                    
                    response = requests.post(webhook_url, json=alert_payload, timeout=30)
                    return response.status_code == 200
            
            return True  # No alert needed
            
        except Exception as e:
            self.logger.error(f"Error triggering high-priority alert: {str(e)}")
            return False
    
    def trigger_quality_check_workflow(self, pipeline_results: Dict) -> bool:
        """
        Trigger quality check workflow for analysis validation
        """
        try:
            webhook_url = self.webhook_urls.get('quality_check')
            if not webhook_url:
                return True  # Quality check is optional
            
            quality_payload = {
                'analysis_id': pipeline_results.get('analysis_id'),
                'quality_score': pipeline_results.get('quality_score', 0.0),
                'genome_completeness': pipeline_results.get('genome_completeness', {}),
                'alignment_quality': pipeline_results.get('alignment_quality', {}),
                'timestamp': datetime.now().isoformat()
            }
            
            response = requests.post(webhook_url, json=quality_payload, timeout=30)
            return response.status_code == 200
            
        except Exception as e:
            self.logger.error(f"Error triggering quality check workflow: {str(e)}")
            return False
    
    def process_pipeline_completion(self, pipeline_results: Dict) -> Dict:
        """
        Main method to process pipeline completion and trigger all relevant workflows
        """
        results = {
            'analysis_complete': False,
            'high_priority_alert': False,
            'quality_check': False,
            'errors': []
        }
        
        try:
            # Trigger main analysis complete workflow
            results['analysis_complete'] = self.trigger_analysis_complete_workflow(pipeline_results)
            
            # Trigger high-priority alert if needed
            results['high_priority_alert'] = self.trigger_high_priority_alert(pipeline_results)
            
            # Trigger quality check workflow
            results['quality_check'] = self.trigger_quality_check_workflow(pipeline_results)
            
            self.logger.info(f"n8n workflow processing complete: {results}")
            
        except Exception as e:
            error_msg = f"Error in n8n workflow processing: {str(e)}"
            self.logger.error(error_msg)
            results['errors'].append(error_msg)
        
        return results

def integrate_with_pipeline(pipeline_results: Dict, config_file: str = "config/n8n_config.yaml") -> Dict:
    """
    Convenience function to integrate n8n workflows with pipeline results
    """
    n8n_manager = N8nIntegrationManager(config_file)
    return n8n_manager.process_pipeline_completion(pipeline_results)

# Example usage in pipeline
if __name__ == "__main__":
    # Example pipeline results
    example_results = {
        'genome_count': 10,
        'resistance_gene_count': 25,
        'high_priority_count': 7,
        'quality_score': 0.95,
        'user_email': 'researcher@university.edu',
        'html_report_path': 'reports/analysis_report.html',
        'csv_summary_path': 'reports/summary.csv',
        'processing_time': 45
    }
    
    # Trigger n8n workflows
    results = integrate_with_pipeline(example_results)
    print(f"n8n Integration Results: {results}")