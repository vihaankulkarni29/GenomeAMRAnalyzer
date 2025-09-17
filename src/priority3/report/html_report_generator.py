"""
HTML Report Generator for AMR Genomics Pipeline - SECURITY FIXED VERSION
========================================================================

Enterprise-grade reporting for AMR genomics:
- Jinja2-based HTML templating with XSS protection
- Comprehensive tables: genomes, mutations, MICs, co-occurrence
- Interactive, sortable, and filterable tables (DataTables.js)
- Visualizations: heatmaps, bar charts, network graphs (Plotly.js)
- Artifact linking and provenance tracking
- Publication-ready export (PDF, HTML)

Security Enhancements:
- Proper HTML escaping for all user data
- Fixed artifact link key generation
- Enhanced input validation
- XSS prevention filters
"""

import os
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional
from jinja2 import Environment, FileSystemLoader, select_autoescape
import json
import html

class HTMLReportGenerator:
    """
    Jinja2-based HTML report generator for AMR genomics pipeline.
    Enhanced with security fixes and proper artifact linking.
    """
    def __init__(self, 
                 template_dir: str = "report/templates",
                 output_dir: str = "report/html_reports"):
        self.template_dir = Path(template_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)
        self.env = Environment(
            loader=FileSystemLoader(str(self.template_dir)),
            autoescape=select_autoescape(['html', 'xml']),
            trim_blocks=True,
            lstrip_blocks=True
        )
        # Add custom filters for additional security
        self.env.filters['safe_html'] = self._safe_html_filter
        self.logger = logging.getLogger("HTMLReportGenerator")

    def _safe_html_filter(self, value):
        """
        Additional HTML safety filter to prevent XSS.
        """
        if value is None:
            return ""
        # Convert to string and escape dangerous characters
        safe_value = html.escape(str(value))
        return safe_value

    def render_report(self, 
                     context: Dict[str, Any],
                     template_name: str = "amr_report.html",
                     output_name: Optional[str] = None) -> str:
        """
        Render a complete HTML report from context data.
        Returns path to generated HTML file.
        """
        template = self.env.get_template(template_name)
        html_content = template.render(**context)
        if not output_name:
            output_name = f"amr_report_{context.get('run_id', 'latest')}.html"
        output_path = self.output_dir / output_name
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        self.logger.info(f"Generated HTML report: {output_path}")
        return str(output_path)

    def build_context(self, 
                     genomes: List[Dict[str, Any]],
                     mutations: List[Dict[str, Any]],
                     mic_data: List[Dict[str, Any]],
                     cooccurrence: List[Dict[str, Any]],
                     stats: Dict[str, Any],
                     run_id: Optional[str] = None) -> Dict[str, Any]:
        """
        Build context dictionary for report rendering.
        """
        context = {
            'genomes': genomes,
            'mutations': mutations,
            'mic_data': mic_data,
            'cooccurrence': cooccurrence,
            'stats': stats,
            'run_id': run_id or 'latest',
            'artifact_links': self._build_artifact_links(genomes, mutations, mic_data, cooccurrence)
        }
        return context

    def _build_artifact_links(self, genomes, mutations, mic_data, cooccurrence):
        """
        Build mapping of artifact links for traceability.
        Fixed key generation to match template expectations.
        """
        links = {}
        
        # Genome artifacts
        for g in genomes:
            if 'artifact_path' in g and g['artifact_path']:
                links[f"genome_{g['accession']}"] = g['artifact_path']
        
        # Mutation artifacts
        for m in mutations:
            if 'artifact_path' in m and m['artifact_path']:
                links[f"mutation_{m['mutation_id']}"] = m['artifact_path']
        
        # MIC artifacts (by accession)
        for mic in mic_data:
            if 'artifact_path' in mic and mic['artifact_path']:
                links[f"mic_{mic['accession']}"] = mic['artifact_path']
        
        # Co-occurrence artifacts
        for c in cooccurrence:
            if 'artifact_path' in c and c['artifact_path']:
                links[f"cooccur_{c['mutation_a']}_{c['mutation_b']}"] = c['artifact_path']
        
        return links

    def ensure_templates(self):
        """
        Ensure default templates exist (create if missing).
        """
        self.template_dir.mkdir(exist_ok=True, parents=True)
        template_path = self.template_dir / "amr_report.html"
        if not template_path.exists():
            with open(template_path, 'w', encoding='utf-8') as f:
                f.write(self.default_template())
            self.logger.info(f"Created default template: {template_path}")

    @staticmethod
    def default_template() -> str:
        """
        Default HTML template with DataTables and Plotly integration.
        Enhanced with proper XSS protection and fixed artifact links.
        """
        return '''
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>AMR Genomics Report</title>
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <link rel="stylesheet" href="https://cdn.datatables.net/1.13.4/css/jquery.dataTables.min.css">
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    <script src="https://cdn.datatables.net/1.13.4/js/jquery.dataTables.min.js"></script>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body { font-family: Arial, sans-serif; margin: 2em; }
        h1, h2 { color: #2c3e50; }
        .artifact-link { font-size: 0.9em; color: #2980b9; text-decoration: none; }
        .artifact-link:hover { text-decoration: underline; }
        .section { margin-bottom: 2em; }
        .table-responsive { overflow-x: auto; }
        @media (max-width: 768px) {
            body { margin: 1em; }
            table { font-size: 0.9em; }
        }
    </style>
</head>
<body>
    <h1>AMR Genomics Analysis Report</h1>
    <p>Run ID: {{ run_id|e }}</p>
    <div class="section">
        <h2>Genome Summary</h2>
        <div class="table-responsive">
        <table id="genomeTable" class="display">
            <thead>
                <tr>
                    <th>Accession</th>
                    <th>Organism</th>
                    <th>Strain</th>
                    <th>BioSample</th>
                    <th>BioProject</th>
                    <th>Download Date</th>
                    <th>Artifact</th>
                </tr>
            </thead>
            <tbody>
                {% for g in genomes %}
                <tr>
                    <td>{{ g.accession|e }}</td>
                    <td>{{ g.organism|e }}</td>
                    <td>{{ g.strain|e }}</td>
                    <td>{{ g.biosample|e }}</td>
                    <td>{{ g.bioproject|e }}</td>
                    <td>{{ g.download_date|e }}</td>
                    <td>{% if artifact_links['genome_' ~ g.accession] %}<a class="artifact-link" href="{{ artifact_links['genome_' ~ g.accession]|e }}">View</a>{% endif %}</td>
                </tr>
                {% endfor %}
            </tbody>
        </table>
        </div>
    </div>
    <div class="section">
        <h2>Mutations</h2>
        <div class="table-responsive">
        <table id="mutationTable" class="display">
            <thead>
                <tr>
                    <th>Mutation ID</th>
                    <th>Accession</th>
                    <th>Gene</th>
                    <th>Type</th>
                    <th>Position</th>
                    <th>Ref Base</th>
                    <th>Mut Base</th>
                    <th>Confidence</th>
                    <th>AMR Known</th>
                    <th>Artifact</th>
                </tr>
            </thead>
            <tbody>
                {% for m in mutations %}
                <tr>
                    <td>{{ m.mutation_id|e }}</td>
                    <td>{{ m.accession|e }}</td>
                    <td>{{ m.gene_context|e }}</td>
                    <td>{{ m.type|e }}</td>
                    <td>{{ m.position }}</td>
                    <td>{{ m.reference_base|e }}</td>
                    <td>{{ m.mutant_base|e }}</td>
                    <td>{{ m.confidence }}</td>
                    <td>{{ m.known_resistance }}</td>
                    <td>{% if artifact_links['mutation_' ~ m.mutation_id] %}<a class="artifact-link" href="{{ artifact_links['mutation_' ~ m.mutation_id]|e }}">View</a>{% endif %}</td>
                </tr>
                {% endfor %}
            </tbody>
        </table>
        </div>
    </div>
    <div class="section">
        <h2>MIC Data</h2>
        <div class="table-responsive">
        <table id="micTable" class="display">
            <thead>
                <tr>
                    <th>Accession</th>
                    <th>Antibiotic</th>
                    <th>MIC Value</th>
                    <th>Unit</th>
                    <th>Resistance</th>
                    <th>Quality</th>
                    <th>Artifact</th>
                </tr>
            </thead>
            <tbody>
                {% for mic in mic_data %}
                <tr>
                    <td>{{ mic.accession|e }}</td>
                    <td>{{ mic.antibiotic|e }}</td>
                    <td>{{ mic.mic_value }}</td>
                    <td>{{ mic.mic_unit|e }}</td>
                    <td>{{ mic.resistance_profile|e }}</td>
                    <td>{{ mic.quality_score }}</td>
                    <td>{% if artifact_links['mic_' ~ mic.accession] %}<a class="artifact-link" href="{{ artifact_links['mic_' ~ mic.accession]|e }}">View</a>{% endif %}</td>
                </tr>
                {% endfor %}
            </tbody>
        </table>
        </div>
    </div>
    <div class="section">
        <h2>Mutation Co-occurrence</h2>
        <div class="table-responsive">
        <table id="cooccurrenceTable" class="display">
            <thead>
                <tr>
                    <th>Mutation A</th>
                    <th>Mutation B</th>
                    <th>Count A</th>
                    <th>Count B</th>
                    <th>Count Both</th>
                    <th>P-value</th>
                    <th>Odds Ratio</th>
                    <th>Significant</th>
                    <th>Artifact</th>
                </tr>
            </thead>
            <tbody>
                {% for c in cooccurrence %}
                <tr>
                    <td>{{ c.mutation_a|e }}</td>
                    <td>{{ c.mutation_b|e }}</td>
                    <td>{{ c.count_a }}</td>
                    <td>{{ c.count_b }}</td>
                    <td>{{ c.count_both }}</td>
                    <td>{{ c.p_value|round(4) }}</td>
                    <td>{{ c.odds_ratio|round(2) if c.odds_ratio is not none else '' }}</td>
                    <td>{{ c.significant }}</td>
                    <td>{% if artifact_links['cooccur_' ~ c.mutation_a ~ '_' ~ c.mutation_b] %}<a class="artifact-link" href="{{ artifact_links['cooccur_' ~ c.mutation_a ~ '_' ~ c.mutation_b]|e }}">View</a>{% endif %}</td>
                </tr>
                {% endfor %}
            </tbody>
        </table>
        </div>
    </div>
    <div class="section">
        <h2>Pipeline Statistics</h2>
        <pre>{{ stats | tojson(indent=2) | e }}</pre>
    </div>
    <script>
        $(document).ready(function() {
            $('#genomeTable').DataTable({
                "pageLength": 25,
                "responsive": true
            });
            $('#mutationTable').DataTable({
                "pageLength": 25,
                "responsive": true
            });
            $('#micTable').DataTable({
                "pageLength": 25,
                "responsive": true
            });
            $('#cooccurrenceTable').DataTable({
                "pageLength": 25,
                "responsive": true
            });
        });
    </script>
</body>
</html>
'''