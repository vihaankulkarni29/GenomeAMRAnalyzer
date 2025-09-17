"""
Export Manager
--------------
Handles standardized export of processed data in multiple formats (CSV, JSON, FASTA, Parquet).
Ensures robust, metadata-rich, and researcher-friendly outputs.
"""

import os
import pandas as pd
import json
from typing import List, Dict, Any

class ExportManagerError(Exception):
    pass

class ExportManager:
    def __init__(self, output_dir: str):
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

    def export_csv(self, data: List[Dict[str, Any]], filename: str):
        try:
            df = pd.DataFrame(data)
            df.to_csv(os.path.join(self.output_dir, filename), index=False)
        except Exception as e:
            raise ExportManagerError(f"CSV export failed: {e}")

    def export_json(self, data: List[Dict[str, Any]], filename: str):
        try:
            with open(os.path.join(self.output_dir, filename), 'w') as f:
                json.dump(data, f, indent=2)
        except Exception as e:
            raise ExportManagerError(f"JSON export failed: {e}")

    def export_fasta(self, records: List[Dict[str, str]], filename: str):
        try:
            with open(os.path.join(self.output_dir, filename), 'w') as f:
                for rec in records:
                    f.write(f">{rec['id']}\n{rec['seq']}\n")
        except Exception as e:
            raise ExportManagerError(f"FASTA export failed: {e}")

    def export_parquet(self, data: List[Dict[str, Any]], filename: str):
        try:
            df = pd.DataFrame(data)
            df.to_parquet(os.path.join(self.output_dir, filename))
        except Exception as e:
            raise ExportManagerError(f"Parquet export failed: {e}")
