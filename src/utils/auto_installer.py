"""
Automatic Tool Installation and Setup for GenomeAMRAnalyzer
Handles RGI, CARD database, and other dependencies automatically.
"""

import os
import sys
import subprocess
import shutil
import urllib.request
import json
import zipfile
import tarfile
from pathlib import Path
from typing import Optional, Dict, Any
import logging

logger = logging.getLogger(__name__)

class AutoInstaller:
    """Automatically install and configure required tools."""
    
    def __init__(self, install_dir: Optional[Path] = None):
        self.install_dir = install_dir or Path.home() / ".genomeamr"
        self.install_dir.mkdir(parents=True, exist_ok=True)
        self.tools_dir = self.install_dir / "tools"
        self.data_dir = self.install_dir / "data"
        self.tools_dir.mkdir(parents=True, exist_ok=True)
        self.data_dir.mkdir(parents=True, exist_ok=True)
        
        # Status file to track installations
        self.status_file = self.install_dir / "installation_status.json"
        
    def get_installation_status(self) -> Dict[str, Any]:
        """Get current installation status."""
        if self.status_file.exists():
            with self.status_file.open() as f:
                return json.load(f)
        return {}
    
    def update_installation_status(self, tool: str, status: Dict[str, Any]) -> None:
        """Update installation status for a tool."""
        current = self.get_installation_status()
        current[tool] = status
        with self.status_file.open("w") as f:
            json.dump(current, f, indent=2)
    
    def check_conda_available(self) -> bool:
        """Check if conda is available."""
        return shutil.which("conda") is not None
    
    def install_rgi_conda(self) -> bool:
        """Install RGI using conda if available."""
        if not self.check_conda_available():
            return False
        
        try:
            logger.info("Installing RGI via conda...")
            subprocess.run([
                "conda", "install", "-c", "bioconda", "-c", "conda-forge", 
                "rgi", "-y"
            ], check=True, capture_output=True)
            
            # Verify installation
            result = subprocess.run(["rgi", "--version"], capture_output=True, text=True)
            if result.returncode == 0:
                logger.info(f"RGI installed successfully: {result.stdout.strip()}")
                return True
        except subprocess.CalledProcessError as e:
            logger.warning(f"Conda RGI installation failed: {e}")
        
        return False
    
    def install_rgi_pip(self) -> bool:
        """Install RGI using pip as fallback."""
        try:
            logger.info("Installing RGI via pip...")
            subprocess.run([
                sys.executable, "-m", "pip", "install", "rgi"
            ], check=True, capture_output=True)
            
            # Verify installation
            result = subprocess.run(["rgi", "--version"], capture_output=True, text=True)
            if result.returncode == 0:
                logger.info(f"RGI installed successfully: {result.stdout.strip()}")
                return True
        except subprocess.CalledProcessError as e:
            logger.warning(f"Pip RGI installation failed: {e}")
        
        return False
    
    def download_card_database(self) -> bool:
        """Download and setup CARD database."""
        try:
            card_dir = self.data_dir / "card"
            card_dir.mkdir(parents=True, exist_ok=True)
            
            # Download CARD database
            logger.info("Downloading CARD database...")
            card_url = "https://card.mcmaster.ca/latest/data"
            card_file = card_dir / "card.tar.bz2"
            
            urllib.request.urlretrieve(card_url, card_file)
            
            # Extract database
            logger.info("Extracting CARD database...")
            with tarfile.open(card_file, "r:bz2") as tar:
                tar.extractall(card_dir)
            
            # Find the JSON file
            json_files = list(card_dir.glob("**/*.json"))
            if json_files:
                card_json = json_files[0]
                # Create symlink for easy access
                easy_path = card_dir / "card.json"
                if easy_path.exists():
                    easy_path.unlink()
                easy_path.symlink_to(card_json)
                
                logger.info(f"CARD database ready: {easy_path}")
                return True
            else:
                logger.error("No JSON file found in CARD database")
                return False
                
        except Exception as e:
            logger.error(f"CARD database download failed: {e}")
            return False
    
    def setup_rgi_database(self) -> bool:
        """Setup RGI with CARD database."""
        try:
            card_json = self.data_dir / "card" / "card.json"
            if not card_json.exists():
                logger.error("CARD database not found")
                return False
            
            logger.info("Setting up RGI database...")
            subprocess.run([
                "rgi", "card_annotation", "--input", str(card_json)
            ], check=True, capture_output=True)
            
            logger.info("RGI database setup complete")
            return True
            
        except subprocess.CalledProcessError as e:
            logger.error(f"RGI database setup failed: {e}")
            return False
    
    def install_rgi_complete(self) -> bool:
        """Complete RGI installation including database."""
        status = self.get_installation_status()
        
        # Check if already installed
        if status.get("rgi", {}).get("installed", False):
            if shutil.which("rgi"):
                logger.info("RGI already installed and available")
                return True
        
        # Try conda first, then pip
        rgi_installed = self.install_rgi_conda() or self.install_rgi_pip()
        
        if not rgi_installed:
            logger.error("Failed to install RGI")
            return False
        
        # Download and setup CARD database
        if not self.download_card_database():
            logger.error("Failed to download CARD database")
            return False
        
        if not self.setup_rgi_database():
            logger.error("Failed to setup RGI database")
            return False
        
        # Update status
        self.update_installation_status("rgi", {
            "installed": True,
            "version": subprocess.run(["rgi", "--version"], capture_output=True, text=True).stdout.strip(),
            "database_path": str(self.data_dir / "card" / "card.json")
        })
        
        logger.info("RGI installation complete!")
        return True
    
    def ensure_tools_installed(self) -> Dict[str, bool]:
        """Ensure all required tools are installed."""
        results = {}
        
        # Install RGI if not available
        if not shutil.which("rgi"):
            logger.info("RGI not found, installing automatically...")
            results["rgi"] = self.install_rgi_complete()
        else:
            logger.info("RGI already available")
            results["rgi"] = True
        
        # Note: We'll handle EMBOSS internally, so no installation needed
        results["emboss"] = True  # Internal implementation
        
        return results
    
    def get_tool_paths(self) -> Dict[str, str]:
        """Get paths to installed tools."""
        return {
            "rgi": shutil.which("rgi") or "rgi",
            "card_database": str(self.data_dir / "card" / "card.json"),
            "install_dir": str(self.install_dir)
        }


def auto_setup() -> bool:
    """Perform automatic setup of all tools."""
    installer = AutoInstaller()
    
    logger.info("Starting automatic tool installation...")
    results = installer.ensure_tools_installed()
    
    if all(results.values()):
        logger.info("All tools installed successfully!")
        logger.info(f"Installation directory: {installer.install_dir}")
        return True
    else:
        logger.error("Some tools failed to install:")
        for tool, success in results.items():
            status = "✓" if success else "✗"
            logger.error(f"  {status} {tool}")
        return False


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    
    if len(sys.argv) > 1 and sys.argv[1] == "--check":
        # Just check status
        installer = AutoInstaller()
        status = installer.get_installation_status()
        paths = installer.get_tool_paths()
        
        print("Installation Status:")
        print(json.dumps(status, indent=2))
        print("\nTool Paths:")
        print(json.dumps(paths, indent=2))
    else:
        # Perform installation
        success = auto_setup()
        sys.exit(0 if success else 1)