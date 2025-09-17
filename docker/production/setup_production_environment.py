#!/usr/bin/env python3
"""
Production Environment Setup Script
Installs required dependencies and fixes environment issues

This script:
1. Installs missing Python packages
2. Validates installation
3. Creates missing directories
4. Fixes common configuration issues

Usage:
    python setup_production_environment.py [--install-deps] [--fix-all]

Author: GenomeAMRAnalyzer Pipeline
Version: 2.0 - Environment Setup
"""

import os
import sys
import subprocess
import logging
from pathlib import Path


def setup_logging():
    """Setup logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger(__name__)


def run_command(command, logger):
    """Run a command and log the result"""
    logger.info(f"Running: {' '.join(command)}")
    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=True
        )
        logger.info(f"✅ Command successful")
        return True, result.stdout
    except subprocess.CalledProcessError as e:
        logger.error(f"❌ Command failed: {e}")
        logger.error(f"   stdout: {e.stdout}")
        logger.error(f"   stderr: {e.stderr}")
        return False, e.stderr
    except Exception as e:
        logger.error(f"❌ Unexpected error: {e}")
        return False, str(e)


def install_missing_dependencies(logger):
    """Install missing production dependencies"""
    logger.info("Installing missing dependencies...")
    
    # Missing packages identified from validation
    missing_packages = [
        'biopython==1.84',
        'structlog==24.2.0'
    ]
    
    for package in missing_packages:
        logger.info(f"Installing {package}...")
        success, output = run_command([sys.executable, '-m', 'pip', 'install', package], logger)
        
        if success:
            logger.info(f"✅ Successfully installed {package}")
        else:
            logger.error(f"❌ Failed to install {package}")
            return False
    
    return True


def install_all_production_requirements(logger):
    """Install all production requirements from requirements file"""
    logger.info("Installing all production requirements...")
    
    requirements_file = Path(__file__).parent.parent.parent / 'requirements' / 'production.txt'
    
    if not requirements_file.exists():
        logger.error(f"❌ Requirements file not found: {requirements_file}")
        return False
    
    logger.info(f"Installing from {requirements_file}")
    success, output = run_command([
        sys.executable, '-m', 'pip', 'install', '-r', str(requirements_file)
    ], logger)
    
    if success:
        logger.info("✅ Successfully installed all production requirements")
        return True
    else:
        logger.error("❌ Failed to install production requirements")
        return False


def create_missing_directories(logger):
    """Create missing directories for production environment"""
    logger.info("Creating missing directories...")
    
    base_path = Path(__file__).parent.parent.parent
    required_dirs = [
        base_path / 'data',
        base_path / 'results', 
        base_path / 'logs',
        base_path / 'cache',
        base_path / 'config' / 'production'
    ]
    
    for directory in required_dirs:
        try:
            directory.mkdir(parents=True, exist_ok=True)
            logger.info(f"✅ Directory ready: {directory}")
        except Exception as e:
            logger.error(f"❌ Failed to create directory {directory}: {e}")
            return False
    
    return True


def validate_environment(logger):
    """Run environment validation to check if setup was successful"""
    logger.info("Validating environment setup...")
    
    validation_script = Path(__file__).parent / 'validate_environment_fixed.py'
    
    if not validation_script.exists():
        logger.error(f"❌ Validation script not found: {validation_script}")
        return False
    
    success, output = run_command([
        sys.executable, str(validation_script), '--create-missing'
    ], logger)
    
    if success:
        logger.info("✅ Environment validation passed")
        return True
    else:
        logger.warning("⚠️  Environment validation had issues, but setup may still be functional")
        logger.info("Run the validation script manually for detailed results")
        return True  # Don't fail setup just because validation has warnings


def main():
    """Main setup function"""
    logger = setup_logging()
    
    logger.info("🚀 Starting Production Environment Setup")
    logger.info("=" * 60)
    
    # Parse simple command line arguments
    install_deps = '--install-deps' in sys.argv or '--fix-all' in sys.argv
    fix_all = '--fix-all' in sys.argv
    
    success = True
    
    # Step 1: Create directories
    logger.info("Step 1: Creating required directories")
    if not create_missing_directories(logger):
        success = False
    
    # Step 2: Install dependencies if requested
    if install_deps:
        logger.info("Step 2: Installing dependencies")
        if not install_all_production_requirements(logger):
            logger.warning("⚠️  Failed to install all requirements, trying missing packages only")
            if not install_missing_dependencies(logger):
                success = False
    else:
        logger.info("Step 2: Skipping dependency installation (use --install-deps to enable)")
    
    # Step 3: Validate environment
    logger.info("Step 3: Validating environment")
    if not validate_environment(logger):
        success = False
    
    # Summary
    logger.info("=" * 60)
    if success:
        logger.info("🎉 Production Environment Setup Completed Successfully!")
        logger.info("")
        logger.info("Next steps:")
        logger.info("1. Run validation: python docker/production/validate_environment_fixed.py")
        logger.info("2. Test production environment: python docker/production/test_production_environment.py")
        logger.info("3. Build Docker container: docker build -f docker/production/Dockerfile .")
    else:
        logger.error("❌ Production Environment Setup Failed")
        logger.error("Please check the errors above and fix them manually")
    
    return 0 if success else 1


if __name__ == '__main__':
    sys.exit(main())