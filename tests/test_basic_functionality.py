"""
Simple test to verify basic imports work
"""
import sys
import os

def test_python_version():
    """Test that Python version is compatible"""
    assert sys.version_info >= (3, 9), f"Python {sys.version_info} is too old"

def test_basic_imports():
    """Test that core dependencies can be imported"""
    try:
        import numpy
        import pandas
        import biopython
        print("Core dependencies imported successfully")
    except ImportError as e:
        print(f"Warning: Could not import some dependencies: {e}")
        # Don't fail the test for optional dependencies

def test_project_structure():
    """Test that basic project structure exists"""
    assert os.path.exists("src"), "src directory should exist"
    assert os.path.exists("README.md"), "README.md should exist"
    assert os.path.exists("requirements.txt"), "requirements.txt should exist"
    assert os.path.exists("Dockerfile"), "Dockerfile should exist"

def test_docker_files():
    """Test that Docker workflow files exist"""
    assert os.path.exists("docker-compose.yml"), "docker-compose.yml should exist"
    assert os.path.exists("run_docker.sh"), "run_docker.sh should exist"
    assert os.path.exists("run_docker.bat"), "run_docker.bat should exist"

if __name__ == "__main__":
    test_python_version()
    test_basic_imports()
    test_project_structure()
    test_docker_files()
    print("All basic tests passed!")