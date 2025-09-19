#!/bin/bash

# GenomeAMRAnalyzer v1.0.0 Release Script
# Automated Docker Hub publishing and Git tagging for production release

set -e  # Exit on any error

# Configuration
DOCKER_IMAGE="vihaankulkarni29/genomeamranalyzer"
VERSION="1.0.0"
LATEST_TAG="latest"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to check prerequisites
check_prerequisites() {
    log_info "Checking prerequisites..."
    
    # Check if Docker is installed and running
    if ! command -v docker &> /dev/null; then
        log_error "Docker is not installed or not in PATH"
        exit 1
    fi
    
    if ! docker info &> /dev/null; then
        log_error "Docker daemon is not running"
        exit 1
    fi
    
    # Check if Git is installed
    if ! command -v git &> /dev/null; then
        log_error "Git is not installed or not in PATH"
        exit 1
    fi
    
    # Check if we're in a Git repository
    if ! git rev-parse --git-dir &> /dev/null; then
        log_error "Not in a Git repository"
        exit 1
    fi
    
    # Check if Docker Hub credentials are available
    log_info "Checking Docker Hub authentication..."
    if ! docker info | grep -q "Username:"; then
        log_warning "Not logged into Docker Hub. Please run 'docker login' first"
        read -p "Do you want to login now? (y/n): " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            docker login
        else
            log_error "Docker Hub login required for publishing"
            exit 1
        fi
    fi
    
    log_success "Prerequisites check passed"
}

# Function to validate current state
validate_state() {
    log_info "Validating current repository state..."
    
    # Check if working directory is clean
    if ! git diff-index --quiet HEAD --; then
        log_error "Working directory is not clean. Please commit or stash changes"
        git status --short
        exit 1
    fi
    
    # Check if we're on main branch
    current_branch=$(git branch --show-current)
    if [ "$current_branch" != "main" ]; then
        log_warning "Not on main branch (currently on: $current_branch)"
        read -p "Do you want to continue anyway? (y/n): " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            log_error "Aborted by user"
            exit 1
        fi
    fi
    
    # Check if tag already exists
    if git tag -l | grep -q "^v${VERSION}$"; then
        log_error "Tag v${VERSION} already exists"
        exit 1
    fi
    
    log_success "Repository state validation passed"
}

# Function to build Docker image
build_docker_image() {
    log_info "Building Docker image..."
    
    # Build image with version tag
    docker build -t "${DOCKER_IMAGE}:${VERSION}" -t "${DOCKER_IMAGE}:${LATEST_TAG}" .
    
    if [ $? -eq 0 ]; then
        log_success "Docker image built successfully"
    else
        log_error "Docker image build failed"
        exit 1
    fi
}

# Function to test Docker image
test_docker_image() {
    log_info "Testing Docker image..."
    
    # Basic functionality test
    docker run --rm "${DOCKER_IMAGE}:${VERSION}" python -c "
import sys
sys.path.append('/app/src')
try:
    from simplified_card_integrator import SimplifiedCARDIntegrator
    from fasta_aa_extractor_integration import FastaAAExtractorIntegration
    print('✓ Core modules import successfully')
except ImportError as e:
    print(f'✗ Import error: {e}')
    sys.exit(1)
"
    
    if [ $? -eq 0 ]; then
        log_success "Docker image test passed"
    else
        log_error "Docker image test failed"
        exit 1
    fi
}

# Function to publish to Docker Hub
publish_docker_image() {
    log_info "Publishing Docker image to Docker Hub..."
    
    # Push version tag
    log_info "Pushing version tag: ${VERSION}"
    docker push "${DOCKER_IMAGE}:${VERSION}"
    
    # Push latest tag
    log_info "Pushing latest tag"
    docker push "${DOCKER_IMAGE}:${LATEST_TAG}"
    
    if [ $? -eq 0 ]; then
        log_success "Docker image published successfully"
    else
        log_error "Docker image publishing failed"
        exit 1
    fi
}

# Function to create Git tag
create_git_tag() {
    log_info "Creating Git tag v${VERSION}..."
    
    # Create annotated tag
    git tag -a "v${VERSION}" -m "Release v${VERSION} - Production ready GenomeAMRAnalyzer

Features:
- Comprehensive AMR gene analysis pipeline
- RGI to Abricate migration for improved performance
- Docker containerization for easy deployment
- Full integration test suite
- Production-ready v1.0.0 release

Docker Hub: ${DOCKER_IMAGE}:${VERSION}"
    
    # Push tag to remote
    git push origin "v${VERSION}"
    
    if [ $? -eq 0 ]; then
        log_success "Git tag created and pushed successfully"
    else
        log_error "Git tag creation/push failed"
        exit 1
    fi
}

# Function to generate release summary
generate_release_summary() {
    log_info "Generating release summary..."
    
    echo ""
    echo "=================================="
    echo "   GenomeAMRAnalyzer v${VERSION} RELEASE SUMMARY"
    echo "=================================="
    echo ""
    echo "🚀 Release Status: SUCCESS"
    echo ""
    echo "📦 Docker Image:"
    echo "   - Repository: ${DOCKER_IMAGE}"
    echo "   - Version Tag: ${VERSION}"
    echo "   - Latest Tag: ${LATEST_TAG}"
    echo "   - Pull Command: docker pull ${DOCKER_IMAGE}:${VERSION}"
    echo ""
    echo "🏷️  Git Tag:"
    echo "   - Tag: v${VERSION}"
    echo "   - Branch: $(git branch --show-current)"
    echo "   - Commit: $(git rev-parse --short HEAD)"
    echo ""
    echo "🔗 Links:"
    echo "   - Docker Hub: https://hub.docker.com/r/${DOCKER_IMAGE/\//-}"
    echo "   - GitHub: https://github.com/vihaankulkarni29/GenomeAMRAnalyzer"
    echo ""
    echo "📋 Next Steps:"
    echo "   1. Update Docker Hub repository description"
    echo "   2. Create GitHub release from tag v${VERSION}"
    echo "   3. Update project documentation"
    echo "   4. Announce release to users"
    echo ""
    echo "=================================="
}

# Main execution flow
main() {
    echo ""
    echo "🧬 GenomeAMRAnalyzer v${VERSION} Release Automation"
    echo "=================================================="
    echo ""
    
    # Confirmation prompt
    log_warning "This will:"
    echo "  • Build Docker image ${DOCKER_IMAGE}:${VERSION}"
    echo "  • Test the image functionality"
    echo "  • Publish to Docker Hub"
    echo "  • Create and push Git tag v${VERSION}"
    echo ""
    read -p "Do you want to proceed with the release? (y/n): " -n 1 -r
    echo ""
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        log_info "Release aborted by user"
        exit 0
    fi
    
    # Execute release steps
    check_prerequisites
    validate_state
    build_docker_image
    test_docker_image
    publish_docker_image
    create_git_tag
    generate_release_summary
    
    log_success "🎉 GenomeAMRAnalyzer v${VERSION} released successfully!"
}

# Run main function
main "$@"