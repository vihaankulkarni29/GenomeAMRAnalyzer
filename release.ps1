# GenomeAMRAnalyzer v1.0.0 Release Script (PowerShell)
# Automated Docker Hub publishing and Git tagging for production release

param(
    [switch]$Force,
    [string]$Version = "1.0.0"
)

# Configuration
$DOCKER_IMAGE = "vihaankulkarni29/genomeamranalyzer"
$LATEST_TAG = "latest"

# Color functions for PowerShell
function Write-ColorOutput {
    param(
        [string]$Message,
        [string]$Color = "White"
    )
    Write-Host $Message -ForegroundColor $Color
}

function Write-Info { param([string]$Message) Write-ColorOutput "[INFO] $Message" -Color Cyan }
function Write-Success { param([string]$Message) Write-ColorOutput "[SUCCESS] $Message" -Color Green }
function Write-Warning { param([string]$Message) Write-ColorOutput "[WARNING] $Message" -Color Yellow }
function Write-Error { param([string]$Message) Write-ColorOutput "[ERROR] $Message" -Color Red }

# Function to check prerequisites
function Test-Prerequisites {
    Write-Info "Checking prerequisites..."
    
    # Check if Docker is installed and running
    try {
        $null = Get-Command docker -ErrorAction Stop
        $null = docker info 2>$null
    }
    catch {
        Write-Error "Docker is not installed or not running"
        exit 1
    }
    
    # Check if Git is installed
    try {
        $null = Get-Command git -ErrorAction Stop
    }
    catch {
        Write-Error "Git is not installed or not in PATH"
        exit 1
    }
    
    # Check if we're in a Git repository
    try {
        $null = git rev-parse --git-dir 2>$null
    }
    catch {
        Write-Error "Not in a Git repository"
        exit 1
    }
    
    # Check Docker Hub authentication
    Write-Info "Checking Docker Hub authentication..."
    $dockerInfo = docker info 2>$null | Out-String
    if (-not $dockerInfo.Contains("Username:")) {
        Write-Warning "Not logged into Docker Hub. Please run 'docker login' first"
        $response = Read-Host "Do you want to login now? (y/n)"
        if ($response -eq "y" -or $response -eq "Y") {
            docker login
            if ($LASTEXITCODE -ne 0) {
                Write-Error "Docker login failed"
                exit 1
            }
        }
        else {
            Write-Error "Docker Hub login required for publishing"
            exit 1
        }
    }
    
    Write-Success "Prerequisites check passed"
}

# Function to validate current state
function Test-RepositoryState {
    Write-Info "Validating current repository state..."
    
    # Check if working directory is clean
    $gitStatus = git status --porcelain
    if ($gitStatus) {
        Write-Error "Working directory is not clean. Please commit or stash changes"
        git status --short
        exit 1
    }
    
    # Check if we're on main branch
    $currentBranch = git branch --show-current
    if ($currentBranch -ne "main" -and -not $Force) {
        Write-Warning "Not on main branch (currently on: $currentBranch)"
        $response = Read-Host "Do you want to continue anyway? (y/n)"
        if ($response -ne "y" -and $response -ne "Y") {
            Write-Error "Aborted by user"
            exit 1
        }
    }
    
    # Check if tag already exists
    $existingTag = git tag -l | Where-Object { $_ -eq "v$Version" }
    if ($existingTag -and -not $Force) {
        Write-Error "Tag v$Version already exists"
        exit 1
    }
    
    Write-Success "Repository state validation passed"
}

# Function to build Docker image
function Build-DockerImage {
    Write-Info "Building Docker image..."
    
    # Build image with version tag
    docker build -t "${DOCKER_IMAGE}:${Version}" -t "${DOCKER_IMAGE}:${LATEST_TAG}" .
    
    if ($LASTEXITCODE -eq 0) {
        Write-Success "Docker image built successfully"
    }
    else {
        Write-Error "Docker image build failed"
        exit 1
    }
}

# Function to test Docker image
function Test-DockerImage {
    Write-Info "Testing Docker image..."
    
    # Basic functionality test
    $testScript = @"
import sys
sys.path.append('/app/src')
try:
    from simplified_card_integrator import SimplifiedCARDIntegrator
    from fasta_aa_extractor_integration import FastaAAExtractorIntegration
    print('✓ Core modules import successfully')
except ImportError as e:
    print(f'✗ Import error: {e}')
    sys.exit(1)
"@
    
    docker run --rm "${DOCKER_IMAGE}:${Version}" python -c $testScript
    
    if ($LASTEXITCODE -eq 0) {
        Write-Success "Docker image test passed"
    }
    else {
        Write-Error "Docker image test failed"
        exit 1
    }
}

# Function to publish to Docker Hub
function Publish-DockerImage {
    Write-Info "Publishing Docker image to Docker Hub..."
    
    # Push version tag
    Write-Info "Pushing version tag: $Version"
    docker push "${DOCKER_IMAGE}:${Version}"
    
    if ($LASTEXITCODE -ne 0) {
        Write-Error "Failed to push version tag"
        exit 1
    }
    
    # Push latest tag
    Write-Info "Pushing latest tag"
    docker push "${DOCKER_IMAGE}:${LATEST_TAG}"
    
    if ($LASTEXITCODE -eq 0) {
        Write-Success "Docker image published successfully"
    }
    else {
        Write-Error "Docker image publishing failed"
        exit 1
    }
}

# Function to create Git tag
function New-GitTag {
    Write-Info "Creating Git tag v$Version..."
    
    # Create annotated tag
    $tagMessage = @"
Release v$Version - Production ready GenomeAMRAnalyzer

Features:
- Comprehensive AMR gene analysis pipeline
- RGI to Abricate migration for improved performance
- Docker containerization for easy deployment
- Full integration test suite
- Production-ready v1.0.0 release

Docker Hub: ${DOCKER_IMAGE}:${Version}
"@
    
    git tag -a "v$Version" -m $tagMessage
    
    if ($LASTEXITCODE -ne 0) {
        Write-Error "Git tag creation failed"
        exit 1
    }
    
    # Push tag to remote
    git push origin "v$Version"
    
    if ($LASTEXITCODE -eq 0) {
        Write-Success "Git tag created and pushed successfully"
    }
    else {
        Write-Error "Git tag push failed"
        exit 1
    }
}

# Function to generate release summary
function Show-ReleaseSummary {
    Write-Info "Generating release summary..."
    
    $currentBranch = git branch --show-current
    $commitHash = git rev-parse --short HEAD
    $dockerHubUrl = "https://hub.docker.com/r/$($DOCKER_IMAGE.Replace('/', '-'))"
    
    Write-Host ""
    Write-Host "==================================" -ForegroundColor Magenta
    Write-Host "   GenomeAMRAnalyzer v$Version RELEASE SUMMARY" -ForegroundColor Magenta
    Write-Host "==================================" -ForegroundColor Magenta
    Write-Host ""
    Write-Host "🚀 Release Status: " -NoNewline
    Write-Host "SUCCESS" -ForegroundColor Green
    Write-Host ""
    Write-Host "📦 Docker Image:" -ForegroundColor Yellow
    Write-Host "   - Repository: $DOCKER_IMAGE"
    Write-Host "   - Version Tag: $Version"
    Write-Host "   - Latest Tag: $LATEST_TAG"
    Write-Host "   - Pull Command: " -NoNewline
    Write-Host "docker pull ${DOCKER_IMAGE}:${Version}" -ForegroundColor Cyan
    Write-Host ""
    Write-Host "🏷️  Git Tag:" -ForegroundColor Yellow
    Write-Host "   - Tag: v$Version"
    Write-Host "   - Branch: $currentBranch"
    Write-Host "   - Commit: $commitHash"
    Write-Host ""
    Write-Host "🔗 Links:" -ForegroundColor Yellow
    Write-Host "   - Docker Hub: $dockerHubUrl"
    Write-Host "   - GitHub: https://github.com/vihaankulkarni29/GenomeAMRAnalyzer"
    Write-Host ""
    Write-Host "📋 Next Steps:" -ForegroundColor Yellow
    Write-Host "   1. Update Docker Hub repository description"
    Write-Host "   2. Create GitHub release from tag v$Version"
    Write-Host "   3. Update project documentation"
    Write-Host "   4. Announce release to users"
    Write-Host ""
    Write-Host "==================================" -ForegroundColor Magenta
}

# Main execution flow
function Start-Release {
    Write-Host ""
    Write-Host "🧬 GenomeAMRAnalyzer v$Version Release Automation" -ForegroundColor Magenta
    Write-Host "==================================================" -ForegroundColor Magenta
    Write-Host ""
    
    # Confirmation prompt
    if (-not $Force) {
        Write-Warning "This will:"
        Write-Host "  • Build Docker image ${DOCKER_IMAGE}:${Version}"
        Write-Host "  • Test the image functionality"
        Write-Host "  • Publish to Docker Hub"
        Write-Host "  • Create and push Git tag v$Version"
        Write-Host ""
        $response = Read-Host "Do you want to proceed with the release? (y/n)"
        if ($response -ne "y" -and $response -ne "Y") {
            Write-Info "Release aborted by user"
            exit 0
        }
    }
    
    # Execute release steps
    Test-Prerequisites
    Test-RepositoryState
    Build-DockerImage
    Test-DockerImage
    Publish-DockerImage
    New-GitTag
    Show-ReleaseSummary
    
    Write-Success "🎉 GenomeAMRAnalyzer v$Version released successfully!"
}

# Run main function
Start-Release