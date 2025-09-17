#!/bin/bash
# Production Deployment Script for ProductionWildTypeAligner
# Automated deployment with validation and rollback capabilities

set -e

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
COMPOSE_FILE="$PROJECT_ROOT/docker/production/docker-compose.yml"
BACKUP_DIR="/var/backups/wildtype-aligner"
LOG_FILE="/var/log/wildtype-aligner-deployment.log"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Logging function
log() {
    echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# Error handling
error_exit() {
    log "${RED}ERROR: $1${NC}"
    exit 1
}

# Success message
success() {
    log "${GREEN}SUCCESS: $1${NC}"
}

# Warning message
warning() {
    log "${YELLOW}WARNING: $1${NC}"
}

# Check prerequisites
check_prerequisites() {
    log "Checking deployment prerequisites..."
    
    # Check if running as root or with sufficient privileges
    if [[ $EUID -ne 0 ]] && ! groups $USER | grep -q docker; then
        error_exit "Script must be run as root or user must be in docker group"
    fi
    
    # Check Docker installation
    if ! command -v docker >/dev/null 2>&1; then
        error_exit "Docker is not installed"
    fi
    
    # Check Docker Compose installation
    if ! command -v docker-compose >/dev/null 2>&1; then
        error_exit "Docker Compose is not installed"
    fi
    
    # Check if compose file exists
    if [[ ! -f "$COMPOSE_FILE" ]]; then
        error_exit "Docker Compose file not found: $COMPOSE_FILE"
    fi
    
    # Check disk space (minimum 10GB)
    available_space=$(df / | awk 'NR==2 {print $4}')
    if [[ $available_space -lt 10485760 ]]; then
        error_exit "Insufficient disk space. At least 10GB required."
    fi
    
    success "All prerequisites met"
}

# Backup existing deployment
backup_existing() {
    log "Creating backup of existing deployment..."
    
    # Create backup directory
    mkdir -p "$BACKUP_DIR/$(date +%Y%m%d_%H%M%S)"
    CURRENT_BACKUP="$BACKUP_DIR/$(date +%Y%m%d_%H%M%S)"
    
    # Backup volumes if they exist
    if docker volume ls | grep -q wildtype-aligner; then
        log "Backing up Docker volumes..."
        docker run --rm -v wildtype-aligner-redis-data:/data -v "$CURRENT_BACKUP":/backup alpine tar czf /backup/redis-data.tar.gz -C /data .
        docker run --rm -v wildtype-aligner-postgres-data:/data -v "$CURRENT_BACKUP":/backup alpine tar czf /backup/postgres-data.tar.gz -C /data .
    fi
    
    # Backup configuration
    if [[ -d "$PROJECT_ROOT/config/production" ]]; then
        cp -r "$PROJECT_ROOT/config/production" "$CURRENT_BACKUP/"
    fi
    
    success "Backup created at $CURRENT_BACKUP"
}

# Build Docker images
build_images() {
    log "Building Docker images..."
    
    cd "$PROJECT_ROOT"
    
    # Build production image
    docker-compose -f "$COMPOSE_FILE" build --no-cache wildtype-aligner
    
    # Tag with timestamp for rollback capability
    docker tag genome-amr-analyzer/wildtype-aligner:production \
               genome-amr-analyzer/wildtype-aligner:production-$(date +%Y%m%d_%H%M%S)
    
    success "Docker images built successfully"
}

# Deploy services
deploy_services() {
    log "Deploying services..."
    
    cd "$PROJECT_ROOT"
    
    # Start services
    docker-compose -f "$COMPOSE_FILE" up -d
    
    # Wait for services to be healthy
    log "Waiting for services to become healthy..."
    timeout=300
    elapsed=0
    
    while [[ $elapsed -lt $timeout ]]; do
        if docker-compose -f "$COMPOSE_FILE" ps | grep -q "healthy"; then
            success "Services are healthy"
            return 0
        fi
        sleep 10
        elapsed=$((elapsed + 10))
        log "Waiting for health checks... ($elapsed/$timeout seconds)"
    done
    
    error_exit "Services failed to become healthy within $timeout seconds"
}

# Validate deployment
validate_deployment() {
    log "Validating deployment..."
    
    # Check if containers are running
    if ! docker-compose -f "$COMPOSE_FILE" ps | grep -q "Up"; then
        error_exit "Containers are not running properly"
    fi
    
    # Test application functionality
    log "Testing application functionality..."
    
    # Run health check command inside container
    if docker-compose -f "$COMPOSE_FILE" exec -T wildtype-aligner python -c "
from src.production_wildtype_aligner import ProductionWildTypeAligner
from src.sepi_configuration_manager import SEPIConfigurationManager
print('Application validation successful')
"; then
        success "Application validation passed"
    else
        error_exit "Application validation failed"
    fi
    
    # Check resource usage
    log "Checking resource usage..."
    docker stats --no-stream --format "table {{.Container}}\t{{.CPUPerc}}\t{{.MemUsage}}" \
        $(docker-compose -f "$COMPOSE_FILE" ps -q)
    
    success "Deployment validation completed"
}

# Cleanup old images and containers
cleanup() {
    log "Performing cleanup..."
    
    # Remove old images (keep last 3 versions)
    old_images=$(docker images genome-amr-analyzer/wildtype-aligner --format "{{.Tag}}" | \
                grep "production-" | sort -r | tail -n +4)
    
    if [[ -n "$old_images" ]]; then
        echo "$old_images" | while read tag; do
            docker rmi "genome-amr-analyzer/wildtype-aligner:$tag" || true
        done
        log "Cleaned up old images"
    fi
    
    # Remove unused volumes and networks
    docker system prune -f --volumes
    
    success "Cleanup completed"
}

# Rollback function
rollback() {
    log "Performing rollback..."
    
    # Stop current deployment
    docker-compose -f "$COMPOSE_FILE" down
    
    # Get latest backup
    latest_backup=$(ls -1t "$BACKUP_DIR" | head -1)
    
    if [[ -n "$latest_backup" ]]; then
        log "Restoring from backup: $latest_backup"
        
        # Restore volumes
        if [[ -f "$BACKUP_DIR/$latest_backup/redis-data.tar.gz" ]]; then
            docker run --rm -v wildtype-aligner-redis-data:/data -v "$BACKUP_DIR/$latest_backup":/backup alpine \
                sh -c "cd /data && tar xzf /backup/redis-data.tar.gz"
        fi
        
        if [[ -f "$BACKUP_DIR/$latest_backup/postgres-data.tar.gz" ]]; then
            docker run --rm -v wildtype-aligner-postgres-data:/data -v "$BACKUP_DIR/$latest_backup":/backup alpine \
                sh -c "cd /data && tar xzf /backup/postgres-data.tar.gz"
        fi
        
        success "Rollback completed"
    else
        error_exit "No backup found for rollback"
    fi
}

# Main deployment function
main() {
    log "Starting ProductionWildTypeAligner deployment"
    
    # Parse command line arguments
    case "${1:-deploy}" in
        "deploy")
            check_prerequisites
            backup_existing
            build_images
            deploy_services
            validate_deployment
            cleanup
            success "Deployment completed successfully"
            ;;
        "rollback")
            rollback
            ;;
        "validate")
            validate_deployment
            ;;
        "cleanup")
            cleanup
            ;;
        *)
            echo "Usage: $0 {deploy|rollback|validate|cleanup}"
            exit 1
            ;;
    esac
}

# Trap signals for graceful shutdown
trap 'error_exit "Deployment interrupted"' INT TERM

# Run main function
main "$@"