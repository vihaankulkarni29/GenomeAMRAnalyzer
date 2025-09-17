#!/bin/bash
# Production entrypoint script for ProductionWildTypeAligner
# Handles initialization, validation, and graceful startup

set -e

# Function to log messages with timestamp
log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

# Function to validate environment
validate_environment() {
    log "Validating production environment..."
    
    # Check Python imports
    python -c "
import sys
sys.path.insert(0, '/app/src')
try:
    from production_wildtype_aligner import ProductionWildTypeAligner
    from sepi_configuration_manager import SEPIConfigurationManager
    print('✅ Core imports successful')
except ImportError as e:
    print(f'❌ Import failed: {e}')
    sys.exit(1)
"
    
    # Check required directories
    for dir in /app/data /app/logs /app/results /app/cache; do
        if [ ! -d "$dir" ]; then
            log "Creating missing directory: $dir"
            mkdir -p "$dir"
        fi
    done
    
    # Check EMBOSS installation
    if command -v water >/dev/null 2>&1; then
        log "✅ EMBOSS WATER found"
    else
        log "⚠️  EMBOSS WATER not found - will use BioPython fallback"
    fi
    
    log "Environment validation completed"
}

# Function to setup configuration
setup_configuration() {
    log "Setting up production configuration..."
    
    # Set default environment variables if not provided
    export PYTHONPATH="${PYTHONPATH:-/app/src:/app}"
    export ALIGNER_LOG_LEVEL="${ALIGNER_LOG_LEVEL:-INFO}"
    export ALIGNER_MAX_CONCURRENT="${ALIGNER_MAX_CONCURRENT:-4}"
    export ALIGNER_OUTPUT_DIR="${ALIGNER_OUTPUT_DIR:-/app/results}"
    
    # Create configuration directory if it doesn't exist
    mkdir -p /app/config/production
    
    log "Configuration setup completed"
}

# Function to handle graceful shutdown
cleanup() {
    log "Received shutdown signal, performing cleanup..."
    # Add any cleanup operations here
    log "Cleanup completed"
    exit 0
}

# Trap signals for graceful shutdown
trap cleanup SIGTERM SIGINT

# Main execution
main() {
    log "Starting ProductionWildTypeAligner container"
    
    # Validate environment
    validate_environment
    
    # Setup configuration
    setup_configuration
    
    # If no arguments provided, show help
    if [ $# -eq 0 ]; then
        log "No arguments provided, showing help"
        exec python /app/src/production_wildtype_aligner.py --help
    fi
    
    # Execute the provided command
    log "Executing command: $*"
    exec "$@"
}

# Run main function
main "$@"