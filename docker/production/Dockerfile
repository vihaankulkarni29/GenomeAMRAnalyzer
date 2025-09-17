# Production Dockerfile for ProductionWildTypeAligner
# Multi-stage build for optimized production deployment
# Base image: Python 3.11 slim for security and size optimization

# Stage 1: Base dependencies
FROM python:3.11-slim-bullseye as base

# Set environment variables for production
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    DEBIAN_FRONTEND=noninteractive

# Create non-root user for security
RUN groupadd --gid 1001 appuser && \
    useradd --uid 1001 --gid appuser --shell /bin/bash --create-home appuser

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    wget \
    git \
    && rm -rf /var/lib/apt/lists/*

# Stage 2: Dependencies installation
FROM base as dependencies

# Copy requirements files
COPY requirements/base.txt requirements/production.txt /tmp/

# Install Python dependencies
RUN pip install --no-cache-dir -r /tmp/production.txt

# Stage 3: EMBOSS installation (optional but recommended)
FROM dependencies as emboss

# Install EMBOSS suite for optimal alignment performance
RUN apt-get update && apt-get install -y --no-install-recommends \
    emboss \
    emboss-data \
    && rm -rf /var/lib/apt/lists/*

# Verify EMBOSS installation
RUN water -help || echo "EMBOSS WATER installed but may need configuration"

# Stage 4: Application setup
FROM emboss as application

# Set working directory
WORKDIR /app

# Copy application code
COPY src/ /app/src/
COPY config/ /app/config/
COPY MetaDataHarvester/ /app/MetaDataHarvester/

# Create necessary directories with proper permissions
RUN mkdir -p /app/data /app/logs /app/results /app/cache && \
    chown -R appuser:appuser /app

# Copy entrypoint script
COPY docker/production/entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh && chown appuser:appuser /entrypoint.sh

# Switch to non-root user
USER appuser

# Set Python path
ENV PYTHONPATH=/app/src:/app

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD python -c "from src.production_wildtype_aligner import ProductionWildTypeAligner; print('Health check passed')" || exit 1

# Expose port (if running as service)
EXPOSE 8000

# Default command
ENTRYPOINT ["/entrypoint.sh"]
CMD ["python", "src/production_wildtype_aligner.py", "--help"]