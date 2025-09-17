# Contributing to GenomeAMRAnalyzer

Thank you for your interest in contributing to GenomeAMRAnalyzer! This document provides guidelines for contributing to the project.

## 🤝 Code of Conduct

We are committed to providing a welcoming and inclusive environment for all contributors. Please read and follow our [Code of Conduct](CODE_OF_CONDUCT.md).

## 🚀 Getting Started

### Development Environment Setup

1. **Fork and Clone**
   ```bash
   git clone https://github.com/yourusername/GenomeAMRAnalyzer.git
   cd GenomeAMRAnalyzer
   ```

2. **Create Virtual Environment**
   ```bash
   python -m venv dev-env
   source dev-env/bin/activate  # Linux/macOS
   # dev-env\Scripts\activate  # Windows
   ```

3. **Install Development Dependencies**
   ```bash
   pip install -e ".[dev,visualization,performance]"
   ```

4. **Verify Setup**
   ```bash
   pytest tests/ -v
   black --check src/
   flake8 src/
   ```

### Development Workflow

1. **Create Feature Branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make Changes**
   - Follow coding standards (see below)
   - Add tests for new functionality
   - Update documentation as needed

3. **Run Quality Checks**
   ```bash
   # Code formatting
   black src/ tests/
   
   # Linting
   flake8 src/ tests/
   
   # Type checking
   mypy src/
   
   # Run tests
   pytest tests/ -v --cov=src/
   ```

4. **Commit Changes**
   ```bash
   git add .
   git commit -m "feat: Add new feature description"
   ```

5. **Push and Create PR**
   ```bash
   git push origin feature/your-feature-name
   ```

## 📝 Contribution Types

### 🐛 Bug Reports

When reporting bugs, please include:

- **Environment Details**: OS, Python version, package versions
- **Reproduction Steps**: Minimal example to reproduce the issue
- **Expected vs Actual Behavior**: Clear description of the problem
- **Error Messages**: Complete error traces and log outputs
- **Configuration**: Relevant configuration files or settings

**Template:**
```markdown
## Bug Description
Brief description of the issue

## Environment
- OS: [e.g., Windows 10, Ubuntu 20.04]
- Python: [e.g., 3.11.2]
- GenomeAMRAnalyzer: [e.g., 2.0.0]

## Reproduction Steps
1. Run command: `...`
2. With configuration: `...`
3. Expected: `...`
4. Actual: `...`

## Error Output
```
[paste error message here]
```

## Additional Context
Any other relevant information
```

### ✨ Feature Requests

For new features, please provide:

- **Use Case**: Scientific or technical justification
- **Proposed Implementation**: High-level design approach
- **Alternative Solutions**: Other approaches considered
- **Impact Assessment**: Performance, compatibility considerations

### 🔧 Code Contributions

#### Coding Standards

**Python Style**
- Follow [PEP 8](https://pep8.org/) style guidelines
- Use [Black](https://black.readthedocs.io/) for code formatting
- Maximum line length: 88 characters (Black default)
- Use type hints for all function signatures

**Documentation**
- Docstrings for all public functions/classes
- Follow [Google docstring format](https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings)
- Include examples in docstrings for complex functions
- Update README.md for user-facing changes

**Testing**
- Write tests for all new functionality
- Maintain >90% code coverage
- Use pytest fixtures for test data
- Include integration tests for major features

**Error Handling**
- Use specific exception types
- Provide informative error messages
- Include context in error logs
- Implement graceful degradation where possible

#### Code Review Process

1. **Automated Checks**: All PRs must pass CI/CD pipeline
2. **Peer Review**: At least one reviewer approval required
3. **Documentation Review**: Ensure documentation is updated
4. **Testing Validation**: Verify tests cover new functionality

## 🧪 Testing Guidelines

### Test Categories

**Unit Tests**
- Test individual functions/methods
- Mock external dependencies
- Fast execution (<1 second per test)

**Integration Tests**
- Test module interactions
- Use real data when possible
- Validate end-to-end workflows

**Performance Tests**
- Benchmark critical functions
- Memory usage validation
- Scalability testing

### Test Data

- Use synthetic data for unit tests
- Provide small real datasets for integration tests
- Document data sources and licensing
- Avoid large files in repository

### Running Tests

```bash
# Run all tests
pytest tests/ -v

# Run with coverage
pytest tests/ -v --cov=src/ --cov-report=html

# Run specific test categories
pytest tests/unit/ -v
pytest tests/integration/ -v

# Run performance tests
pytest tests/performance/ -v --benchmark-only
```

## 📚 Documentation

### Types of Documentation

**User Documentation**
- README.md: Installation and basic usage
- User guides: Detailed tutorials
- API reference: Module documentation
- Examples: Common use cases

**Developer Documentation**
- Architecture overview
- Module design documents
- Testing strategies
- Deployment guides

### Documentation Standards

- Use Markdown for all documentation
- Include code examples with expected outputs
- Provide links to external resources
- Keep documentation current with code changes

## 🚀 Release Process

### Version Numbering

We follow [Semantic Versioning](https://semver.org/):
- **MAJOR**: Incompatible API changes
- **MINOR**: New functionality, backward compatible
- **PATCH**: Bug fixes, backward compatible

### Release Checklist

1. **Pre-release**
   - Update version numbers
   - Update CHANGELOG.md
   - Run full test suite
   - Update documentation

2. **Release**
   - Create release branch
   - Tag release version
   - Build and test packages
   - Create GitHub release

3. **Post-release**
   - Announce release
   - Update installation docs
   - Monitor for issues

## 🏗️ Architecture Guidelines

### Module Design

**Principles**
- Single responsibility principle
- Clear, documented interfaces
- Minimal dependencies between modules
- Consistent error handling patterns

**Structure**
```
src/
├── core/                 # Core functionality
├── analysis/            # Analysis modules
├── visualization/       # Plotting and reports
├── utils/              # Utility functions
└── tests/              # Test suites
```

### Configuration Management

- Use YAML for configuration files
- Validate configuration on startup
- Provide sensible defaults
- Document all configuration options

### Performance Considerations

- Profile code before optimization
- Use appropriate data structures
- Implement lazy loading where beneficial
- Consider memory usage for large datasets

## 🎯 Specific Contribution Areas

### High-Priority Areas

1. **Performance Optimization**
   - Large genome dataset handling
   - Memory usage improvements
   - Parallel processing enhancements

2. **Testing Expansion**
   - Edge case coverage
   - Integration test scenarios
   - Performance benchmarks

3. **Documentation Improvements**
   - Tutorial videos
   - Advanced usage examples
   - Troubleshooting guides

4. **New Analysis Features**
   - Additional resistance mechanisms
   - Statistical methods
   - Visualization options

### Beginner-Friendly Issues

Look for issues labeled:
- `good first issue`
- `documentation`
- `help wanted`
- `beginner`

## 🔍 Review Criteria

### Code Quality
- [ ] Follows coding standards
- [ ] Includes appropriate tests
- [ ] Has clear documentation
- [ ] Handles errors gracefully

### Functionality
- [ ] Solves stated problem
- [ ] Doesn't break existing features
- [ ] Performance is acceptable
- [ ] Integrates well with existing code

### Documentation
- [ ] Updated relevant documentation
- [ ] Includes usage examples
- [ ] Clear commit messages
- [ ] Updated CHANGELOG if needed

## 📞 Getting Help

### Communication Channels

- **GitHub Issues**: Bug reports and feature requests
- **GitHub Discussions**: General questions and community support
- **Email**: [maintainer@genomeamranalyzer.org] for sensitive issues

### Development Support

- **Architecture Questions**: Create GitHub discussion
- **Code Review**: Submit draft PR for early feedback
- **Design Decisions**: Open issue with `design` label

## 🏆 Recognition

### Contributors

All contributors are recognized in:
- Contributors section of README.md
- CONTRIBUTORS.md file
- Release notes for significant contributions

### Contribution Types

We recognize various contribution types:
- Code contributions
- Documentation improvements
- Bug reports and testing
- Community support and outreach
- Scientific validation and feedback

## 📄 Legal

### License Agreement

By contributing to GenomeAMRAnalyzer, you agree that your contributions will be licensed under the [MIT License](LICENSE).

### Copyright

- Retain copyright on your contributions
- Grant project license to use contributions
- Ensure you have right to license contributed code

---

Thank you for contributing to GenomeAMRAnalyzer! Your efforts help advance antimicrobial resistance research and improve global health outcomes.

For questions about contributing, please open a GitHub discussion or contact the maintainers.