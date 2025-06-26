# Contributing to CoinfectionSimulator.jl

Thank you for your interest in contributing to CoinfectionSimulator.jl! We welcome contributions from the community.

## Ways to Contribute

- Report bugs or issues
- Suggest new features or improvements
- Submit code improvements or bug fixes
- Improve documentation
- Add examples or tutorials

## Getting Started

1. Fork the repository on GitHub
2. Clone your fork locally:
   ```bash
   git clone https://github.com/yourusername/CoinfectionSimulator.jl.git
   ```
3. Create a new branch for your feature:
   ```bash
   git checkout -b feature/your-feature-name
   ```

## Development Setup

1. Navigate to the package directory and activate the environment:
   ```julia
   using Pkg
   Pkg.activate(".")
   Pkg.instantiate()
   ```

2. Run the tests to ensure everything works:
   ```julia
   Pkg.test()
   ```

## Code Style

- Follow standard Julia coding conventions
- Use descriptive variable and function names
- Add docstrings to all public functions
- Include type annotations where helpful
- Write tests for new functionality

## Testing

- All new code should include appropriate tests
- Ensure all existing tests continue to pass
- Aim for comprehensive test coverage
- Tests should be fast and deterministic when possible

## Documentation

- Update docstrings for any modified functions
- Add examples to demonstrate new features
- Update the README.md if needed
- Consider adding to the examples/ directory

## Pull Request Process

1. Ensure your code follows the style guidelines
2. Add or update tests as needed
3. Update documentation
4. Commit your changes with clear, descriptive messages
5. Push to your fork and submit a pull request
6. Respond to any feedback during the review process

## Questions?

If you have questions about contributing, please open an issue or reach out to the maintainers.

Thank you for contributing!
