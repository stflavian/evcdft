# GitHub Actions Workflows for EVC_DFT

This directory contains GitHub Actions workflow files for automated testing and CI/CD.

## Available Workflows

### 1. `test.yml` - Comprehensive Test Suite

**Trigger**: Push or pull request to `main` branch

**Purpose**: Runs the complete test suite across multiple Julia versions

**Features**:
- Tests on Julia 1.9 and 1.10
- Caches Julia packages for faster runs
- Runs all unit tests, integration tests, and phase 1 tests
- Runs basic validation tests

**Usage**: Automatically runs on every push and pull request to main

---

### 2. `ci.yml` - Full CI Pipeline

**Trigger**: Push or pull request to `main` branch

**Purpose**: Comprehensive CI pipeline testing across multiple operating systems

**Features**:
- Tests on Julia 1.9 and 1.10
- Runs on Ubuntu, macOS, and Windows
- Caches Julia packages
- Precompiles packages for faster runs
- Runs comprehensive test suite

**Usage**: Automatically runs on every push and pull request to main

---

### 3. `quick-test.yml` - Quick Validation

**Trigger**: Push to `main`, `dev`, or `feature/*` branches, and pull requests to `main`

**Purpose**: Quick validation of basic functionality

**Features**:
- Runs on Julia 1.9 (single version)
- Only runs basic validation tests
- Faster than full test suite
- Good for quick feedback during development

**Usage**: Automatically runs on every push to development branches

---

## How to Use

### Running Workflows Manually

1. Go to the **Actions** tab in your GitHub repository
2. Select the workflow you want to run
3. Click **Run workflow** button
4. Select the branch and click **Run workflow**

### Viewing Results

1. Go to the **Actions** tab in your GitHub repository
2. Click on the workflow run you want to view
3. Expand the steps to see detailed output
4. Check for any failures or errors

### Debugging Failures

If a workflow fails:

1. Check the error message in the workflow output
2. Look at which test failed
3. Run the failing test locally to reproduce the issue
4. Fix the issue and push a new commit

## Customizing Workflows

### Adding Julia Versions

To test on additional Julia versions, add them to the `matrix.version` list in the workflow file:

```yaml
matrix:
  version:
    - '1.9'
    - '1.10'
    - '1.11'  # Add new version
```

### Adding Test Files

To add new test files to the test suite, include them in the workflow:

```yaml
include("test/unit_tests.jl")
include("test/integration_tests.jl")
include("test/test_phase1.jl")
include("test/new_tests.jl")  # Add new test file
```

### Changing Trigger Events

To change when workflows run, modify the `on` section:

```yaml
on:
  push:
    branches: [ main, dev ]  # Add dev branch
  pull_request:
    branches: [ main ]
  schedule:
    - cron: '0 0 * * *'  # Run daily at midnight
```

## Best Practices

1. **Keep workflows fast**: Quick validation workflows should complete in under 5 minutes
2. **Use caching**: Cache Julia packages to speed up workflows
3. **Test on multiple versions**: Ensure compatibility across Julia versions
4. **Test on multiple OS**: Catch platform-specific issues early
5. **Separate concerns**: Use different workflows for different purposes (quick validation vs. comprehensive testing)

## Environment Variables

The workflows use the following environment variables:

- `JULIA_NUM_THREADS: auto` - Uses all available threads
- `JULIA_PROJECT: @.` - Uses the project in the current directory

## Dependencies

The workflows require:

- Julia (specified in the matrix)
- FFTW.jl package (installed automatically)
- Test.jl package (part of Julia standard library)

## Troubleshooting

### "Package not found" errors

Make sure:
1. The package is listed in `Project.toml`
2. The package is installed with `Pkg.add()` or `Pkg.instantiate()`
3. The correct Julia version is being used

### "Module not found" errors

Make sure:
1. The module is properly exported in the main module file
2. The module file is included in the main module
3. The correct `using` statement is used in tests

### Slow workflows

To speed up workflows:
1. Use caching for Julia packages
2. Reduce the number of Julia versions tested
3. Split tests into separate workflows
4. Use `Pkg.precompile()` to precompile packages
