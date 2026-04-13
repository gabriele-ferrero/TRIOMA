# Pre-commit Setup Guide

This document explains how to set up and use pre-commit hooks for TRIOMA.

## Installation

### 1. Install pre-commit and tools
```bash
pip install -e ".[dev]"
```

### 2. Install git hooks
```bash
pre-commit install
```

This creates `.git/hooks/pre-commit` that runs before each commit.

## Usage

### Automatic (on every commit)
Once installed, pre-commit hooks run automatically before each commit:
```bash
git commit -m "Your commit message"
# Hooks run automatically, fix issues if needed
# If files were modified, stage them and commit again
```

### Manual (on demand)
Run all hooks on all files:
```bash
pre-commit run --all-files
```

Run specific hook:
```bash
pre-commit run black --all-files
pre-commit run ruff --all-files
```

## Hooks Included

| Hook | Purpose | Auto-fix |
|------|---------|----------|
| `trailing-whitespace` | Remove trailing spaces | ✓ Yes |
| `end-of-file-fixer` | Fix missing newline at EOF | ✓ Yes |
| `check-yaml` | Validate YAML syntax | ✗ No |
| `check-json` | Validate JSON syntax | ✗ No |
| `check-toml` | Validate TOML syntax | ✗ No |
| `check-added-large-files` | Prevent large files (>500KB) | ✗ No |
| `check-merge-conflict` | Detect merge conflicts | ✗ No |
| `detect-private-key` | Detect hardcoded secrets | ✗ No |
| `black` | Format Python code | ✓ Yes |
| `ruff` | Lint and fix Python code | ✓ Yes |
| `mypy` | Type checking | ✗ No |
| `pydocstyle` | Check docstring style | ✗ No |
| `bandit` | Security scanning | ✗ No |
| `autopep8` | Python syntax check | ✓ Yes (with args) |
| `markdownlint` | Lint Markdown files | ✓ Yes |

## Workflow

1. **Make changes** to your code
2. **Stage files**: `git add .`
3. **Commit**: `git commit -m "Your message"`
4. **Hooks run**:
   - If all pass → commit succeeds
   - If some fail with auto-fix → files are modified, re-stage and commit
   - If some fail without auto-fix → fix manually, re-stage, and commit

Example:
```bash
# Make changes
echo "print('hello')" > test.py

# Try to commit
git add test.py
git commit -m "add test"

# Black reformats it automatically
# Ruff may suggest fixes
# Need to stage changes and commit again
git add test.py
git commit -m "add test"
# Now it passes
```

## Skipping Hooks (not recommended)

If absolutely necessary, skip hooks with:
```bash
git commit --no-verify
# or
SKIP=hook-id git commit -m "message"
```

## Troubleshooting

### Hooks not running
```bash
# Reinstall hooks
pre-commit install --install-hooks
```

### Update hooks to latest versions
```bash
pre-commit autoupdate
```

### Debugging a specific hook
```bash
pre-commit run <hook-id> --all-files --verbose
```

## Configuration

Edit `.pre-commit-config.yaml` to:
- Add/remove hooks
- Change hook versions
- Modify hook arguments
- Exclude directories

Example: to skip docs folder, edit:
```yaml
exclude: |
  (?x)^(
    docs/|
    Examples/|
    build/
  )
```

## Disabling Pre-commit

If you need to disable pre-commit:
```bash
pre-commit uninstall
```

To re-enable:
```bash
pre-commit install
```
