#!/usr/bin/env python
"""Verify pyproject.toml configuration is correct."""

import tomllib
from pathlib import Path

def verify_pyproject():
    """Verify that pyproject.toml has correct structure."""
    with open("pyproject.toml", "rb") as f:
        config = tomllib.load(f)
    
    # Test 1: No pytest in production dependencies
    prod_deps = config["project"]["dependencies"]
    for dep in prod_deps:
        assert "pytest" not in dep.lower(), f"pytest should not be in production: {dep}"
        assert "==" not in dep or dep.startswith("string=="), f"Use ranges, not ==: {dep}"
    
    # Test 2: pytest in dev dependencies
    dev_deps = config["project"]["optional-dependencies"]["dev"]
    assert any("pytest" in d for d in dev_deps), "pytest should be in dev deps"
    assert any("pytest-cov" in d for d in dev_deps), "pytest-cov should be in dev deps"
    
    # Test 3: Version ranges for production
    assert any("scipy" in d for d in prod_deps), "scipy dependency missing"
    assert any("numpy" in d for d in prod_deps), "numpy dependency missing"
    assert any("matplotlib" in d for d in prod_deps), "matplotlib dependency missing"
    
    print("All pyproject.toml checks passed!")
    print(f"\nProduction dependencies: {len(prod_deps)}")
    for dep in prod_deps:
        print(f"  * {dep}")
    print(f"\nDev dependencies: {len(dev_deps)}")
    for dep in dev_deps:
        print(f"  * {dep}")

if __name__ == "__main__":
    verify_pyproject()