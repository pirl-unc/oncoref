#!/bin/bash
set -o errexit

ruff check cancerdata tests
ruff format --check cancerdata tests
echo "All checks passed!"
