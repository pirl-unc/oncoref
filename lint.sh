#!/bin/bash
set -o errexit

ruff check oncoref tests scripts/*cta*.py scripts/audit_other_*.py scripts/build_other_literature_report.py
ruff format --check oncoref tests scripts/*cta*.py scripts/audit_other_*.py scripts/build_other_literature_report.py
echo "All checks passed!"
