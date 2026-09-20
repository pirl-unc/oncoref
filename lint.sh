#!/bin/bash
set -o errexit

ruff check oncoref tests scripts/*cta*.py
ruff format --check oncoref tests scripts/*cta*.py
echo "All checks passed!"
