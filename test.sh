#!/bin/bash
set -o errexit

./lint.sh
python -m pytest tests/ --cov=cancerdata --cov-report=term-missing
