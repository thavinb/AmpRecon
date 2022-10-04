#!/usr/bin/env bash

coverage run --source './bin' -m pytest -v python_tests/count_reads_per_region/test_*.py \
&& coverage html -d ./htmlcov \
&& coverage report --fail-under 85

