#!/usr/bin/env bash

coverage run --source './bin' -m pytest -v /ampseq-pipeline/tests/python-tests/count_reads_per_region/test_count_reads_per_region.py \
&& coverage html -d ./htmlcov \
&& coverage report --fail-under 80

