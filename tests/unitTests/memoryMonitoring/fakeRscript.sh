#!/usr/bin/env bash

# Minimal PNG/PDF renderer used only when the unit-test host lacks R/ggplot2.
# The real renderer is exercised automatically when those dependencies exist.
set -euo pipefail

[[ ${1:-} == "--vanilla" ]] && shift
[[ $# -eq 6 ]]
trace_file=$2
summary_file=$3
output_png=$4
output_pdf=$5
run_id=$6

[[ -s "$trace_file" && -s "$summary_file" && -n "$run_id" ]]

# Valid 1 x 1 PNG fixture.
printf '%s' 'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk+A8AAQUBAScY42YAAAAASUVORK5CYII=' |
  base64 --decode > "$output_png"

# The fake PDF only supplies the signature/trailer needed by shell-only tests.
# Hosts with R exercise and validate the real multi-page PDF instead.
printf '%%PDF-1.4\n%% fake INSPIIRED2 memory report for %s\n%%%%EOF\n' "$run_id" > "$output_pdf"
