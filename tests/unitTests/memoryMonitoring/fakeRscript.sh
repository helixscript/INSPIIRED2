#!/usr/bin/env bash

# Minimal renderer used only when the unit-test host lacks R/ggplot2. The real
# renderer is exercised automatically when those dependencies are available.
set -euo pipefail

[[ ${1:-} == "--vanilla" ]] && shift
[[ $# -eq 4 ]]
trace_file=$2
output_svg=$3
run_id=$4

awk -F '\t' -v target_run="$run_id" '
  NR == 1 {
    for(i = 1; i <= NF; i++) col[$i] = i
    next
  }
  $(col["run_id"]) == target_run {
    modules[$(col["module"])] = 1
    working = $(col["working_set_bytes"]) + 0
    if(working > peak) peak = working
  }
  END {
    print "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
    print "<svg xmlns=\"http://www.w3.org/2000/svg\">"
    printf "<text>%s</text>\n", target_run
    printf "<text>%.2f GiB</text>\n", peak / 1073741824
    for(module in modules) printf "<text>%s</text>\n", module
    print "</svg>"
  }
' "$trace_file" > "$output_svg"
