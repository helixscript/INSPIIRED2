#!/usr/bin/env bash

set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)
sampler="$repo_root/tools/memoryMonitor.sh"
using_fake_renderer=0
if ! command -v Rscript >/dev/null 2>&1 ||
   ! Rscript --vanilla -e 'quit(status = !all(vapply(c("data.table", "ggplot2"), requireNamespace, logical(1), quietly = TRUE)))' >/dev/null 2>&1; then
  export INSPIIRED2_RSCRIPT="$repo_root/tests/unitTests/memoryMonitoring/fakeRscript.sh"
  using_fake_renderer=1
fi
tmp_base=${TMPDIR:-/tmp}
[[ -d "$tmp_base" && -w "$tmp_base" ]] || tmp_base=$(dirname "${BASH_SOURCE[0]}")
test_root=$(mktemp -d "$tmp_base/memoryMonitor.XXXXXX")
if [[ ${KEEP_MEMORY_TEST_OUTPUT:-0} == 1 ]]; then
  trap 'printf "test output retained at %s\n" "$test_root"' EXIT
else
  trap 'rm -rf "$test_root"' EXIT
fi

fake_cgroup="$test_root/cgroup"
output_dir="$test_root/output with spaces"
control_dir="$output_dir/control"
mkdir -p "$fake_cgroup" "$control_dir"

printf '%s\n' 209715200 > "$fake_cgroup/memory.current"
printf '%s\n' 536870912 > "$fake_cgroup/memory.max"
printf '%s\n' 0 > "$fake_cgroup/memory.swap.current"
printf '%s\n' \
  'anon 104857600' \
  'file 52428800' \
  'shmem 10485760' \
  'inactive_file 20971520' > "$fake_cgroup/memory.stat"
printf '%s\n' \
  'low 0' \
  'high 0' \
  'max 0' \
  'oom 0' \
  'oom_kill 0' > "$fake_cgroup/memory.events"

trace="$output_dir/INSPIIRED2_memoryUsage.tsv"
svg="$output_dir/INSPIIRED2_memoryUsage.svg"
html="$output_dir/INSPIIRED2_memoryUsage.html"
summary="$output_dir/INSPIIRED2_memoryUsage_summary.tsv"
ready="$control_dir/test.ready"
stop="$control_dir/test.stop"
status_file="$control_dir/test.status"
done="$control_dir/test.done"
error_file="$control_dir/test.error"

INSPIIRED2_CGROUP_DIR="$fake_cgroup" bash "$sampler" \
  --trace "$trace" --svg "$svg" --html "$html" --summary "$summary" \
  --run-id run_test --invocation-id invocation_test --module buildFragments \
  --file-tag buildFragments --sample-seconds 0.2 --plot-seconds 1 \
  --parent-pid $$ --ready-file "$ready" --stop-file "$stop" \
  --status-file "$status_file" --done-file "$done" --error-file "$error_file" &
sampler_pid=$!

for _ in $(seq 1 100); do
  [[ -e "$ready" ]] && break
  sleep 0.05
done
[[ -e "$ready" ]]

printf '%s\n' 314572800 > "$fake_cgroup/memory.current"
sleep 1.3
[[ -s "$svg" ]]
((using_fake_renderer == 0)) || grep -q '0.27 GiB' "$svg"

# A competing sampler for the same output directory must fail softly without
# touching the shared trace. The biological module would continue unmonitored.
locked_ready="$control_dir/locked.ready"
locked_done="$control_dir/locked.done"
locked_error="$control_dir/locked.error"
INSPIIRED2_CGROUP_DIR="$fake_cgroup" bash "$sampler" \
  --trace "$trace" --svg "$svg" --html "$html" --summary "$summary" \
  --run-id competing_run --invocation-id competing_invocation --module prepReads \
  --file-tag prepReads --sample-seconds 0.2 --plot-seconds 1 --parent-pid $$ \
  --ready-file "$locked_ready" --stop-file "$control_dir/locked.stop" \
  --status-file "$control_dir/locked.status" --done-file "$locked_done" \
  --error-file "$locked_error"
[[ -e "$locked_done" && -s "$locked_error" && ! -e "$locked_ready" ]]
grep -q 'another memory monitor is already active' "$locked_error"
! grep -q 'competing_invocation' "$trace"

printf '%s\n' 7 > "$status_file"
: > "$stop"
wait "$sampler_pid"

[[ -s "$trace" && -s "$svg" && -s "$html" && -s "$summary" && -e "$done" ]]
[[ -z $(find "$output_dir" -type f -name '*.tmp_*' -print -quit) ]]
((using_fake_renderer == 0)) || grep -q 'buildFragments' "$svg"
grep -q 'This page refreshes every 1 seconds' "$html"
if command -v python3 >/dev/null 2>&1; then
  python3 -c 'import sys, xml.etree.ElementTree as ET; ET.parse(sys.argv[1])' "$svg"
fi

awk -F '\t' '
  NR == 1 { for(i = 1; i <= NF; i++) col[$i] = i; next }
  {
    rows++
    if($(col["module"]) != "buildFragments") exit 1
    if($(col["working_set_bytes"]) != $(col["memory_bytes"]) - $(col["inactive_file_bytes"])) exit 1
    if($(col["event"]) == "end" && $(col["exit_status"]) == 7) found_end = 1
    if($(col["memory_bytes"]) == 314572800) found_update = 1
  }
  END { if(rows < 3 || !found_end || !found_update) exit 1 }
' "$trace"

grep -q $'run_test\tinvocation_test\tbuildFragments' "$summary"
grep -q $'\t7$' "$summary"

# Private cgroup namespaces can expose 0::/ while mountinfo retains a
# host-side root. Discovery must then use the cgroup2 mountpoint itself.
discovery_dir="$test_root/discovery_cgroup"
discovery_proc="$test_root/discovery.proc.cgroup"
discovery_mountinfo="$test_root/discovery.proc.mountinfo"
mkdir -p "$discovery_dir"
cp "$fake_cgroup"/* "$discovery_dir/"
printf '%s\n' '0::/' > "$discovery_proc"
printf '29 23 0:26 /docker/example %s rw,nosuid,nodev,noexec,relatime - cgroup2 cgroup rw\n' \
  "$discovery_dir" > "$discovery_mountinfo"
discovery_control="$test_root/discovery-control"
mkdir -p "$discovery_control"
discovery_ready="$discovery_control/discovery.ready"
discovery_stop="$discovery_control/discovery.stop"
discovery_status="$discovery_control/discovery.status"
discovery_done="$discovery_control/discovery.done"
discovery_error="$discovery_control/discovery.error"
INSPIIRED2_TEST_PROC_CGROUP="$discovery_proc" \
INSPIIRED2_TEST_PROC_MOUNTINFO="$discovery_mountinfo" bash "$sampler" \
  --trace "$output_dir/discovery.tsv" --svg "$output_dir/discovery.svg" \
  --html "$output_dir/discovery.html" --summary "$output_dir/discovery_summary.tsv" \
  --run-id run_discovery --invocation-id invocation_discovery --module demultiplex \
  --file-tag demultiplex --sample-seconds 0.2 --plot-seconds 1 --parent-pid $$ \
  --ready-file "$discovery_ready" --stop-file "$discovery_stop" \
  --status-file "$discovery_status" --done-file "$discovery_done" \
  --error-file "$discovery_error" &
sampler_pid=$!
for _ in $(seq 1 100); do
  [[ -e "$discovery_ready" ]] && break
  sleep 0.05
done
[[ -e "$discovery_ready" && ! -e "$discovery_error" ]]
printf '%s\n' 0 > "$discovery_status"
: > "$discovery_stop"
wait "$sampler_pid"
grep -q $'\t2\t' "$output_dir/discovery.tsv"

# A second module invocation must append to the same run and update the plot.
ready2="$control_dir/test2.ready"
stop2="$control_dir/test2.stop"
status2="$control_dir/test2.status"
done2="$control_dir/test2.done"
error2="$control_dir/test2.error"
printf '%s\n' 419430400 > "$fake_cgroup/memory.current"
INSPIIRED2_CGROUP_DIR="$fake_cgroup" bash "$sampler" \
  --trace "$trace" --svg "$svg" --html "$html" --summary "$summary" \
  --run-id run_test --invocation-id invocation_test2 --module prepReads \
  --file-tag prepReads --sample-seconds 0.2 --plot-seconds 1 \
  --parent-pid $$ --ready-file "$ready2" --stop-file "$stop2" \
  --status-file "$status2" --done-file "$done2" --error-file "$error2" &
sampler_pid=$!
for _ in $(seq 1 100); do
  [[ -e "$ready2" ]] && break
  sleep 0.05
done
[[ -e "$ready2" ]]
printf '%s\n' 0 > "$status2"
: > "$stop2"
wait "$sampler_pid"

if ((using_fake_renderer == 1)); then
  grep -q 'buildFragments' "$svg"
  grep -q 'prepReads' "$svg"
fi
[[ $(grep -c '^run_id' "$trace") -eq 1 ]]
[[ $(grep -c '^run_test' "$summary") -eq 2 ]]

# A new run ID keeps old measurements in the TSV but excludes them from the live plot.
run2_ready="$control_dir/run2.ready"
run2_stop="$control_dir/run2.stop"
run2_status="$control_dir/run2.status"
run2_done="$control_dir/run2.done"
run2_error="$control_dir/run2.error"
INSPIIRED2_CGROUP_DIR="$fake_cgroup" bash "$sampler" \
  --trace "$trace" --svg "$svg" --html "$html" --summary "$summary" \
  --run-id run_test2 --invocation-id invocation_run2 --module demultiplex \
  --file-tag demultiplex --sample-seconds 0.2 --plot-seconds 1 \
  --parent-pid $$ --ready-file "$run2_ready" --stop-file "$run2_stop" \
  --status-file "$run2_status" --done-file "$run2_done" --error-file "$run2_error" &
sampler_pid=$!
for _ in $(seq 1 100); do
  [[ -e "$run2_ready" ]] && break
  sleep 0.05
done
[[ -e "$run2_ready" ]]
printf '%s\n' 0 > "$run2_status"
: > "$run2_stop"
wait "$sampler_pid"
if ((using_fake_renderer == 1)); then
  grep -q 'demultiplex' "$svg"
  ! grep -q 'buildFragments' "$svg"
  ! grep -q 'prepReads' "$svg"
fi
[[ $(grep -c '^run_test2' "$summary") -eq 1 ]]

# The cgroup-v1 fallback must normalize its unlimited sentinel and mem+swap value.
fake_v1="$test_root/cgroup_v1"
mkdir -p "$fake_v1"
printf '%s\n' 104857600 > "$fake_v1/memory.usage_in_bytes"
printf '%s\n' 9223372036854771712 > "$fake_v1/memory.limit_in_bytes"
printf '%s\n' 125829120 > "$fake_v1/memory.memsw.usage_in_bytes"
printf '%s\n' 2 > "$fake_v1/memory.failcnt"
printf '%s\n' \
  'total_rss 73400320' \
  'total_cache 31457280' \
  'total_shmem 5242880' \
  'total_inactive_file 10485760' > "$fake_v1/memory.stat"

v1_trace="$output_dir/v1.tsv"
v1_ready="$control_dir/v1.ready"
v1_stop="$control_dir/v1.stop"
v1_status="$control_dir/v1.status"
v1_done="$control_dir/v1.done"
v1_error="$control_dir/v1.error"
INSPIIRED2_CGROUP_DIR="$fake_v1" bash "$sampler" \
  --trace "$v1_trace" --svg "$output_dir/v1.svg" --html "$output_dir/v1.html" \
  --summary "$output_dir/v1_summary.tsv" --run-id run_v1 --invocation-id invocation_v1 \
  --module alignReads --file-tag alignReads --sample-seconds 0.2 --plot-seconds 1 \
  --parent-pid $$ --ready-file "$v1_ready" --stop-file "$v1_stop" \
  --status-file "$v1_status" --done-file "$v1_done" --error-file "$v1_error" &
sampler_pid=$!
for _ in $(seq 1 100); do
  [[ -e "$v1_ready" ]] && break
  sleep 0.05
done
[[ -e "$v1_ready" ]]
printf '%s\n' 0 > "$v1_status"
: > "$v1_stop"
wait "$sampler_pid"

awk -F '\t' '
  NR == 1 { for(i = 1; i <= NF; i++) col[$i] = i; next }
  NR == 2 {
    if($(col["cgroup_version"]) != 1) exit 1
    if($(col["memory_limit_bytes"]) != "") exit 1
    if($(col["swap_bytes"]) != 20971520) exit 1
    if($(col["oom_events"]) != 2) exit 1
  }
' "$v1_trace"

# A report-rendering failure is diagnostic only: sampling still reaches its
# end row and records a nonfatal warning for the launcher.
render_ready="$control_dir/render.ready"
render_stop="$control_dir/render.stop"
render_status="$control_dir/render.status"
render_done="$control_dir/render.done"
render_error="$control_dir/render.error"
render_trace="$output_dir/render_failure.tsv"
INSPIIRED2_CGROUP_DIR="$fake_cgroup" bash "$sampler" \
  --trace "$render_trace" --svg "/proc/INSPIIRED2_unwritable.svg" \
  --html "$output_dir/render_failure.html" --summary "$output_dir/render_failure_summary.tsv" \
  --run-id run_render --invocation-id invocation_render --module prepReads \
  --file-tag prepReads --sample-seconds 0.2 --plot-seconds 1 --parent-pid $$ \
  --ready-file "$render_ready" --stop-file "$render_stop" \
  --status-file "$render_status" --done-file "$render_done" --error-file "$render_error" \
  2> "$output_dir/render_failure.stderr" &
sampler_pid=$!
for _ in $(seq 1 100); do
  [[ -e "$render_ready" ]] && break
  sleep 0.05
done
[[ -e "$render_ready" ]]
printf '%s\n' 0 > "$render_status"
: > "$render_stop"
wait "$sampler_pid"
[[ -e "$render_done" && -s "$render_error" ]]
grep -q 'SVG plot' "$render_error"
awk -F '\t' '$8 == "end" && $9 == 0 { found = 1 } END { exit !found }' "$render_trace"

missing="$test_root/missing"
missing_done="$control_dir/missing.done"
missing_error="$control_dir/missing.error"
INSPIIRED2_CGROUP_DIR="$missing" bash "$sampler" \
  --trace "$output_dir/missing.tsv" --svg "$output_dir/missing.svg" \
  --html "$output_dir/missing.html" --summary "$output_dir/missing_summary.tsv" \
  --run-id run_missing --invocation-id invocation_missing --module prepReads \
  --file-tag prepReads --sample-seconds 1 --plot-seconds 1 --parent-pid $$ \
  --ready-file "$control_dir/missing.ready" --stop-file "$control_dir/missing.stop" \
  --status-file "$control_dir/missing.status" --done-file "$missing_done" \
  --error-file "$missing_error"

[[ -e "$missing_done" && -s "$missing_error" ]]
grep -q 'could not be identified' "$missing_error"

printf '%s\n' 'memory monitoring tests passed'
