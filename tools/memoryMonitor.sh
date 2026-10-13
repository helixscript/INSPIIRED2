#!/usr/bin/env bash

# Lightweight container-wide memory sampler used by the INSPIIRED2 launcher.
# Report updates use a short-lived R/ggplot2 process between memory samples.

set -u

trace_file=""
png_file=""
pdf_file=""
summary_file=""
run_id=""
invocation_id=""
module=""
file_tag=""
sample_seconds="2"
plot_seconds="30"
parent_pid=""
ready_file=""
stop_file=""
status_file=""
done_file=""
error_file=""

while (($#)); do
  case "$1" in
    --trace) trace_file=$2; shift 2 ;;
    --png) png_file=$2; shift 2 ;;
    --pdf) pdf_file=$2; shift 2 ;;
    --summary) summary_file=$2; shift 2 ;;
    --run-id) run_id=$2; shift 2 ;;
    --invocation-id) invocation_id=$2; shift 2 ;;
    --module) module=$2; shift 2 ;;
    --file-tag) file_tag=$2; shift 2 ;;
    --sample-seconds) sample_seconds=$2; shift 2 ;;
    --plot-seconds) plot_seconds=$2; shift 2 ;;
    --parent-pid) parent_pid=$2; shift 2 ;;
    --ready-file) ready_file=$2; shift 2 ;;
    --stop-file) stop_file=$2; shift 2 ;;
    --status-file) status_file=$2; shift 2 ;;
    --done-file) done_file=$2; shift 2 ;;
    --error-file) error_file=$2; shift 2 ;;
    *) shift ;;
  esac
done

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
plotter="$script_dir/renderMemoryPlot.R"
summarizer="$script_dir/renderMemorySummary.awk"
rscript_bin=${INSPIIRED2_RSCRIPT:-}
[[ -n "$rscript_bin" ]] || rscript_bin=$(command -v Rscript 2>/dev/null || true)

fail_monitor(){
  local message=$1
  [[ -n "$error_file" ]] && printf '%s\n' "$message" > "$error_file"
  [[ -n "$done_file" ]] && : > "$done_file"
  exit 0
}

for required in trace_file png_file pdf_file summary_file run_id invocation_id module \
                parent_pid ready_file stop_file status_file done_file error_file; do
  [[ -n "${!required}" ]] || fail_monitor "missing required sampler argument: $required"
done
[[ -r "$summarizer" ]] || fail_monitor "memory summary helper was not found"
[[ "$parent_pid" =~ ^[0-9]+$ ]] || fail_monitor "invalid launcher process identifier"
command -v flock >/dev/null 2>&1 || fail_monitor "the flock command required for memory monitoring was not found"

# Only one sampler may own an output directory. A competing sampler exits
# cleanly and the launcher continues the biological module without monitoring.
monitor_dir=$(dirname "$ready_file")
lock_file="$monitor_dir/active.lock"
run_file="$monitor_dir/currentRunID.txt"
exec 9>"$lock_file" || fail_monitor "the memory-monitor lock file could not be opened"
flock -n 9 || fail_monitor "another memory monitor is already active for this output directory"


resolve_cgroup_dir(){
  local mount_root=$1 mount_point=$2 cgroup_path=$3 relative
  if [[ "$cgroup_path" == "$mount_root" ]]; then
    relative=""
  elif [[ "$mount_root" == "/" ]]; then
    relative=$cgroup_path
  elif [[ "$cgroup_path" == "$mount_root/"* ]]; then
    relative=${cgroup_path#"$mount_root"}
  else
    return 1
  fi
  printf '%s%s\n' "$mount_point" "$relative"
}


discover_cgroup(){
  local fake_dir=${INSPIIRED2_CGROUP_DIR:-}
  local proc_cgroup=${INSPIIRED2_TEST_PROC_CGROUP:-/proc/self/cgroup}
  local proc_mountinfo=${INSPIIRED2_TEST_PROC_MOUNTINFO:-/proc/self/mountinfo}
  local cgroup_path mount_info mount_root mount_point candidate container_hint=0

  if [[ -n "$fake_dir" ]]; then
    if [[ -r "$fake_dir/memory.current" ]]; then
      cgroup_version=2; cgroup_dir=$fake_dir; return 0
    elif [[ -r "$fake_dir/memory.usage_in_bytes" ]]; then
      cgroup_version=1; cgroup_dir=$fake_dir; return 0
    fi
    return 1
  fi

  [[ -r "$proc_cgroup" && -r "$proc_mountinfo" ]] || return 1
  [[ -e /.dockerenv || -e /run/.containerenv ]] && container_hint=1
  grep -Eq '(docker|containerd|kubepods|libpod|podman)' "$proc_cgroup" 2>/dev/null && container_hint=1

  cgroup_path=$(awk -F: '$1 == "0" && $2 == "" { print $3; exit }' "$proc_cgroup")
  mount_info=$(awk 'index($0, " - cgroup2 ") { print $4 "\t" $5; exit }' "$proc_mountinfo")
  if [[ -n "$cgroup_path" && -n "$mount_info" ]]; then
    IFS=$'\t' read -r mount_root mount_point <<< "$mount_info"
    candidate=$(resolve_cgroup_dir "$mount_root" "$mount_point" "$cgroup_path") || candidate=""
    [[ "$mount_root" != "/" ]] && container_hint=1
    if [[ -n "$candidate" && -r "$candidate/memory.current" && $container_hint -eq 1 ]]; then
      cgroup_version=2; cgroup_dir=$candidate; return 0
    fi
    # A private cgroup namespace may report 0::/ while mountinfo retains the
    # host-side mount root. In that layout the namespace root is the mountpoint.
    if [[ "$cgroup_path" == "/" && -r "$mount_point/memory.current" && $container_hint -eq 1 ]]; then
      cgroup_version=2; cgroup_dir=$mount_point; return 0
    fi
  fi

  cgroup_path=$(awk -F: '$2 ~ /(^|,)memory(,|$)/ { print $3; exit }' "$proc_cgroup")
  mount_info=$(awk 'index($0, " - cgroup ") && $0 ~ /(^|[ ,])memory([ ,]|$)/ { print $4 "\t" $5; exit }' "$proc_mountinfo")
  if [[ -n "$cgroup_path" && -n "$mount_info" ]]; then
    IFS=$'\t' read -r mount_root mount_point <<< "$mount_info"
    candidate=$(resolve_cgroup_dir "$mount_root" "$mount_point" "$cgroup_path") || candidate=""
    [[ "$mount_root" != "/" ]] && container_hint=1
    if [[ -n "$candidate" && -r "$candidate/memory.usage_in_bytes" && $container_hint -eq 1 ]]; then
      cgroup_version=1; cgroup_dir=$candidate; return 0
    fi
    if [[ "$cgroup_path" == "/" && -r "$mount_point/memory.usage_in_bytes" && $container_hint -eq 1 ]]; then
      cgroup_version=1; cgroup_dir=$mount_point; return 0
    fi
  fi

  return 1
}


discover_cgroup || fail_monitor "a container memory cgroup could not be identified"

# Select the shared run identifier only after acquiring the output-directory
# lock. Demultiplex begins a new pipeline run; later modules reuse that ID.
if [[ "$module" != "demultiplex" && -r "$run_file" ]]; then
  IFS= read -r recorded_run_id < "$run_file" || recorded_run_id=""
  [[ -n "$recorded_run_id" ]] && run_id=$recorded_run_id
fi
if [[ "$module" == "demultiplex" || ! -s "$run_file" ]]; then
  run_tmp="${run_file}.tmp_${invocation_id}"
  if ! printf '%s\n' "$run_id" > "$run_tmp" || ! mv -f "$run_tmp" "$run_file"; then
    rm -f "$run_tmp"
    fail_monitor "the memory-monitor run identifier could not be recorded"
  fi
fi

if [[ $cgroup_version -eq 2 ]]; then
  usage_file="$cgroup_dir/memory.current"
  limit_file="$cgroup_dir/memory.max"
  stat_file="$cgroup_dir/memory.stat"
  events_file="$cgroup_dir/memory.events"
  swap_file="$cgroup_dir/memory.swap.current"
else
  usage_file="$cgroup_dir/memory.usage_in_bytes"
  limit_file="$cgroup_dir/memory.limit_in_bytes"
  stat_file="$cgroup_dir/memory.stat"
  events_file="$cgroup_dir/memory.failcnt"
  swap_file="$cgroup_dir/memory.memsw.usage_in_bytes"
fi

header=$'run_id\tinvocation_id\ttimestamp_utc\tepoch_seconds\tinvocation_elapsed_seconds\tmodule\tfile_tag\tevent\texit_status\tcgroup_version\tmemory_bytes\tinactive_file_bytes\tworking_set_bytes\tanon_bytes\tfile_bytes\tshmem_bytes\tswap_bytes\tmemory_limit_bytes\tpercent_limit\toom_events\toom_kill_events'
if [[ -s "$trace_file" ]]; then
  IFS= read -r existing_header < "$trace_file"
  [[ "$existing_header" == "$header" ]] || fail_monitor "existing memory trace has an incompatible header"
else
  printf '%s\n' "$header" > "$trace_file" || fail_monitor "memory trace could not be created"
fi

start_epoch=$(date +%s)
last_report_epoch=0
termination_requested=0
finished=0


sample_memory(){
  local event=$1 exit_status=${2:-} now_epoch timestamp elapsed memory limit
  local anon=0 file_cache=0 shmem=0 inactive_file=0 swap=0 oom=0 oom_kill=0
  local rss=0 cache=0 total_rss=0 total_cache=0 total_shmem=0 total_inactive_file=0
  local key value working_set percent="" memsw=0

  memory=$(<"$usage_file") || return 1
  [[ "$memory" =~ ^[0-9]+$ ]] || return 1
  now_epoch=$(date +%s)
  timestamp=$(date -u +'%Y-%m-%dT%H:%M:%SZ')
  elapsed=$((now_epoch - start_epoch))

  if [[ -r "$stat_file" ]]; then
    while read -r key value _; do
      case "$key" in
        anon) anon=$value ;;
        file) file_cache=$value ;;
        shmem) shmem=$value ;;
        inactive_file) inactive_file=$value ;;
        rss) rss=$value ;;
        cache) cache=$value ;;
        total_rss) total_rss=$value ;;
        total_cache) total_cache=$value ;;
        total_shmem) total_shmem=$value ;;
        total_inactive_file) total_inactive_file=$value ;;
      esac
    done < "$stat_file"
  fi

  if [[ $cgroup_version -eq 1 ]]; then
    ((total_rss > 0)) && anon=$total_rss || anon=$rss
    ((total_cache > 0)) && file_cache=$total_cache || file_cache=$cache
    ((total_shmem > 0)) && shmem=$total_shmem
    ((total_inactive_file > 0)) && inactive_file=$total_inactive_file
  fi

  working_set=$((memory - inactive_file))
  ((working_set < 0)) && working_set=0

  limit=""
  if [[ -r "$limit_file" ]]; then
    limit=$(<"$limit_file")
    if [[ "$limit" == "max" || ! "$limit" =~ ^[0-9]+$ ]]; then
      limit=""
    elif [[ $cgroup_version -eq 1 && ${#limit} -ge 18 && $limit -ge 900000000000000000 ]]; then
      limit=""
    fi
  fi

  if [[ -r "$swap_file" ]]; then
    swap=$(<"$swap_file")
    [[ "$swap" =~ ^[0-9]+$ ]] || swap=0
    if [[ $cgroup_version -eq 1 ]]; then
      memsw=$swap; swap=$((memsw - memory)); ((swap < 0)) && swap=0
    fi
  fi

  if [[ $cgroup_version -eq 2 && -r "$events_file" ]]; then
    while read -r key value _; do
      [[ "$key" == "oom" ]] && oom=$value
      [[ "$key" == "oom_kill" ]] && oom_kill=$value
    done < "$events_file"
  elif [[ -r "$events_file" ]]; then
    oom=$(<"$events_file"); [[ "$oom" =~ ^[0-9]+$ ]] || oom=0
  fi

  [[ -n "$limit" ]] && percent=$(awk -v m="$memory" -v l="$limit" 'BEGIN { printf "%.4f", 100 * m / l }')

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$run_id" "$invocation_id" "$timestamp" "$now_epoch" "$elapsed" "$module" "$file_tag" \
    "$event" "$exit_status" "$cgroup_version" "$memory" "$inactive_file" "$working_set" \
    "$anon" "$file_cache" "$shmem" "$swap" "$limit" "$percent" "$oom" "$oom_kill" >> "$trace_file"
}


render_outputs(){
  local png_tmp pdf_tmp summary_tmp render_message
  local render_failures=""
  png_tmp="${png_file}.tmp_${invocation_id}"
  pdf_tmp="${pdf_file}.tmp_${invocation_id}"
  summary_tmp="${summary_file}.tmp_${invocation_id}"

  if awk -v runID="$run_id" -f "$summarizer" "$trace_file" > "$summary_tmp" &&
     [[ -s "$summary_tmp" ]] && mv -f "$summary_tmp" "$summary_file"; then
    :
  else
    rm -f "$summary_tmp"
    render_failures="summary"
  fi

  if [[ -z "$render_failures" && -n "$rscript_bin" && -r "$plotter" ]] &&
     "$rscript_bin" --vanilla "$plotter" "$trace_file" "$summary_file" \
       "$png_tmp" "$pdf_tmp" "$run_id" >/dev/null &&
     [[ -s "$png_tmp" && -s "$pdf_tmp" ]] &&
     mv -f "$pdf_tmp" "$pdf_file" && mv -f "$png_tmp" "$png_file"; then
    :
  else
    rm -f "$png_tmp" "$pdf_tmp"
    [[ -n "$render_failures" ]] && render_failures+=", "
    render_failures+="PNG/PDF report"
  fi

  if [[ -n "$render_failures" ]]; then
    render_message="memory report update failed for: $render_failures"
    if [[ ! -r "$error_file" ]] || ! grep -Fxq "$render_message" "$error_file"; then
      printf '%s\n' "$render_message" >> "$error_file"
    fi
  elif [[ -r "$error_file" ]] && [[ $(wc -l < "$error_file") -eq 1 ]] &&
       grep -q '^memory report update failed for:' "$error_file"; then
    rm -f "$error_file"
  fi
  last_report_epoch=$(date +%s)
}


finish_monitor(){
  ((finished == 1)) && return
  finished=1
  render_outputs
  : > "$done_file"
}

trap 'termination_requested=1' TERM INT
trap finish_monitor EXIT

sample_memory "start" "" || fail_monitor "container memory could not be sampled"
: > "$ready_file"
render_outputs

while :; do
  if [[ -e "$stop_file" || $termination_requested -eq 1 ]] || ! kill -0 "$parent_pid" 2>/dev/null; then
    exit_status=""
    [[ -r "$status_file" ]] && exit_status=$(<"$status_file")
    sample_memory "end" "$exit_status" || true
    break
  fi

  sleep "$sample_seconds" || true

  if [[ -e "$stop_file" || $termination_requested -eq 1 ]] || ! kill -0 "$parent_pid" 2>/dev/null; then
    exit_status=""
    [[ -r "$status_file" ]] && exit_status=$(<"$status_file")
    sample_memory "end" "$exit_status" || true
    break
  fi

  if ! sample_memory "sample" ""; then
    printf '%s\n' "container memory could no longer be sampled" > "$error_file"
    break
  fi
  now_epoch=$(date +%s)
  ((now_epoch - last_report_epoch >= plot_seconds)) && render_outputs
done

exit 0
