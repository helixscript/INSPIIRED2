BEGIN { FS = OFS = "\t" }

NR == 1 {
  for(i = 1; i <= NF; i++) col[$i] = i
  next
}

$(col["run_id"]) == runID {
  id = $(col["invocation_id"])
  if(!(id in seen)){
    seen[id] = ++n
    order[n] = id
    module[id] = $(col["module"])
    fileTag[id] = $(col["file_tag"])
    startStamp[id] = $(col["timestamp_utc"])
    startEpoch[id] = $(col["epoch_seconds"]) + 0
  }

  endStamp[id] = $(col["timestamp_utc"])
  endEpoch[id] = $(col["epoch_seconds"]) + 0
  samples[id]++

  memory = $(col["memory_bytes"]) + 0
  working = $(col["working_set_bytes"]) + 0
  anon = $(col["anon_bytes"]) + 0
  shmem = $(col["shmem_bytes"]) + 0
  limit = $(col["memory_limit_bytes"])

  if(memory > peakMemory[id]) peakMemory[id] = memory
  if(working > peakWorking[id]) peakWorking[id] = working
  if(anon > peakAnon[id]) peakAnon[id] = anon
  if(shmem > peakShmem[id]) peakShmem[id] = shmem
  if(limit ~ /^[0-9]+$/ && limit + 0 > 0) memoryLimit[id] = limit
  if($(col["event"]) == "end") exitStatus[id] = $(col["exit_status"])
}

END {
  print "run_id", "invocation_id", "module", "file_tag", "start_utc", "end_utc", \
        "duration_seconds", "samples", "peak_memory_bytes", "peak_working_set_bytes", \
        "peak_anon_bytes", "peak_shmem_bytes", "memory_limit_bytes", "peak_percent_limit", \
        "exit_status"

  for(i = 1; i <= n; i++){
    id = order[i]
    percent = ""
    if(memoryLimit[id] + 0 > 0) percent = sprintf("%.4f", 100 * peakMemory[id] / memoryLimit[id])
    print runID, id, module[id], fileTag[id], startStamp[id], endStamp[id], \
          endEpoch[id] - startEpoch[id], samples[id], sprintf("%.0f", peakMemory[id]), \
          sprintf("%.0f", peakWorking[id]), sprintf("%.0f", peakAnon[id]), \
          sprintf("%.0f", peakShmem[id]), memoryLimit[id], percent, exitStatus[id]
  }
}
