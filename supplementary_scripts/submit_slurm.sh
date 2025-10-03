#!/bin/bash

# The name of the job:
#SBATCH --job-name="freyr"

# The partition in which to submit the job:
#SBATCH --partition="batch"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# The total amount of memory in megabytes in the job:
#SBATCH --mem-per-cpu=4GB

# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL
# begins
#SBATCH --mail-type=BEGIN
# ends successfully
#SBATCH --mail-type=END

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=7-0:0:00

# Output errors and logs into same file
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out

# Fail loudly and propagate errors everywhere
set -Eeuo pipefail
shopt -s inherit_errexit || true   # bash ≥4.4; ok if it fails

# Log failing command + line number (also inside funcs/subshells)
trap 's=$?; echo "$(date "+%F %T") [ERR] exit=$s line=$LINENO cmd=${BASH_COMMAND@Q}" >&2' ERR
trap 's=$?; echo "$(date "+%F %T") [EXIT] status=$s" >&2' EXIT
set -o errtrace

# Optional: turn on xtrace when DEBUG=1
if [[ "${DEBUG:-0}" == "1" ]]; then
  export PS4='+ $(date "+%T") ${BASH_SOURCE##*/}:${LINENO}:${FUNCNAME[0]:-main}: '
  set -x
fi

# ---------- helpers ----------
log()       { echo "$(date '+%Y-%m-%d %H:%M:%S') [LOG] $*"; }
error_exit(){ echo "$(date '+%Y-%m-%d %H:%M:%S') [ERROR] $1"; exit "${2:-1}"; }

[[ -n "${SLURM_JOB_ID:-}" ]] || error_exit "Submit with sbatch" 1

# ---------- arg parsing ----------
is_allowed() {
  local param="$1"; shift
  local allowed_list=("$@")
  for allowed in "${allowed_list[@]}"; do
    [[ "$param" == "$allowed" ]] && return 0
  done
  return 1
}

allowed_args=(
  "primer_params"
  "read_dir"         # comma-separated list or can repeat --read_dir argument
  "primers"
)

declare -A vals
for a in "${allowed_args[@]}"; do vals["$a"]=""; done
nextflow_args=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --nextflow-*)
      param="${1#--nextflow-}"
      if [[ $# -gt 1 && ! "$2" =~ ^-- ]]; then nextflow_args+=("-$param" "$2"); shift
      else nextflow_args+=("-$param"); fi
      ;;
    --*)
      param="${1#--}"
      if is_allowed "$param" "${allowed_args[@]}"; then
        [[ $# -gt 1 && ! "$2" =~ ^-- ]] || error_exit "Missing value for --$param" 2
        # allow repeating --read_dir; append with comma
        if [[ "$param" == "read_dir" && -n "${vals[$param]}" ]]; then
          vals["$param"]+=",${2}"
        else
          vals["$param"]="$2"
        fi
        shift
      else
        error_exit "Invalid Freyr argument --$param" 2
      fi
      ;;
    *)
      error_exit "Unrecognized positional arg: $1" 2
      ;;
  esac
  shift
done

# ---------- locate read_dir(s) if not provided ----------
echo 'locating read directory'
if [[ -z "${vals[read_dir]}" ]]; then
  log "read_dir not provided; searching under ${SLURM_SUBMIT_DIR}/data for directories with R1 FASTQs..."
  mapfile -t auto_dirs < <(
    find "${SLURM_SUBMIT_DIR}/data" -type f -name '*.fastq.gz' \
      \( -name '*_R1_*' -o -name '*_R1.*' -o -name '*_1.fastq.gz' -o -name '*_1.fq.gz' \) \
      -printf '%h\n' | sort -u
  )
  [[ ${#auto_dirs[@]} -gt 0 ]] || error_exit "No directories with R1 FASTQs found under ${SLURM_SUBMIT_DIR}/data" 4
  IFS=',' vals[read_dir]="${auto_dirs[*]}"
  IFS=$' \t\n'
  log "Discovered read dirs: ${vals[read_dir]}"
fi

# ---------- required inputs ----------
[[ -n "${vals[primers]}" ]]        || error_exit "No primers provided (--primers)" 6
[[ -f "${vals[primer_params]}" ]]  || error_exit "Primer params file not found: ${vals[primer_params]}" 6

mkdir -p ./inputs
cp -f "${vals[primer_params]}" ./inputs/primer_params.csv

# ---------- sample discovery (from data dirs) ----------
output_csv=./inputs/Sample_info.csv
echo "sample,primers,fwd,rev" > "$output_csv"

# Function for mate detection 
find_mate() {
  local r1="$1" r2=""
  r2="${r1/_R1_/_R2_}"; [[ -f "$r2" ]] && { echo "$r2"; return; }
  r2="${r1/_R1./_R2.}"; [[ -f "$r2" ]] && { echo "$r2"; return; }
  r2="${r1/_1.fastq/_2.fastq}"; [[ -f "$r2" ]] && { echo "$r2"; return; }
  r2="${r1/_1.fq/_2.fq}";       [[ -f "$r2" ]] && { echo "$r2"; return; }
  r2="${r1/_1.fastq.gz/_2.fastq.gz}"; [[ -f "$r2" ]] && { echo "$r2"; return; }
  r2="${r1/_1.fq.gz/_2.fq.gz}";       [[ -f "$r2" ]] && { echo "$r2"; return; }
  echo ""
}

# Function for extracting sample names
extract_sample() {
  local core="$1"

  # Strip read indicator
  if [[ "$core" == *_R1_* ]]; then
    core="${core%%_R1_*}"
  elif [[ "$core" == *_R1.* ]]; then
    core="${core%%_R1.*}"
  elif [[ "$core" == *_1.fastq* || "$core" == *_1.fq* ]]; then
    core="${core%_1.fastq.gz}"; core="${core%_1.fq.gz}"
    core="${core%_1.fastq}";    core="${core%_1.fq}"
  fi

  # Strip typical lane/sample suffix tokens only (DO NOT drop any leading FCID)
  core="${core%%_L[0-9][0-9][0-9]*}"
  core="${core%%_S[0-9]*}"
  core="${core%%_001}"

  echo "$core"
}

IFS=',' read -r -a READDIRS <<< "${vals[read_dir]}"
IFS=$' \t\n'
maxdepth=6

log "Finding samples in read directories..."
for dir in "${READDIRS[@]}"; do
  [[ -d "$dir" ]] || { log "WARN: read_dir not found: $dir"; continue; }

  # enumerate R1s in this dir (NUL-safe)
  while IFS= read -r -d '' fwd; do
    rev="$(find_mate "$fwd")"
    if [[ -z "$rev" ]]; then
      log "WARN: missing R2 for $fwd"
      continue
    fi

    base="$(basename "$fwd")"
    sample="$(extract_sample "$base")"

    # write one row per pair
    printf "%s,%s,%s,%s\n" "$sample" "${vals[primers]}" "$fwd" "$rev" >> "$output_csv"
  done < <(
    find "$dir" -maxdepth "$maxdepth" -type f -name '*.fastq.gz' \
      \( -name '*_R1_*' -o -name '*_R1.*' -o -name '*_1.fastq.gz' -o -name '*_1.fq.gz' -o -name '*_1.fastq' -o -name '*_1.fq' \) \
      ! -iname '*undetermined*' \
      -print0 | sort -z -u
  )
done

rows=$(( $(wc -l < "$output_csv") - 1 ))
[[ $rows -gt 0 ]] || error_exit "No valid FASTQ pairs discovered in the provided read_dir(s)" 7
log "Discovered $rows read pairs. Wrote $output_csv"

# ---------- run freyr ----------
module load Java/17.0.6
module load shifter/22.02.1

# Print nextflow args
if ((${#nextflow_args[@]})); then
  log "Running Nextflow with args:"
  for a in "${nextflow_args[@]}"; do printf '  %q\n' "$a"; done
fi

# Run nextflow
NXF_VER=23.05.0-edge \
  nextflow run . \
    --samplesheet "$output_csv" \
    --primer_params ./inputs/primer_params.csv \
    --slurm_account "${SLURM_JOB_ACCOUNT}" \
    -profile basc_slurm,error1_backoff \
    "${nextflow_args[@]}"

# ---------- job stats ----------
/usr/local/bin/showJobStats.scr