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

set -e

# Helper Functions
log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') [LOG] $1"
}

error_exit() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') [ERROR] $1"
    exit "$2"
}

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi


###### PARSE AND VALIDATE ARGUMENTS
# Function to check if a parameter is allowed
is_allowed() {
    local param="$1"
    shift
    local allowed_list=("$@")
    for allowed in "${allowed_list[@]}"; do
        if [[ "$param" == "$allowed" ]]; then
            return 0
        fi
    done
    return 1
}

# Define allowed arguments
allowed_samplesheet_args=( \
    "samplesheet" \
    "primer_params" \
    "read_dir" \
    "primers"
)

# Initialize samplesheet values (associative array)
declare -A samplesheet_values

# Initialize all allowed tool arguments with empty values
for arg in "${allowed_samplesheet_args[@]}"; do
    samplesheet_values["$arg"]=""
done



while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --nextflow-*)
            param="${1#--nextflow-}"  # Strip --nextflow- prefix
            if [[ "$#" -gt 1 && ! "$2" == --* ]]; then
                nextflow_args+=("--$param" "$2")
                shift
            else
                nextflow_args+=("--$param")
            fi
            ;;
        *)
            param="${1#--}"
            if is_allowed "$param" "${allowed_samplesheet_args[@]}"; then
                if [[ "$#" -gt 1 && ! "$2" == --* ]]; then
                    samplesheet_values["$param"]="$2"
                    shift
                else
                    echo "Error: Missing value for $param"
                    exit 1
                fi
            else
                echo "Error: Invalid Freyr argument --$param"
                exit 1
            fi
            ;;
    esac
    shift
done

###### LOCATE SAMPLESHEETS AND READ DIR IF NOT PROVIDED

# Find samplesheet if not provided
if [[ -z "${samplesheet_values["samplesheet"]}" ]]; then
  log "samplesheet not provided; searching for files..."
  sample_sheets=$(find "${SLURM_SUBMIT_DIR}/data" -maxdepth 2 -name '*SampleSheet*.csv' -type f)
  if [[ -n "$sample_sheets" ]]; then
      samplesheet_values["samplesheet"]="$(echo $sample_sheets | tr ' ' ',')"
  else
      error_exit "No suitable sample sheets found" 4
  fi
fi

# Find read_dir if not provided
if [[ -z "${samplesheet_values["read_dir"]}" ]]; then
  log "read_dir not provided; searching for directories..."
  read_dir=""
  for sheet in ${samplesheet_values["samplesheet"]//,/ }; do
      dir=$(dirname "$sheet")
      read_dir="${read_dir}${dir} "
  done
  
  if [[ -n "$read_dir" ]]; then
      samplesheet_values["read_dir"]="$(echo $read_dir | tr ' ' ',')"
  else
      error_exit "No suitable read_dir found" 5
  fi
fi

# Break if primers have not been provided
if [[ -z "${samplesheet_values["primers"]}" ]]; then
  error_exit "No primers provided (--primers)" 6
fi

# Copy primer params to inputs and break if not provided
if [[ -f "${samplesheet_values["primer_params"]}" ]]; then
    cp "${samplesheet_values["primer_params"]}" ./inputs/primer_params.csv
else
    error_exit "Primer params file not found: ${samplesheet_values["primer_params"]}" 6
fi


###### PARSE SAMPLESHEETS
echo "Creating freyr samplesheet from miseq samplesheet"
parse_miseq_samplesheet() {
  local SampleSheet="$1"
  local primers="$2"
  local read_dir="$3"

  if [ -z "$SampleSheet" ]; then
    echo "Error: need to provide a MiSeq SampleSheet.csv file" >&2
    return 1
  fi

  # Locate [Data] section
  data_line=$(grep -n "^\[Data\]" "$SampleSheet" | cut -d: -f1)
  if [ -z "$data_line" ]; then
    echo "Error: [Data] section not found in $SampleSheet" >&2
    return 1
  fi

  header_line=$((data_line + 1))
  data_start=$((header_line + 1))

  # Iterate through samples, remove fastqs if they dont have matching files
  tail -n +"$data_start" "$SampleSheet" | cut -d',' -f1 | grep -v '^$' | while read -r sample; do
    fwd_pattern="${read_dir}/*${sample}*_R1*.fastq.gz"
    rev_pattern="${read_dir}/*${sample}*_R2*.fastq.gz"

    fwd_file=$(ls $fwd_pattern 2>/dev/null | head -n 1)
    rev_file=$(ls $rev_pattern 2>/dev/null | head -n 1)

    if [[ -f "$fwd_file" && -f "$rev_file" ]]; then
      # Extract actual FASTQ prefix from filename (before _R1)
      fq_basename=$(basename "$fwd_file")
      actual_sample=$(echo "$fq_basename" | sed -E 's/_R1.*//')
      echo "${actual_sample},${primers},${read_dir}"
    else
      >&2 echo "$(date '+%Y-%m-%d %H:%M:%S') [WARN] Skipping sample '${sample}' – FASTQ file(s) missing"
    fi
  done
}

# Convert values values to ordered list
samplesheet_values_list=()
for arg in "${allowed_samplesheet_args[@]}"; do
    samplesheet_values_list+=("${samplesheet_values[$arg]}")
done

mkdir -p ./inputs
output_csv=./inputs/Sample_info.csv
echo "sample,primers,read_dir" > "$output_csv"

# Apply to all detected SampleSheets
IFS=',' read -ra SHEETS <<< "${samplesheet_values["samplesheet"]}"
IFS=',' read -ra READDIRS <<< "${samplesheet_values["read_dir"]}"
for i in "${!SHEETS[@]}"; do
  sheet="${SHEETS[$i]}"
  dir="${READDIRS[$i]}"
  primers="${samplesheet_values["primers"]}"

  parse_miseq_samplesheet "$sheet" "$primers" "$dir" >> "$output_csv"
done

if [[ $(wc -l < "$output_csv") -le 1 ]]; then
    error_exit "No valid samples found after FASTQ validation. Exiting." 7
fi

# Load modules
module load Java/17.0.6
module load shifter/22.02.1

log "Running Nextflow with arguments: ${nextflow_args[*]}"
NXF_VER=23.05.0-edge \
    nextflow run . \
    --samplesheet ./inputs/Sample_info.csv \
    --primer_params ./inputs/primer_params.csv \
    --slurm_account "${SLURM_JOB_ACCOUNT}" \
    -profile basc_slurm,error1_backoff \
    "${nextflow_args[@]}"

# Output useful job stats
/usr/local/bin/showJobStats.scr 