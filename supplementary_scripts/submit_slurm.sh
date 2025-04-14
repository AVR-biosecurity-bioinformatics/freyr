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

# Define allowed arguments
allowed_tool_args=( \
    "sample_sheet" \
    "run_parameters" \
    "read_dir" \
    "pcr_primers" \
    "for_primer_seq" \
    "rev_primer_seq" \
    "target_gene" \
    "max_primer_mismatch" \
    "read_min_length" \
    "read_max_length" \
    "read_max_ee" \
    "read_trunc_length" \
    "read_trim_left" \
    "read_trim_right" \
    "asv_min_length" \
    "asv_max_length" \
    "high_sensitivity" \
    "concat_unmerged" \
    "genetic_code" \
    "coding" \
    "phmm" \
    "idtaxa_db" \
    "ref_fasta" \
    "idtaxa_confidence" \
    "run_blast" \
    "blast_min_identity" \
    "blast_min_coverage" \
    "target_kingdom" \
    "target_phylum" \
    "target_class" \
    "target_order" \
    "target_family" \
    "target_genus" \
    "target_species" \
    "min_sample_reads" \
    "min_taxa_reads" \
    "min_taxa_ra"
)

allowed_nextflow_args=( \
    "profile" \
    "resume" \
    "with-trace"
)

# Initialize tool values (associative array)
declare -A tool_values
nextflow_args=()

# Initialize all allowed tool arguments with empty values
for arg in "${allowed_tool_args[@]}"; do
    tool_values["$arg"]=""
done

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

# Parse and validate arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --nextflow-*)
            # Parse Nextflow-specific arguments
            param="${1#--nextflow-}"  # Strip the prefix
            if is_allowed "$param" "${allowed_nextflow_args[@]}"; then
                if [[ "$#" -gt 1 && ! "$2" == --* ]]; then
                    # Add parameter and value
                    nextflow_args+=("--$param" "$2")
                    shift  # Skip the value
                else
                    # Add parameter only
                    nextflow_args+=("--$param")
                fi
            else
                echo "Error: Invalid Nextflow argument --nextflow-$param"
                exit 1
            fi
            ;;
        *)
            # Parse arguments for the tool
            param="${1#--}"  # Strip any leading "--" (optional)
            if is_allowed "$param" "${allowed_tool_args[@]}"; then
                if [[ "$#" -gt 1 && ! "$2" == --* ]]; then
                    # Assign the value to the tool parameter
                    tool_values["$param"]="$2"
                    shift  # Skip the value
                else
                    echo "Error: Missing value for $param"
                    exit 1
                fi
            else
                echo "Error: Invalid sample sheet creation argument $param"
                exit 1
            fi
            ;;
    esac
    shift
done

# If samplesheet, run_parameters, and 
if [[ -z "${tool_values["sample_sheet"]}" ]]; then
  
  log "sample_sheet not provided; searching for files..."
  sample_sheets=$(find "${SLURM_SUBMIT_DIR}/data" -maxdepth 2 -name '*SampleSheet*.csv' -type f)
  if [[ -n "$sample_sheets" ]]; then
      tool_values["sample_sheet"]="$(echo $sample_sheets | tr ' ' ',')"
  else
      error_exit "No suitable sample sheets found" 4
  fi
fi

if [[ -z "${tool_values["run_parameters"]}" ]]; then
  log "run_parameters not provided; searching for files..."
  run_parameters=""
  for sheet in ${tool_values["sample_sheet"]//,/ }; do
      dir=$(dirname "$sheet")
      params=$(find "$dir" -maxdepth 2 -name '*unParameters.xml' -type f)
      run_parameters="${run_parameters}${params} "
  done

  if [[ -n "$run_parameters" ]]; then
      tool_values["run_parameters"]="$(echo $run_parameters | tr ' ' ',')"
  else
      error_exit "No suitable run parameters found" 5
  fi
fi


if [[ -z "${tool_values["read_dir"]}" ]]; then
  log "read_dir not provided; searching for files..."
  data_dirs=""
  for sheet in ${tool_values["sample_sheet"]//,/ }; do
      dir=$(dirname "$sheet")
      read_dir="${read_dir}${dir} "
  done
  
  if [[ -n "read_dir" ]]; then
      tool_values["read_dir"]="$(echo $read_dir | tr ' ' ',')"
  else
      error_exit "No suitable read_dir found" 5
  fi
fi

# Convert tool values to ordered list
tool_values_list=()
for arg in "${allowed_tool_args[@]}"; do
    tool_values_list+=("${tool_values[$arg]}")
done

# Load modules
#module load R/4.4.2-gfbf-2024a
module load Java/17.0.6
module load shifter/22.02.1

# Run the tool
log "Running sample sheet creation with arguments: ${tool_values_list[@]}"
shifter --image=jackscanlan/piperline-multi:0.0.1 -- Rscript supplementary_scripts/create_inputs.R "${tool_values_list[@]}"

# Run Nextflow
log "Running Nextflow with arguments: ${nextflow_args[@]}"
# make sure you replace the square bracketed file names (including the brackets) with the names of the files you made earlier
NXF_VER=23.04.5 \
    nextflow run . \
    --samplesheet ./inputs/Sample_info.csv \
    --loci_params ./inputs/loci_params.csv \
    --slurm_account ${SLURM_JOB_ACCOUNT} \
    -profile basc_slurm

# once the dataset has run, clean up your analysis directory
#rm -rf ${SLURM_SUBMIT_DIR}/output/modules/* ${SLURM_SUBMIT_DIR}/output/*.html ${SLURM_SUBMIT_DIR}/work/*

# Output useful job stats
/usr/local/bin/showJobStats.scr 