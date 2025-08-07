#!/bin/bash
#SBATCH --job-name=step2_hepc
#SBATCH --account=proj-gm0001
#SBATCH --array=0-100%9999  #now there are 36 jobs + PacBio
#SBATCH --output=/stornext/snfs170/next-gen/scratch/adrian/projects/r2t/results_hepc_test2/hepc_test_%A_%a.out
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=adrianalexandro777@gmail.com
#SBATCH --export=ALL


# ---------- ENTORNO ----------
source ~/.bashrc
conda activate my_env


# ---------- RUTAS ----------
PROJECT=/stornext/snfs170/next-gen/scratch/adrian/projects/r2t
READS_DIR=${PROJECT}/results_reads_hepc
OUT_ROOT=${PROJECT}/results_hepc_test2          
LOGDIR=${OUT_ROOT}
PREFIX=HepC
mkdir -p "${LOGDIR}"

TASK=${SLURM_ARRAY_TASK_ID:-0}
TMPDIR=/space1/tmp/r2t_hepc_${SLURM_JOB_ID}_${TASK}
METRICFILE=${LOGDIR}/v2t_step2_%A_%a.metrics
TIMECMD=/stornext/snfs170/next-gen/scratch/adrian/local/bin/time
METRICFILE_RESOLVED="${METRICFILE//%A/${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}}"
METRICFILE_RESOLVED="${METRICFILE_RESOLVED//%a/${SLURM_ARRAY_TASK_ID:-$TASK}}"

RECOVERY_TSV="${OUT_ROOT}/missing_labels.tsv"
mapfile -t PACBIO_SAMPLES < <(grep -vE '^\s*(#|$)' pacbio_samples)
in_list() { local needle="$1"; shift; for x in "$@"; do [[ "$x" == "$needle" ]] && return 0; done; return 1; }


# Find R1/R2 in READS_DIR from LABEL if metrics didn't record them or files moved.
find_r1_r2_from_label () {
  local label="$1"
  local out_r1 out_r2
  out_r1=""; out_r2=""

  # Paired-end patterns
  for r1c in \
    "${READS_DIR}/${label}_1.fastq.gz" "${READS_DIR}/${label}_1.fastq" \
    "${READS_DIR}/${label}_R1.fastq.gz" "${READS_DIR}/${label}_R1.fastq"; do
    [[ -f "$r1c" ]] || continue

    # base without _1/_R1
    local base="$r1c"
    [[ "$base" == *_R1.fastq.gz ]] && base="${base%_R1.fastq.gz}" || true
    [[ "$base" == *_R1.fastq    ]] && base="${base%_R1.fastq}"    || true
    [[ "$base" == *_1.fastq.gz  ]] && base="${base%_1.fastq.gz}"  || true
    [[ "$base" == *_1.fastq     ]] && base="${base%_1.fastq}"     || true

    for r2c in "${base}_2.fastq.gz" "${base}_2.fastq" "${base}_R2.fastq.gz" "${base}_R2.fastq"; do
      if [[ -f "$r2c" ]]; then
        out_r1="$r1c"; out_r2="$r2c"
        echo "$out_r1|$out_r2"
        return 0
      fi
    done
  done

  # Single-end patterns
  for sec in "${READS_DIR}/${label}.fastq.gz" "${READS_DIR}/${label}.fastq"; do
    if [[ -f "$sec" ]]; then
      echo "$sec|"
      return 0
    fi
  done

  echo "|"
  return 1
}


# ---------- BUILD JOB LIST ----------

shopt -s nullglob extglob

declare -a LABELS TYPES R1S R2S
declare -A SEEN


# Expect header "FILE  LABEL  SAMPLE  JOB_ID"
mapfile -t METRICS_MISSING < <(awk -v root="$OUT_ROOT" 'NR>1 {
  file=$1
  if (file ~ /^\//) {
    print file          # absolute path stays
  } else {
    print root "/" file # make relative names absolute
  }
}' "$RECOVERY_TSV")



for mf in "${METRICS_MISSING[@]}"; do
  [[ -f "$mf" ]] || { echo "[WARN] Missing metrics file listed in TSV: $mf" >&2; continue; }

  # Extract LABEL / READ_TYPE / R1 / R2 from the metrics header
  label=$(awk -F': ' '/^LABEL[[:space:]]*:/     {print $2; exit}' "$mf")
  read_type=$(awk -F': ' '/^READ_TYPE[[:space:]]*:/ {print $2; exit}' "$mf")
  r1=$(awk -F': ' '/^R1[[:space:]]*:/        {print $2; exit}' "$mf" | head -n1)
  r2=$(awk -F': ' '/^R2[[:space:]]*:/        {print $2; exit}' "$mf" | head -n1)

  # Sanity: skip if no LABEL
  [[ -n "${label:-}" ]] || { echo "[WARN] No LABEL in $mf; skipping" >&2; continue; }
  [[ -n "${SEEN[$label]:-}" ]] && continue

  # If files from metrics header no longer exist, try to rediscover in READS_DIR
  if [[ -n "${r1:-}" ]] && [[ ! -f "$r1" ]]; then
    echo "[INFO] R1 not found from metrics; rediscovering for $label" >&2
    r1r2=$(find_r1_r2_from_label "$label")
    r1="${r1r2%%|*}"; r2="${r1r2#*|}"
  fi
  if [[ -z "${r1:-}" ]]; then
    echo "[ERROR] Could not resolve R1 for $label; skipping" >&2
    continue
  fi

  # Infer type if missing/incorrect based on presence of R2
  if [[ -z "${read_type:-}" ]]; then
    read_type=$([[ -n "${r2:-}" ]] && echo "pe_short" || echo "se_short")
  fi

  sample="${label#${PREFIX}_}"
  if in_list "$sample" "${PACBIO_SAMPLES[@]}"; then
    read_type="pacbio"
    r2=""
    echo "[INFO] OVERRIDE (PACBIO set): $label -> READ_TYPE=pacbio" >&2
  fi

  SEEN[$label]=1

  LABELS+=("$label")
  TYPES+=("$read_type")
  R1S+=("$r1")
  R2S+=("${r2:-}")
done

## Add pacbio that are not in missings
for s in "${PACBIO_SAMPLES[@]}"; do
  label="${PREFIX}_${s}"
  [[ -n "${SEEN[$label]:-}" ]] && continue  # ya agregado vía missing + override

  r1r2=$(find_r1_r2_from_label "$label")
  r1="${r1r2%%|*}"; r2=""  # ignorar R2 para pacbio
  if [[ -z "$r1" ]]; then
    echo "[WARN] PACBIO sin FASTQ encontrado: $label; skip" >&2
    continue
  fi
  SEEN[$label]=1
  LABELS+=("$label")
  TYPES+=("pacbio")  # pacbio read type
  R1S+=("$r1")
  R2S+=("")
  echo "[INFO] ADD (PACBIO-only): $label" >&2
done


TOTAL=${#LABELS[@]}



# ---------- SAVES ----------
if (( TOTAL == 0 )); then
  echo "[ERROR] No jobs were built (missing + pacbio). Check TSV/paths." >&2
  exit 1
fi

if (( TASK >= TOTAL )); then
  echo "[INFO] TASK ${TASK} out of range (TOTAL=${TOTAL}). Nothing to do." >&2
  exit 0
fi

LABEL="${LABELS[$TASK]}"
READ_TYPE="${TYPES[$TASK]}"
R1="${R1S[$TASK]}"
R2="${R2S[$TASK]}"
SAMPLE="${LABEL#${PREFIX}_}"

# HERE INFORMATIVE LOG
{
  echo "=== $(date '+%F %T') ==="
  echo "JOB_ID        : ${SLURM_JOB_ID}_${TASK}"
  echo "LABEL         : ${LABEL}"
  echo "SAMPLE        : ${SAMPLE}"
  echo "READ_TYPE     : ${READ_TYPE}"
  echo "R1            : ${R1}"
  [[ -n "$R2" ]] && echo "R2            : ${R2}" || true
  echo "THREADS       : ${SLURM_CPUS_PER_TASK}"
  echo "TMPDIR        : ${TMPDIR}"
} >> "${METRICFILE_RESOLVED}"

# ---------- The execution ----------
if [[ "$READ_TYPE" == "pe_short" ]]; then
  "${TIMECMD}" -v \
  srun -n 1 -c "${SLURM_CPUS_PER_TASK}" --cpu-bind=cores --mem-bind=local \
  virus2tree_step2.sh \
      --read_type pe_short \
      --reads "$R1" "$R2" \
      --threads "${SLURM_CPUS_PER_TASK}" \
      --root_dir "${OUT_ROOT}" \
      --temp_dir "${TMPDIR}" \
      --out_dir r2t_ref \
      2>>"${METRICFILE_RESOLVED}"
else
  "${TIMECMD}" -v \
  srun -n 1 -c "${SLURM_CPUS_PER_TASK}" --cpu-bind=cores --mem-bind=local \
  virus2tree_step2.sh \
      --read_type "${READ_TYPE}" \
      --reads "$R1" \
      --threads "${SLURM_CPUS_PER_TASK}" \
      --root_dir "${OUT_ROOT}" \
      --temp_dir "${TMPDIR}" \
      --out_dir r2t_ref \
      2>>"${METRICFILE_RESOLVED}"
fi


sleep 30
sacct -j "${SLURM_JOB_ID}" \
      --format=JobIDRaw,JobName,State,ExitCode,Elapsed,TotalCPU,CPUTime,MaxRSS,MaxVMSize,AveRSS,AveCPU \
      >> "${METRICFILE_RESOLVED}"

