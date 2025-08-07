#!/bin/bash
#SBATCH --job-name=step2_hepc
#SBATCH --account=proj-gm0001
#SBATCH --output=/stornext/snfs170/next-gen/scratch/adrian/projects/r2t/results_hepc_test2/snd_recov_%A.out
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
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
# Read sample IDs from the exterior, to ensure diferent job ids for each sample
: "${SAMPLE:?Set SAMPLE via sbatch --export=ALL,SAMPLE=<id>}"
PREFIX=HepC
mkdir -p "${LOGDIR}"

TASK=${SLURM_ARRAY_TASK_ID:-0}
TMPDIR=/space1/tmp/r2t_hepc_${SLURM_JOB_ID}_${TASK}
METRICFILE=${LOGDIR}/snd_recov_%A.metrics
TIMECMD=/stornext/snfs170/next-gen/scratch/adrian/local/bin/time
METRICFILE_RESOLVED="${METRICFILE//%A/${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}}"
METRICFILE_RESOLVED="${METRICFILE_RESOLVED//%a/${SLURM_ARRAY_TASK_ID:-$TASK}}"

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


PACBIO=/stornext/snfs170/next-gen/scratch/adrian/projects/r2t/results_hepc_test2/pacbio_samples

# Read first column of TSV, skipping comments/empty lines
# Read PacBio file
mapfile -t PACBIO_SAMPLES < <(grep -vE '^\s*(#|$)' $PACBIO)
label="${PREFIX}_${SAMPLE}"
# find read sets
r1r2=$(find_r1_r2_from_label "$label")
R1="${r1r2%%|*}"; R2="${r1r2#*|}"
[[ -n "$R1" ]] || { echo "[ERR] No R1 found for label=$label in $READS_DIR"; exit 2; }

if in_list "$SAMPLE" "${PACBIO_SAMPLES[@]}"; then
  READ_TYPE="pacbio"
  R2=""
else
  READ_TYPE=$([[ -n "$R2" ]] && echo "pe_short" || echo "se_short")
fi
LABEL="$label"




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
set +e
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
rc=$?
set -e

echo "ExitCode_cmd  : ${rc}" >> "${METRICFILE_RESOLVED}"

# ---------- MÉTRICAS SLURM (post-job) ----------
sleep 30
sacct -j "${SLURM_JOB_ID}" \
      --format=JobIDRaw,JobName,State,ExitCode,Elapsed,TotalCPU,CPUTime,MaxRSS,MaxVMSize,AveRSS,AveCPU \
      >> "${METRICFILE_RESOLVED}"

exit "$rc"
