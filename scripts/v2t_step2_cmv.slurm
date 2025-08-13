#!/bin/bash
#SBATCH --job-name=step2_cmv
#SBATCH --account=proj-gm0001
#SBATCH --array=0-699%9999  #for 684 jobs
#SBATCH --output=/stornext/snfs170/next-gen/scratch/adrian/projects/r2t/results_cmv_step2/cmv_step2_%A_%a.out
#SBATCH --partition=long
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=adrianalexandro777@gmail.com
#SBATCH --export=ALL


source ~/.bashrc
conda activate my_env

PROJECT=/stornext/snfs170/next-gen/scratch/adrian/projects/r2t
READS_DIR=${PROJECT}/results_reads_cmv
OUT_ROOT=${PROJECT}/results_cmv_step2          # usar ruta absoluta (evita problemas de cwd)
LOGDIR=${OUT_ROOT}
PREFIX=CMV
mkdir -p "${LOGDIR}"

TASK=${SLURM_ARRAY_TASK_ID:-0}
TMPDIR=/space1/tmp/r2t_cmv_${SLURM_JOB_ID}_${TASK}
METRICFILE=${LOGDIR}/cmv_step2_%A_%a.metrics
TIMECMD=/stornext/snfs170/next-gen/scratch/adrian/local/bin/time
METRICFILE_RESOLVED="${METRICFILE//%A/${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}}"
METRICFILE_RESOLVED="${METRICFILE_RESOLVED//%a/${SLURM_ARRAY_TASK_ID:-$TASK}}"


shopt -s nullglob extglob

declare -a LABELS TYPES R1S R2S
declare -A SEEN

# 1) Paired-end: *_1.fastq{,.gz} con su *_2
for r1 in ${READS_DIR}/${PREFIX}_*_1.fastq ${READS_DIR}/${PREFIX}_*_1.fastq.gz \
          ${READS_DIR}/${PREFIX}_*_R1.fastq ${READS_DIR}/${PREFIX}_*_R1.fastq.gz; do
  [[ -f "$r1" ]] || continue

  # Normaliza base sin sufijo de R1
  if [[ "$r1" == *"_R1.fastq.gz" ]]; then base="${r1%_R1.fastq.gz}"
  elif [[ "$r1" == *"_R1.fastq"    ]]; then base="${r1%_R1.fastq}"
  elif [[ "$r1" == *"_1.fastq.gz"  ]]; then base="${r1%_1.fastq.gz}"
  else                                      base="${r1%_1.fastq}"
  fi

  # Busca R2 en cualquiera de las extensiones/estilos
  r2=""
  for cand in "${base}_2.fastq.gz" "${base}_2.fastq" "${base}_R2.fastq.gz" "${base}_R2.fastq"; do
    [[ -f "$cand" ]] && r2="$cand" && break
  done
  if [[ -z "$r2" ]]; then
    echo "[WARN] Falta el par R2 para: $r1; se omite." >&2
    continue
  fi

  label="$(basename "${base}")"   # p.ej., cmv_ERR118961
  [[ -n "${SEEN[$label]:-}" ]] && continue
  SEEN[$label]=1

  LABELS+=("$label")
  TYPES+=("-ax sr")
  R1S+=("$r1")
  R2S+=("$r2")
done

# 2) Single-end: excluye explícitamente *_1/_2 y *_R1/_R2 con extglob
for se in ${READS_DIR}/${PREFIX}_!(*_[12]|*_R[12]).fastq \
          ${READS_DIR}/${PREFIX}_!(*_[12]|*_R[12]).fastq.gz; do
  [[ -f "$se" ]] || continue
  label="$(basename "${se}")"
  label="${label%.fastq.gz}"
  label="${label%.fastq}"

  [[ -n "${SEEN[$label]:-}" ]] && continue
  SEEN[$label]=1

  LABELS+=("$label")
  TYPES+=("-ax sr")
  R1S+=("$se")
  R2S+=("")
done

TOTAL=${#LABELS[@]}

# ---------- GUARDAS ----------
if (( TOTAL == 0 )); then
  echo "[ERROR] No se encontraron FASTQ(s) en ${READS_DIR} con prefijo ${PREFIX}_." >&2
  exit 1
fi

if (( TASK >= TOTAL )); then
  echo "[INFO] TASK ${TASK} fuera de rango (TOTAL=${TOTAL}). Nada que hacer." >&2
  exit 0
fi

# ---------- SELECCIÓN DEL TRABAJO ----------
LABEL="${LABELS[$TASK]}"
READ_TYPE="${TYPES[$TASK]}"
R1="${R1S[$TASK]}"
R2="${R2S[$TASK]}"

SAMPLE="${LABEL#${PREFIX}_}"

PACBIO_SAMPLES=(ERR5052619 SRR13587450 SRR13587451 SRR13587452)
ONT_SAMPLES=(SRR23882357 SRR23882358 ERR2357678 SRR15047131 SRR15047130 SRR15047132 SRR15047133)
in_list() { local needle="$1"; shift; for x in "$@"; do [[ "$x" == "$needle" ]] && return 0; done; return 1; }

### [LONG-READ OVERRIDE] 
if in_list "$SAMPLE" "${PACBIO_SAMPLES[@]}"; then
  READ_TYPE="-ax map-pb"
  R2=""
  echo "[INFO] OVERRIDE: $SAMPLE -> READ_TYPE=-ax map-pb (SE)" >> "${METRICFILE_RESOLVED}"
elif in_list "$SAMPLE" "${ONT_SAMPLES[@]}"; then
  READ_TYPE="-ax map-ont"
  R2=""
  echo "[INFO] OVERRIDE: $SAMPLE -> READ_TYPE=-ax map-ont (SE)" >> "${METRICFILE_RESOLVED}"
fi

# ---------- LOG INICIAL ----------
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

# ---------- EJECUCIÓN + PERFILADO ----------
set +e
if [[ "$READ_TYPE" == "-ax sr" ]]; then
  "${TIMECMD}" -v \
  v2t-step2 \
      --read_type paired \
      -map_op "-ax sr" \
      --reads "$R1" "$R2" \
      --threads "${SLURM_CPUS_PER_TASK}" \
      --root_dir "${OUT_ROOT}" \
      --temp_dir "${TMPDIR}" \
      --out_dir r2t_ref \
      2>>"${METRICFILE_RESOLVED}"
else
  "${TIMECMD}" -v \
  v2t-step2 \
      --read_type single \
      -map_op "$READ_TYPE" \
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
