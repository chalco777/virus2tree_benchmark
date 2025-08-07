#!/bin/bash
#SBATCH --job-name=step2_hepc
#SBATCH --account=proj-gm0001
#SBATCH --array=0-799%9999  #fill all possible nodes
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
OUT_ROOT=${PROJECT}/results_hepc_test2           # usar ruta absoluta (evita problemas de cwd)
LOGDIR=${OUT_ROOT}
PREFIX=HepC
mkdir -p "${LOGDIR}"

TASK=${SLURM_ARRAY_TASK_ID:-0}
TMPDIR=/space1/tmp/r2t_hepc_${SLURM_JOB_ID}_${TASK}
METRICFILE=${LOGDIR}/v2t_step2_%A_%a.metrics
TIMECMD=/stornext/snfs170/next-gen/scratch/adrian/local/bin/time


# ---------- INVENTARIO DE TRABAJOS (PE y SE) ----------
# Construye una lista de trabajos a partir del contenido de READS_DIR.
# Cada entrada tiene: LABEL, TYPE (pe_short|se_short), R1, R2(opt)
# ---------- INVENTARIO DE TRABAJOS (PE y SE) ----------
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

  label="$(basename "${base}")"   # p.ej., HepC_ERR118961
  [[ -n "${SEEN[$label]:-}" ]] && continue
  SEEN[$label]=1

  LABELS+=("$label")
  TYPES+=("pe_short")
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
  TYPES+=("se_short")
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

# Deriva SAMPLE solo para logging (parte después del prefijo y guion bajo)
SAMPLE="${LABEL#${PREFIX}_}"

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
} >> "${METRICFILE}"

# ---------- EJECUCIÓN + PERFILADO ----------
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
      2>>"${METRICFILE}"
else
  "${TIMECMD}" -v \
  srun -n 1 -c "${SLURM_CPUS_PER_TASK}" --cpu-bind=cores --mem-bind=local \
  virus2tree_step2.sh \
      --read_type se_short \
      --reads "$R1" \
      --threads "${SLURM_CPUS_PER_TASK}" \
      --root_dir "${OUT_ROOT}" \
      --temp_dir "${TMPDIR}" \
      --out_dir r2t_ref \
      2>>"${METRICFILE}"
fi
rc=$?
set -e

echo "ExitCode_cmd  : ${rc}" >> "${METRICFILE}"

# ---------- MÉTRICAS SLURM (post-job) ----------
sleep 30
sacct -j "${SLURM_JOB_ID}" \
      --format=JobIDRaw,JobName,State,ExitCode,Elapsed,TotalCPU,CPUTime,MaxRSS,MaxVMSize,AveRSS,AveCPU \
      >> "${METRICFILE}"

exit "$rc"
