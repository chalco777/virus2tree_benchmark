#!/bin/bash
#SBATCH --job-name=step2_hepc
#SBATCH --account=proj-gm0001
#SBATCH --array=0-2%1 #one job at a time
#SBATCH --partition=medium
#SBATCH --output=/stornext/snfs170/next-gen/scratch/adrian/projects/r2t/results_hepc_test/hepc_test_%A_%a.out
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=0
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=adrianalexandro777@gmail.com
#SBATCH --export=ALL          # pasa las vars de entorno a sacct/sstat

# ---------- ENTORNO ----------
source ~/.bashrc
conda activate my_env
#cd /stornext/snfs170/next-gen/scratch/adrian/projects/r2t/results_hepc

# ---------- RUTA LOGS ----------
LOGDIR=/stornext/snfs170/next-gen/scratch/adrian/projects/r2t/results_hepc_test
mkdir -p "${LOGDIR}"
TASK=${SLURM_ARRAY_TASK_ID:-0}
TMPDIR=/space1/tmp/r2t_rsv7_${SLURM_JOB_ID}_${TASK}
METRICFILE=${LOGDIR}/v2t_step2_${SLURM_JOB_ID}_${TASK}.metrics

# ---------- EJECUCIÓN + PERFILADO ----------
# srun hereda los --cpus-per-task del job; /usr/bin/time -v captura ~30 métricas
TIMECMD=/stornext/snfs170/next-gen/scratch/adrian/local/bin/time   # ruta fija
PROJECT=/stornext/snfs170/next-gen/scratch/adrian/projects/r2t
READS_DIR=${PROJECT}/results_reads_hepc
PREFIX=HepC
SAMPLE_IDS=(ERR118961 ERR118962 ERR118963)
SAMPLE=${SAMPLE_IDS[$SLURM_ARRAY_TASK_ID]}

FASTQ_R1="${READS_DIR}/${PREFIX}_${SAMPLE}_1.fastq"
FASTQ_R2="${READS_DIR}/${PREFIX}_${SAMPLE}_2.fastq"


"$TIMECMD" -v \
srun -n 1 -c "${SLURM_CPUS_PER_TASK}" --cpu-bind=cores --mem-bind=local \
virus2tree_step2.sh \
    --read_type pe_short --reads "${FASTQ_R1}" "${FASTQ_R2}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --root_dir "../results_hepc_test" \
    --temp_dir $TMPDIR \
    --out_dir r2t_ref \
    2>>"${METRICFILE}"

# ---------- MÉTRICAS SLURM (post‑job) ----------
# Espera un momento para que sacct tenga los datos finalizados
sleep 30
sacct -j "${SLURM_JOB_ID}" \
      --format=JobIDRaw,JobName,State,ExitCode,Elapsed,TotalCPU,CPUTime,MaxRSS,MaxVMSize,AveRSS,AveCPU \
      >> "${METRICFILE}"
