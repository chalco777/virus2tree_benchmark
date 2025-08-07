#!/bin/bash
#SBATCH --job-name=v2t_step1_hepc
#SBATCH --account=proj-gm0001
#SBATCH --partition=medium
#SBATCH --output=/stornext/snfs170/next-gen/scratch/adrian/projects/r2t/results_hepc_test/hepc_test_%j.out
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
METRICFILE=${LOGDIR}/v2t_step1_${SLURM_JOB_ID}.metrics

# ---------- EJECUCIÓN + PERFILADO ----------
# srun hereda los --cpus-per-task del job; /usr/bin/time -v captura ~30 métricas
TIMECMD=/stornext/snfs170/next-gen/scratch/adrian/local/bin/time   # ruta fija

"$TIMECMD" -v \
srun --cpu-bind=none --mem-bind=none \
     virus2tree_step1.sh \
         -i ../data/hepc_test.csv \
         -g ../data/hepco.csv \
         --use_only_mat_peptides \
         --threads ${SLURM_CPUS_PER_TASK} \
         --root_dir ../results_hepc_test \
         --temp_dir /space1/tmp/r2t_hepc6 \
         --out_dir  r2t_ref \
    2>>"${METRICFILE}"

# ---------- MÉTRICAS SLURM (post‑job) ----------
# Espera un momento para que sacct tenga los datos finalizados
sleep 30
sacct -j "${SLURM_JOB_ID}" \
      --format=JobIDRaw,JobName,State,ExitCode,Elapsed,TotalCPU,CPUTime,MaxRSS,MaxVMSize,AveRSS,AveCPU \
      >> "${METRICFILE}"
