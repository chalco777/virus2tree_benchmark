#!/bin/bash
#SBATCH --job-name=step3_hepc
#SBATCH --account=proj-gm0001
#SBATCH --partition=medium
#SBATCH --output=/stornext/snfs170/next-gen/scratch/adrian/projects/r2t/results_hepc_test2/v2t_step3_${SLURM_JOB_ID}.out
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
LOGDIR=/stornext/snfs170/next-gen/scratch/adrian/projects/r2t/results_hepc_test2
mkdir -p "${LOGDIR}"
METRICFILE=${LOGDIR}/v2t_step3_${SLURM_JOB_ID}.metrics

# ---------- EJECUCIÓN + PERFILADO ----------
# srun hereda los --cpus-per-task del job; /usr/bin/time -v captura ~30 métricas
TIMECMD=/stornext/snfs170/next-gen/scratch/adrian/local/bin/time   # ruta fija
R2T_DIR=/stornext/snfs170/next-gen/scratch/adrian/projects/r2t/results_hepc_test2

"$TIMECMD" -v \
srun -n 1 -c "${SLURM_CPUS_PER_TASK}" --cpu-bind=cores --mem-bind=local \
read2tree --step 3combine \
--standalone_path ${R2T_DIR}/marker_genes \
--dna_reference ${R2T_DIR}/dna_ref.fa \
--output_path ${R2T_DIR}/r2t_ref --tree --debug \
2>>"${METRICFILE}"

# ---------- MÉTRICAS SLURM (post‑job) ----------
# Espera un momento para que sacct tenga los datos finalizados
sleep 30
sacct -j "${SLURM_JOB_ID}" \
      --format=JobIDRaw,JobName,State,ExitCode,Elapsed,TotalCPU,CPUTime,MaxRSS,MaxVMSize,AveRSS,AveCPU \
      >> "${METRICFILE}"
