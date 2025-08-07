#!/bin/bash
#SBATCH --job-name=v2t_step1_rsv             # Nombre del job
#SBATCH --account=proj-gm0001 
#SBATCH --partition=medium
#SBATCH --output=/stornext/snfs170/next-gen/scratch/adrian/projects/r2t/results_rsv/v2t_step1_%j.out    # STDOUT
#SBATCH --time=3-00:00:00        # max time (días‑hh:mm:ss)
#SBATCH --cpus-per-task=24    # threads
#SBATCH --mem=0                  # “0” = all the mem of the node
#SBATCH --mail-type=END,FAIL     
#SBATCH --mail-user=adrianalexandro777@gmail.com

source ~/.bashrc          
conda activate my_env        

LOGDIR=/stornext/snfs170/next-gen/scratch/adrian/projects/r2t/results_rsv
mkdir -p "${LOGDIR}"
METRICFILE=${LOGDIR}/v2t_step1_${SLURM_JOB_ID}.metrics
TIMECMD=/stornext/snfs170/next-gen/scratch/adrian/local/bin/time   # ruta fija

"$TIMECMD" -v \
srun --cpu-bind=none --mem-bind=none \
virus2tree_step1.sh \
    -i ../data/rsv.csv \
    -g ../data/rsvo.csv \
    --threads ${SLURM_CPUS_PER_TASK} \
    --temp_dir /space1/tmp/r2t_rsv7 \
    --root_dir ../results_rsv \
    --out_dir r2t_ref \
    2>>"${METRICFILE}"

sleep 30
sacct -j "${SLURM_JOB_ID}" \
      --format=JobIDRaw,JobName,State,ExitCode,Elapsed,TotalCPU,CPUTime,MaxRSS,MaxVMSize,AveRSS,AveCPU \
      >> "${METRICFILE}"

