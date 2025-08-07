#!/bin/bash
#SBATCH --job-name=v2t_step1_cmv            
#SBATCH --account=proj-gm0001 
#SBATCH --partition=medium
#SBATCH --output=/stornext/snfs170/next-gen/scratch/adrian/projects/r2t/results_cmv/v2t_cmv_%j.out    # STDOUT
#SBATCH --time=3-00:00:00        # max time (días‑hh:mm:ss)
#SBATCH --cpus-per-task=24     # threads
#SBATCH --mem=0                  # “0” = all the mem of the node
#SBATCH --mail-type=END,FAIL     
#SBATCH --mail-user=adrianalexandro777@gmail.com

source ~/.bashrc          
conda activate my_env        
#cd /stornext/snfs170/next-gen/scratch/adrian/projects/r2t/results_cmv

LOGDIR=/stornext/snfs170/next-gen/scratch/adrian/projects/r2t/results_cmv
mkdir -p "${LOGDIR}"
METRICFILE=${LOGDIR}/v2t_cmv_${SLURM_JOB_ID}.metrics
TIMECMD=/stornext/snfs170/next-gen/scratch/adrian/local/bin/time   # ruta fija

"$TIMECMD" -v \
srun --cpu-bind=none --mem-bind=none \
virus2tree_step1.sh \
    -i ../data/cmv.csv \
    -g ../data/cmvo.csv \
    --threads ${SLURM_CPUS_PER_TASK} \
    --temp_dir /space1/tmp/r2t_cmv9 \
    --out_dir r2t_ref \
    --root_dir ../results_cmv \
    --resume_download \
    2>>"${METRICFILE}"

sleep 30
sacct -j "${SLURM_JOB_ID}" \
      --format=JobIDRaw,JobName,State,ExitCode,Elapsed,TotalCPU,CPUTime,MaxRSS,MaxVMSize,AveRSS,AveCPU \
      >> "${METRICFILE}"
