#!/bin/bash
#SBATCH --job-name=v2t_step1_rsv             # Nombre del job
#SBATCH --account=proj-gm0001 
#SBATCH --partition=medium
#SBATCH --output=/stornext/snfs170/next-gen/scratch/adrian/projects/r2t/results_reads_r2t_test/sra_reads_%j.out    # STDOUT
#SBATCH --time=3-00:00:00        # max time (días‑hh:mm:ss)
#SBATCH --cpus-per-task=12    # threads
#SBATCH --mem=0                  # “0” = all the mem of the node
#SBATCH --mail-type=END,FAIL     
#SBATCH --mail-user=adrianalexandro777@gmail.com

source ~/.bashrc          
conda activate my_env        

LOGDIR=/stornext/snfs170/next-gen/scratch/adrian/projects/r2t/results_reads_r2t_test
mkdir -p ${LOGDIR}
METRICFILE=${LOGDIR}/sra_reads_${SLURM_JOB_ID}.metrics
TIMECMD=/stornext/snfs170/next-gen/scratch/adrian/local/bin/time   # ruta fija

"$TIMECMD" -v \
srun --cpu-bind=none --mem-bind=none \
download_sra_reads_copy.sh -i ../data/sra_runs_samples5.csv \
-o ../results_reads_r2t_test -d \
    2>>"${METRICFILE}"
sleep 30
sacct -j "${SLURM_JOB_ID}" \
      --format=JobIDRaw,JobName,State,ExitCode,Elapsed,TotalCPU,CPUTime,MaxRSS,MaxVMSize,AveRSS,AveCPU \
      >> "${METRICFILE}"

