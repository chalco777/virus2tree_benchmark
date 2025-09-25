#!/bin/bash
#SBATCH --job-name=step3_cmv
#SBATCH --account=proj-gm0001
#SBATCH --partition=long
#SBATCH --output=cmv_step3_%j.out
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=adrianalexandro777@gmail.com
#SBATCH --export=ALL         

source ~/.bashrc
conda activate my_env

# ---------- LOGS ----------
METRICFILE=cmv_step3_${SLURM_JOB_ID}.metrics

TIMECMD=/stornext/snfs170/next-gen/scratch/adrian/local/bin/time   # ruta fija
R2T_DIR=/stornext/snfs170/next-gen/scratch/adrian/projects/r2t/results_cmv_step2

"$TIMECMD" -v \
read2tree --step 3combine \
--standalone_path ${R2T_DIR}/marker_genes \
--dna_reference ${R2T_DIR}/dna_ref.fa \
--output_path ${R2T_DIR}/r2t_ref --tree --debug \
2>>"${METRICFILE}"

sleep 30
sacct -j "${SLURM_JOB_ID}" \
      --format=JobIDRaw,JobName,State,ExitCode,Elapsed,TotalCPU,CPUTime,MaxRSS,MaxVMSize,AveRSS,AveCPU \
      >> "${METRICFILE}"
