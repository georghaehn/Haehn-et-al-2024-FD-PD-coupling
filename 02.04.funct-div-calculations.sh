#!/bin/bash
 
#SBATCH --job-name=FD-calculation
#SBATCH --chdir=/work/haehn
#SBATCH --output=/work/%u/%x-%A-%a.log
#SBATCH --error=/work/%u/%x-%A-%a.err
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=16G

module load foss/2020b R/4.0.4-2

input_dir=/data/splot/data_gh
output_dir=/work/$USER/FD-calculation
mkdir -p "$output_dir"
output="$output_dir/$SLURM_ARRAY_TASK_ID.Rds"
x="$SLURM_ARRAY_TASK_MIN"
y="$SLURM_ARRAY_TASK_MAX"
start=$SLURM_ARRAY_TASK_ID
end=$((SLURM_ARRAY_TASK_ID+SLURM_ARRAY_TASK_STEP-1))


Rscript --vanilla $HOME/sPlot-master-thesis/01.scripts/02.04.funct-div-calculations.R \
  "$input_dir" \
  "$output" \
  "$start" \
  "$end"
