#!/bin/bash

#SBATCH --job-name=BRT-calculation
#SBATCH --chdir=/work/haehn
#SBATCH --output=/work/%u/%x-%A-%a.log
#SBATCH --error=/work/%u/%x-%A-%a.err
#SBATCH --time=16-00:00:00
#SBATCH --mem-per-cpu=16G

module load foss/2020b R/4.0.4-2

input_dir=/data/splot/data_gh
output_dir=/work/$USER/BRT-calculation
mkdir -p "$output_dir"
x="$SLURM_ARRAY_TASK_MIN"
y="$SLURM_ARRAY_TASK_MAX"
start=$SLURM_ARRAY_TASK_ID


Rscript --vanilla $HOME/sPlot-project-FD-PD/01.R-scripts/04.02.BRT.R \
  "$input_dir" \
  "$output_dir" \
  "$start"
