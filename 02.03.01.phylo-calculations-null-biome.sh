#!/bin/bash

#SBATCH --job-name=PD-null-biome
#SBATCH --chdir=/work/haehn
#SBATCH --output=/work/%u/%x-%A-%a.log
#SBATCH --error=/work/%u/%x-%A-%a.err
#SBATCH --time=3-00:00:00
#SBATCH --mem-per-cpu=16G

module load foss R/4.2.2

input_dir=/data/splot/data_gh
output_dir=/work/$USER/PD-null-biome
mkdir -p "$output_dir"
output="$output_dir/$SLURM_ARRAY_TASK_ID.Rds"
x="$SLURM_ARRAY_TASK_MIN"
y="$SLURM_ARRAY_TASK_MAX"
start=$SLURM_ARRAY_TASK_ID
end=$((SLURM_ARRAY_TASK_ID+SLURM_ARRAY_TASK_STEP-1))


Rscript --vanilla $HOME/sPlot-project-FD-PD/01.R-scripts/02.03.01.phylo-calculations-null-biome.R \
"$input_dir" \
"$output" \
"$start" \
"$end"
