#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 8G
#SBATCH --time 1:00:00
#SBATCH --array=11-1000
echo STARTING AT `date`
source /work/upcorreia/bin/load_masif_environment.sh
./make_transformations_9A.sh ${SLURM_ARRAY_TASK_ID}
echo FINISHED at `date`
