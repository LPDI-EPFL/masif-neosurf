#!/bin/bash
if [ -f "./prep/batch_submitter.slurm" ]; then
    sbatch ./prep/batch_submitter.slurm
else
    echo "Submission file does not exist. Run preparation script first."
fi
