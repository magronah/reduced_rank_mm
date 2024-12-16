#!/bin/bash
# submit_job.sh

#SBATCH --account=def-bolker

# Job parameters
#SBATCH --job-name=my_job
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G

#SBATCH --array=5,43,66,67,169,170,181,184,193,210,253,292,301,450,451,452,457,473,490
#SBATCH --output=array-%j-%a.out               # log file
#SBATCH --error=array-%j-%a.err                # error file

#SBATCH --mail-user=agronahm@mcmaster.ca    # who to email
#SBATCH --mail-type=FAIL                  # when to email


# Load R module
module load r/4.4.0

# Run R script with array index as argument

Rscript rrzi.R $SLURM_ARRAY_TASK_ID

