#!/bin/bash
# submit_job.sh

#SBATCH --account=def-bolker

# Job parameters
#SBATCH --job-name=my_job
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G

#SBATCH --array=1-7
#SBATCH --output=array-%j-%a.out               # log file
#SBATCH --error=array-%j-%a.err                # error file

#SBATCH --mail-user=agronahm@mcmaster.ca    # who to email
#SBATCH --mail-type=FAIL                  # when to email


# Load R module
module load r/4.3.1

# Run R script with array index as argument

Rscript zero_inflation_est.R $SLURM_ARRAY_TASK_ID
