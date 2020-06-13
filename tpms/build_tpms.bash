#!/bin/bash

# To set a project account for credit charging,
#SBATCH --account=qb48


# Request CPU resource for a serial job
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
# SBATCH --exclusive
#SBATCH --cpus-per-task=1

# Memory usage (MB)
#SBATCH --mem-per-cpu=4000

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=0-03:00:00


# To receive an email when job completes or fails
#SBATCH --mail-user=angus.leung1@monash.edu
#SBATCH --mail-type=END
# SBATCH --mail-type=FAIL


# Set the file for output (stdout)
#SBATCH --output=logs/build_tpms.out

# Set the file for error log (stderr)
#SBATCH --error=logs/build_tpms.err


# Use reserved node to run job when a node reservation is made for you already
# SBATCH --reservation=reservation_name

# Command to run a serial job
module load matlab/r2018a

matlab -nodisplay -nodesktop -r "build_tpms_perSong; exit"
