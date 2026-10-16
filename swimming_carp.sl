#!/bin/bash -l
#SBATCH --job-name=CARP_snkm

# ---------------------------------- Account --------------------------------- #
#SBATCH --account=<your_account>

# ----------------------------- Request resources ---------------------------- #
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=1-00:00:00
#SBATCH --mem=64GB

# ----------------------------------- Queue ---------------------------------- #
#SBATCH -p batch

# ------------------------------- Notifications ------------------------------ #
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<your_email>


echo ==========Running on `hostname`==========
# Date mark
echo ==========Start time is `date`==========

# ---------------------------------------------------------------------------- #
#                              Environment setting                             #
# ---------------------------------------------------------------------------- #
# ---------------------------------- Module ---------------------------------- #
module load Singularity
module list

# ----------------------------------- Conda ---------------------------------- #
# source ~/miniconda3/etc/profile.d/conda.sh
conda activate snkm9
echo "Active Conda environment: $CONDA_DEFAULT_ENV"
snakemake --version

# -------------------------- Working directory setup ------------------------- #

BASE_DIR=<your/path>/carp_snkm
SCRIPT=${BASE_DIR}/Snakefile
CONFIG=${BASE_DIR}/config.yaml
OUTPUT_DIR=${BASE_DIR}/output
mkdir -p $OUTPUT_DIR

echo Working directory is $BASE_DIR
echo Running script is $SCRIPT
echo Output dir is $OUTPUT_DIR


cd $BASE_DIR
snakemake --snakefile $SCRIPT --configfile $CONFIG --cores ${SLURM_CPUS_PER_TASK} --use-singularity --singularity-args "--bind ${TMPDIR}:${TMPDIR}" --printshellcmds

echo
echo "Finished: carp managed to swim to shore at $(date)"