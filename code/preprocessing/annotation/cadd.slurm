#!/bin/bash -l
#SBATCH --job-name=CADD_AFR
#SBATCH --time=24:0:0
#SBATCH --partition=lrgmem
#SBATCH --cpus-per-task=4
#SBATCH --array=1-22

#### load and unload modules you may need

ml anaconda

#### execute code to get CADD annotations

datadir=/scratch/groups/abattle4/victor/WatershedAFR/data/
envloc=~/.conda/envs/cadd
cadd_loc=/scratch/groups/abattle4/jessica/RareVar_AFR/cadd/CADD-scripts
code_loc=/scratch/groups/abattle4/jessica/RareVar_AFR/WatershedAFR/code/preprocessing/annotation/

 ./cadd_slurm.sh $datadir $envloc $cadd_loc $code_loc
echo "Finished with job $SLURM_JOBID"
