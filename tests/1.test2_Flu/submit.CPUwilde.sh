#!/bin/bash
#SBATCH -J jobName
#SBATCH --partition=CPU
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1

source /etc/profile.d/modules.sh

srun hostname -s | sort -u >slurm.hosts
/scratch1/oliver/bin/WiggleToSite/bin/wiggleToSite input.txt moving.pdb overlap.pdb target.pdb > output.txt 
