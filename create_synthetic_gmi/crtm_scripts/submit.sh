#!/bin/sh
#=======================================================
## Below are SLURM (S4) commands
#SBATCH --job-name=crtm_gmi_control
#SBATCH --partition=service
#SBATCH --account=nesdis-rdo1
#SBATCH --ntasks=1
#SBATCH --time=08:00:00
#SBATCH --mem-per-cpu=24000
#SBATCH --output=crtm_gmi_control.%j.out
#SBATCH --error=crtm_gmi_control.%j.err

srun ./crtm_gmi_control.sh
