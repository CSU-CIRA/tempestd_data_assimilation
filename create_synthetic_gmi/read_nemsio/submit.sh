#!/bin/sh
#=======================================================
## Below are SLURM (S4) commands
#SBATCH --job-name=read_fv3gfs_nemsio
#SBATCH --partition=service
#SBATCH --account=nesdis-rdo1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=24000
#SBATCH --time=00:20:00
#SBATCH --output=read_fv3gfs_nemsio.%j.out
#SBATCH --error=read_fv3gfs_nemsio.%j.err

make -f makefile.read_fv3gfs_atm
make -f makefile.read_fv3gfs_sfc

set rtype='anl'
#set rtype='f000'

srun --cpu_bind=core --distribution=block:block --ntasks=1 ./read_fv3gfs_atm_nemsio gfs.t12z.atm${rtype}.nemsio
srun --cpu_bind=core --distribution=block:block --ntasks=1 ./read_fv3gfs_sfc_nemsio gfs.t12z.sfc${rtype}.nemsio
