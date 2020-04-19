#!/usr/bin/bash
#SBATCH --partition=EPIC2
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:2
#SBATCH --job-name="sph_dam_break"
#SBATCH --output=SPH_dam_break.%j.out

cd ${SLURM_SUBMIT_DIR}
. environment.sh idun

export SCALE=1.0
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

make clean
make SCALE=${SCALE}

# If 'plot' doesn't exist, create it
mkdir -p plot
# If 'plot' is full of old files, clean it
rm -f plot/*

# Run the simulation
date
time srun --mpi=pmi2 ./sph
date
echo "SCALE=${SCALE}"
