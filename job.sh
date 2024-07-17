#!/bin/bash
#SBATCH -p gpu_a100_8
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --gres=gpu:1
#SBATCH --qos=qos_gpu_a100_priority
#SBATCH --mem=80G
#SBATCH -t 0-23:00 # time (D-HH:MM)
#SBATCH --job-name="regent-test"
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH --mail-user=f20210781@hyderabad.bits-pilani.ac.in
#SBATCH --mail-type=ALL

apptainer exec --nv --env LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/cuda-12/compat:/legion/bindings/regent --env CPATH=\$CPATH:/usr/local/cuda-12/targets/x86_64-linux/include/ --bind ./data:/data --bind ./build:/build /home/anil/soumitra/regent_env/v24.03.0/v24.03.0.sif /build/meshfree_solver_test.out -ll:csize 77G
