#!/bin/bash

######## Slurm options ########

#### general settings
#SBATCH --job-name=spinglass_test
#SBATCH --output=./logs/slurm-%j.out
#SBATCH --error=./logs/slurm_error-%j.out
#_BATCH --mail-type=ALL #other good setting: END,FAIL
#SBATCH --mail-user=hartmut.schoon@uni-oldenburg.de


#### request resources
#SBATCH --partition=carl.p
#SBATCH --time=0:05:00	#max runtime
#_BATCH --nodes 2
#_BATCH --cpus-per-task=1
#_BATCH --mem=1G
#_BATCH --mem-per-cpu=1G

#### parallel 
#SBATCH --ntasks=10
#_BATCH --array1-10:1%4

######## Slurm options ########

## Path and args
results_path="$WORK/spinglass"
command_args=("no_ui" "path=$results_path/results_{}")
parallel_args=("-N1" "-j$SLURM_NTASKS" "--joblog" "./logs/parallel$SLURM_JOB_ID.log" "--delay" "1")

export RUST_BACKTRACE=1

parallel "${parallel_args[@]}" srun ./target/release/spinglass "${command_args[@]}" ::: {1..10}