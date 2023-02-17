#!/bin/bash

######## Slurm options ########

#### general settings
#SBATCH --job-name=spinglass
#SBATCH --output=./slurm_out/slurm-%j.out

#### request resources
#SBATCH --partition=carl.p
#SBATCH --time=0:05:00	#max runtime

######## Slurm options ########


RESULTS_PATH="$HOME/spinglass/test_results/"
./target/release/spinglass no_ui "path=$RESULTS_PATH"
