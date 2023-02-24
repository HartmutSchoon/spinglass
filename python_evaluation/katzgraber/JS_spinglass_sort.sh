#!/bin/bash

######## Slurm options ########

#### general settings
#SBATCH --job-name=spinglass_sort
#SBATCH --output=./logs/slurm-%j.out
#SBATCH --error=./logs/slurm_error-%j.out
#SBATCH --mail-type=ALL #other good setting: END,FAIL
#SBATCH --mail-user=hartmut.schoon@uni-oldenburg.de


#### request resources
#SBATCH --partition=carl.p
#_BATCH --time=0:05:00	#max runtime
#_BATCH --nodes 2
#_BATCH --cpus-per-task=1
#_BATCH --mem=1G
#_BATCH --mem-per-cpu=1G

#### parallel 
#SBATCH --ntasks=22
#_BATCH --array1-10:1%4

######## Slurm options ########


#Set srun command. Every srun only needs one node but there are $SLURM_NTASKS srun commands running in parallel
srun="srun -n1"
#Set parallel command. One input argument, delay of 0.2 between single execution, j jobs run in paralel and 
#log is saved to logs folder
parallel="parallel -N 1 --delay 0.2 -j $SLURM_NTASKS --joblog ./logs/parallel_$SLURM_JOB_ID.log -a T_list.txt"


$parallel "$srun ./main_new.py {}" 
