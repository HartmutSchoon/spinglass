#!/bin/bash

#Script to start the spinglass simulation

#Read input argument. This should be a integer provided by parallel to specify the task id. 
TASK_ID=$1

#In this path the results are stored
results_path="$WORK/spinglass_L4/results_$TASK_ID"

#Start spinglass simulation
echo "Starting Spinglass with ID $TASK_ID on $(hostname)"
./target/release/spinglass no_ui path="$results_path"