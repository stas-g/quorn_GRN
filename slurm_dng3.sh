#!/bin/bash

#SBATCH --array=1-50

mycommand.exe input_file_$SLURM_ARRAY_TASK_ID
