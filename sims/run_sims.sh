#!/bin/bash  
#SBATCH -J stream_sims -p sched_mit_lnecib -N 1 --ntasks 64 --cpus-per-task 1 -t 96:00:00
python run_sims.py 0  1>run_sims.out 2>run_sims.err
