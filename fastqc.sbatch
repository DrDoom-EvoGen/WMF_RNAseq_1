#!/bin/bash
##
## example-array.slurm.sh: submit an array of jobs with a varying parameter
##
## Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
## All other lines are read by the shell.
##
#SBATCH --account=priority-ryanthum        #specify the account to use
#SBATCH --job-name=fastqc    # job name
#SBATCH --partition=priority              # queue partition to run the job in
#SBATCH --nodes=1                       # number of nodes to allocate
#SBATCH --ntasks-per-node=1             # number of descrete tasks - keep at one except for MPI
#SBATCH --cpus-per-task=30              # number of cores to allocate
#SBATCH --mem=100G                     # 2000 MB of Memory allocated; set --mem with care
#SBATCH --time=0-24:00:00                 # Maximum job run time
##SBATCH --array=1-3                  # Number of jobs in array
#SBATCH --output=fastqc-%j.out
#SBATCH --error=fastqc-%j.err
#SBATCH --mail-user=gregory.chorak@student.montana.edu
#SBATCH --mail-type=ALL

## checking quality of reads after adapter removal and quality trimming. 


cd ./rawDE/


./../programs/FastQC/fastqc *AQtrim.fq

