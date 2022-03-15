#!/bin/bash
##
## example-array.slurm.sh: submit an array of jobs with a varying parameter
##
## Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
## All other lines are read by the shell.
##
#SBATCH --account=priority-ryanthum        #specify the account to use
#SBATCH --job-name=bbtrim    # job name
#SBATCH --partition=priority              # queue partition to run the job in
#SBATCH --nodes=1                       # number of nodes to allocate
#SBATCH --ntasks-per-node=1             # number of descrete tasks - keep at one except for MPI
#SBATCH --cpus-per-task=30              # number of cores to allocate
#SBATCH --mem=100G                     # 2000 MB of Memory allocated; set --mem with care
#SBATCH --time=0-24:00:00                 # Maximum job run time
##SBATCH --array=1-3                  # Number of jobs in array
#SBATCH --output=bbtrim-%j.out
#SBATCH --error=bbtrim-%j.err
#SBATCH --mail-user=gregory.chorak@student.montana.edu
#SBATCH --mail-type=ALL

## adapter trimming using bbduk default settings. 


cd ./rawDE/

for i in `ls -1 *_1.fq.gz | sed 's/_1.fq.gz//'`

do

./../programs/bbmap/bbduk.sh -Xmx100g t=30 in1=$i\_1.fq.gz in2=$i\_2.fq.gz out1=$i\R_1_ARduked.fq out2=$i\R_2_ARduked.fq ref=./../programs/bbmap/resources/adapters.fa ktrim=r mink=11 k=23 hdist=1 tpe tbo


done
