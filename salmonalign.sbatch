#!/bin/bash
##
## example-array.slurm.sh: submit an array of jobs with a varying parameter
##
## Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
## All other lines are read by the shell.
##
#SBATCH --account=priority-ryanthum        #specify the account to use
#SBATCH --job-name=salmonalign    # job name
#SBATCH --partition=priority              # queue partition to run the job in
#SBATCH --nodes=1                       # number of nodes to allocate
#SBATCH --ntasks-per-node=1             # number of descrete tasks - keep at one except for MPI
#SBATCH --cpus-per-task=50              # number of cores to allocate
#SBATCH --mem=200G                     # 2000 MB of Memory allocated; set --mem with care
#SBATCH --time=0-24:00:00                 # Maximum job run time
##SBATCH --array=1-3                  # Number of jobs in array
#SBATCH --output=salmonalign-%j.out
#SBATCH --error=salmonalign-%j.err
#SBATCH --mail-user=gregory.chorak@student.montana.edu
#SBATCH --mail-type=ALL

## alignmnet and quantification of cleaned DE reads using Salmon default settings. 


cd ./rawDE/

for i in `ls -1 *_1_AQtrim.fq | sed 's/_1_AQtrim.fq//'`

do

./../programs/salmon-1.6.0/bin/salmon quant -i ./../milfoil_index -l A \
         -1 $i\_1_AQtrim.fq \
         -2 $i\_2_AQtrim.fq \
         -p 50 --validateMappings -o quants/$i\_quant

done
