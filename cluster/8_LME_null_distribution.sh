#!/bin/bash

#SBATCH --time=1-12   # walltime 
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks) 
#SBATCH --nodes=1   # number of nodes 
#SBATCH --mem-per-cpu=4G  # memory per CPU core 
#SBATCH --account=def-gmitsis  
#SBATCH --job-name=LME-MY_I-MY_MEAS-MY_PIPEL
#SBATCH --mail-user=mary.miedema@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --output=out_LME-MY_I_MY_MEAS_MY_PIPEL.txt

cd /home/miedemam/scratch/MGH_derivatives
input_root="/home/miedemam/scratch/MGH_derivatives/LME_model/inputs/"
output_root="/home/miedemam/scratch/MGH_derivatives/LME_model/outputs/"

# choose parameters
task_ses=MY_I
measure=MY_MEAS
pipel=MY_PIPEL

echo "Estimating significance for LME model for:"
echo $task_ses
echo $pipel
echo $measure

module load matlab/2023b.2

echo "Matlab loaded! Running null model estimation script:"

matlab -nodisplay -r "LME_fit_null_distribution('${measure}','${pipel}',${task_ses},2000,'${input_root}','${output_root}')" -sd /home/miedemam/scratch/MGH_derivatives/matlab_functions

echo "Finished estimating LME significance for all connections!"
echo "This took:"
echo $SECONDS
