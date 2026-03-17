#!/bin/bash

#SBATCH --array=1-25
#SBATCH --time=6:00:00   # walltime 
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks) 
#SBATCH --nodes=1   # number of nodes 
#SBATCH --mem-per-cpu=128G  # memory per CPU core 
#SBATCH --account=def-gmitsis  
#SBATCH --job-name=merge_MNIs
#SBATCH --mail-user=mary.miedema@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --output=out_MNI_merge_ses02_%j.txt

module load StdEnv/2023  
module load fsl/6.0.7.7 

# get the subject number
cd /home/miedemam/scratch/MGH_derivatives
sub_num=$(sed -n "${SLURM_ARRAY_TASK_ID}p" run_subs_25.txt) 

echo "Processing subject:"
echo $sub_num

# specify the pipeline
pipel="wCompCor"
echo $pipel

# set up the directory structure
mkdir ./central/MNI_preproc/${pipel}/sub-${sub_num}

# masks
brainstem_mask=./brainstem_box.nii
cortex_mask=./non_brainstem_box.nii

# choose session and task
declare -a tasks=("rest" "breathing" "coldpressor" "rest" "breathing")  
declare -a sess=("02" "02" "02" "03" "03")
for j in $(seq 3 4)
do
task=${tasks[$j]}
ses_num=${sess[$j]}
out_dir=./central/MNI_preproc/${pipel}/sub-${sub_num}/ses-${ses_num}/func
mkdir ./central/MNI_preproc/${pipel}/sub-${sub_num}/ses-${ses_num}
mkdir ./central/MNI_preproc/${pipel}/sub-${sub_num}/ses-${ses_num}/func
bids_root=sub-${sub_num}_ses-${ses_num}_task-${task}
out_func=$out_dir/${bids_root}_space-MNI152merged_desc-preproc_bold.nii.gz
echo "Now processing: $bids_root"
cortex_func=./central/MNI_func/sub-${sub_num}/ses-${ses_num}/${task}/${pipel}/${bids_root}_MNI.nii.gz
cortex_masked=$out_dir/cortex_masked.nii.gz
brainstem_func=./central/MNI_func/sub-${sub_num}/ses-${ses_num}/${task}/${pipel}/${bids_root}_MNI_brainstem.nii.gz
brainstem_masked=$out_dir/brainstem_masked.nii.gz
# now merge the two files together
fslmaths $cortex_func -mas $cortex_mask $cortex_masked
fslmaths $brainstem_func -mas $brainstem_mask $brainstem_masked
fslmaths $cortex_masked -add $brainstem_masked $out_func
rm $brainstem_masked
rm $cortex_masked
echo "Now finished processing: $bids_root"
done
echo "Completely finished! This took:"
echo $SECONDS
