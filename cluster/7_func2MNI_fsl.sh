#!/bin/bash

#SBATCH --array=1-25
#SBATCH --time=3:30:00   # walltime 
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks) 
#SBATCH --nodes=1   # number of nodes
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=64G  # memory per CPU core 
#SBATCH --account=def-gmitsis  
#SBATCH --job-name=t_MYSES-MYTASK-MYPIP
#SBATCH --mail-user=mary.miedema@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --output=out_t_MYSES-MYTASK-MYPIP_%j.txt

module load StdEnv/2023  
module load fsl/6.0.7.7 

# choose session and task
ses_num=MYSES
task=MYTASK
pipel=MYPIP

# get the subject number
cd /home/miedemam/scratch/MGH_derivatives
sub_num=$(sed -n "${SLURM_ARRAY_TASK_ID}p" run_subs_25.txt) 

echo "Processing subject:"
echo $sub_num

# where files will be saved
out_dir=./central/MNI_func/sub-${sub_num}/ses-${ses_num}/${task}/$pipel
mkdir ./central/MNI_func/sub-${sub_num}
mkdir ./central/MNI_func/sub-${sub_num}/ses-${ses_num}
mkdir ./central/MNI_func/sub-${sub_num}/ses-${ses_num}/${task}
mkdir ./central/MNI_func/sub-${sub_num}/ses-${ses_num}/${task}/$pipel

bids_root=sub-${sub_num}_ses-${ses_num}_task-${task}

# transform the functional file into a standard space using the FSL FEAT transformations
in_func=./central/func_${pipel}/${bids_root}_func.nii.gz
in_reg_root=./feat/sub-${sub_num}/ses-${ses_num}/${task}.feat/reg
# only the session 2 resting scan contains the warp to MNI152 space for each subject
highres2standard_warp=./feat/sub-${sub_num}/ses-02/rest.feat/reg/highres2standard_warp.nii.gz
if [ "$ses_num" = "02" ]; then
func2highres_premat=${in_reg_root}/example_func2highres.mat
fi
if [ "$ses_num" = "03" ]; then
# in this case concatenate the necessary transformations
anat_folder=/home/miedemam/scratch/MGH_derivatives/anat
convert_xfm -omat ${in_reg_root}/example_func2highres2.mat -concat  ${anat_folder}/trans_ses03_2_ses02_sub-${sub_num}.mat ${in_reg_root}/example_func2highres.mat
func2highres_premat=${in_reg_root}/example_func2highres2.mat
fi

applywarp -i $in_func -r /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/CUDA/gcc9/cuda11.0/fsl/6.0.4/fsl/data/standard/MNI152_T1_1mm_brain  -o ${out_dir}/${bids_root}_MNI.nii.gz -w $highres2standard_warp --premat=$func2highres_premat

echo "Finished transforming the functional scan!"
echo "This took:"
echo $SECONDS
