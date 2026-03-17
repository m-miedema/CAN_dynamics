#!/bin/bash

#SBATCH --array=1-25
#SBATCH --time=MYHR:00:00   # walltime 
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks) 
#SBATCH --nodes=1   # number of nodes 
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=64G  # memory per CPU core 
#SBATCH --account=def-gmitsis  
#SBATCH --job-name=t_ants_MYSES-MYTASK-MYPIP
#SBATCH --mail-user=mary.miedema@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --output=out_ants_t_MYSES-MYTASK-MYPIP_%j.txt

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

# first transform using the affine transformation to structural space from FSL FEAT
bids_root=sub-${sub_num}_ses-${ses_num}_task-${task}

in_func=./central/func_${pipel}/${bids_root}_func.nii.gz
in_reg_root=/home/miedemam/scratch/MGH_derivatives/feat/sub-${sub_num}/ses-${ses_num}/${task}.feat/reg
out_struct=${out_dir}/${bids_root}_stuct.nii.gz
if [ "$ses_num" = "02" ]; then
func2highres_premat=${in_reg_root}/example_func2highres.mat
fi
if [ "$ses_num" = "03" ]; then
func2highres_premat=${in_reg_root}/example_func2highres2.mat
fi

flirt -in $in_func -ref ./anat/sub-${sub_num}_ses-02_acq-mprage_T1w_brain.nii.gz -out $out_struct -init $func2highres_premat -applyxfm
echo "Finished transforming the functional scan into structural space."

# now apply the nonlinear ANTS transformation into MNI space
module load ants
MNItemplate=/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/CUDA/gcc9/cuda11.0/fsl/6.0.4/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz
transform_dir=/home/miedemam/scratch/MGH_derivatives/transforms/sub-${sub_num}/
T1toMNItransform_1=${transform_dir}/sub-${sub_num}_ses-02_acq-mprage_T1w_2_MNI152_T1_1mm0GenericAffine.mat
T1toMNItransform_2=${transform_dir}/sub-${sub_num}_ses-02_acq-mprage_T1w_2_MNI152_T1_1mm1Warp.nii.gz

antsApplyTransforms -v -d 3 -e 3 -i $out_struct -r $MNItemplate -o ${out_dir}/${bids_root}_MNI_brainstem.nii.gz -n BSpline -t $T1toMNItransform_2 -t $T1toMNItransform_1
echo "Finished transforming the functional scan into MNI space."

echo "Completely finished! This took:"
echo $SECONDS
