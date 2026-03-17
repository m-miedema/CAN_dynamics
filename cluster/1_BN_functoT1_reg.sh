#!/bin/bash

#SBATCH --array=1-25
#SBATCH --time=1:00:00   # walltime 
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks) 
#SBATCH --nodes=1   # number of nodes 
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G  # memory per CPU core 
#SBATCH --account=def-gmitsis  
#SBATCH --job-name=BN_functoT1_reg
#SBATCH --output=out_BN_reg_%j.txt
#SBATCH --mail-user=mary.miedema@mail.mcgill.ca
#SBATCH --mail-type=ALL

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

# set the script parameters
out_folder=/home/miedemam/scratch/MGH_derivatives/transforms/sub-${sub_num}/
bids_root=sub-${sub_num}_ses-${ses_num}_task-${task}
mkdir ${out_folder}/ses-$ses_num
mkdir ${out_folder}/ses-${ses_num}/$task
mkdir ${out_folder}/ses-${ses_num}/${task}/${pipel}
cp ./anat/sub-${sub_num}_ses-02_acq-mprage_T1w_brain.nii.gz ${out_folder}/ses-${ses_num}/${task}/${pipel}/sub-${sub_num}_ses-02_acq-mprage_T1w_brain.nii.gz

moving=/home/miedemam/scratch/MGH_derivatives/feat/sub-${sub_num}/ses-${ses_num}/${task}.feat/mean_func.nii.gz
fixed_brain=${out_folder}/ses-${ses_num}/${task}/${pipel}/sub-${sub_num}_ses-02_acq-mprage_T1w_brain.nii.gz
dirOUT=$out_folder/ses-${ses_num}/${task}/${pipel}
full_fMRI=/home/miedemam/scratch/MGH_derivatives/central/func_${pipel}/${bids_root}_func.nii.gz

# CANNIBALIZE THE BRAINSTEM NAVIGATOR SCRIPT
cd $dirOUT
## fixedimage=/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/CUDA/gcc9/cuda11.0/fsl/6.0.4/fsl/data/standard/MNI152_T1_1mm

module load ants
module load afni

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=8

betted_anatomy=$fixed_brain
dirout=$dirOUT

moving_nopath=$(basename ${moving})
if [ $moving_nopath = $moving ]; then moving=$(pwd)/$moving;fi
moving_noext=$(basename $moving_nopath .nii)
moving_noext=$(basename $moving_noext .nii.gz)

full_fMRI_nopath=$(basename ${full_fMRI})
if [ $full_fMRI_nopath = $full_fMRI ]; then full_fMRI=$(pwd)/$full_fMRI;fi
full_fMRI_noext=$(basename $full_fMRI_nopath .nii)
full_fMRI_noext=$(basename $full_fMRI_noext .nii.gz)

betted_anatomy_nopath=$(basename ${betted_anatomy})
if [ $betted_anatomy_nopath = $betted_anatomy ]; then betted_anatomy=$(pwd)/$betted_anatomy;fi

mkdir -p $dirout


if [ ! -f ${dirout}/${betted_anatomy_nopath} ]; then
 cp ${betted_anatomy} ${betted_anatomy_nopath}
 betted_anatomy=$(pwd)/${betted_anatomy_nopath}
fi
3drefit -view tlrc -space tlrc ${betted_anatomy}

cd $dirout
echo "############### STARTING AFNI REGISTRATION, SLOWER OPTION ###############"
if [ ! -f ${dirout}/${moving_noext}_AFNIreg2T1w_linear.nii.gz ] | [ ! -f ${dirout}/${moving_noext}epi2anat_mat.aff12.1D ]; then
echo "compute affine coregistration to MPRAGE..."
align_epi_anat.py -volreg off -anat ${betted_anatomy} -master_epi BASE -suffix epi2anat -epi ${moving} -epi_strip None -epi_base 0 -epi2anat -anat_has_skull no -giant_move
3dAFNItoNIFTI -prefix ${dirout}/${moving_noext}_AFNIreg2T1w_linear.nii.gz ${dirout}/${moving_noext}epi2anat+tlrc
fi
echo "compute additional nonlinear coregistration to MPRAGE..."
3dQwarp -source ${dirout}/${moving_noext}epi2anat+tlrc -base ${betted_anatomy} -prefix ${dirout}/${moving_noext}epi2anatQ -lpc -maxlev 4 -blur 0 1 -resample -Qfinal -inedge
echo "create control volume to check linear coregistration quality..."
rm ${dirout}/${moving_noext}epi2anat+tlrc*
3dAFNItoNIFTI ${dirout}/${moving_noext}epi2anatQ_WARP+tlrc
rm ${dirout}/${moving_noext}epi2anatQ_WARP+tlrc*
rm ${dirout}/${moving_noext}epi2anatQ_Allin*
echo "create control volume to check nonlinear coregistration quality..."
3dNwarpApply -source ${moving} -nwarp "${dirout}/${moving_noext}epi2anatQ_WARP.nii ${dirout}/${moving_noext}epi2anat_mat.aff12.1D" -prefix ${dirout}/${moving_noext}_AFNIreg2T1w_slower_nonlinear_QC.nii.gz -master ${betted_anatomy}
echo "DONE COMPUTING TRANSFORMS"
echo "FOR QUALITY CHECK: afni ${betted_anatomy} ${dirout}/${moving_noext}_AFNIreg2T1w_linear.nii.gz ${dirout}/${moving_noext}_AFNIreg2T1w_slower_nonlinear_QC.nii.gz"

echo "applying transforms to the full functional dataset..."
3dNwarpApply -source ${full_fMRI} -nwarp "${dirout}/${moving_noext}epi2anatQ_WARP.nii ${dirout}/${moving_noext}epi2anat_mat.aff12.1D" -prefix ${dirout}/${full_fMRI_noext}_AFNIreg2T1w_slower_nonlinear.nii.gz -master ${betted_anatomy}
echo "DONE"

echo "Finished structural registration!"
echo "This job ran for:"
echo $SECONDS
