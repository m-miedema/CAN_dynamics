#!/bin/bash

#SBATCH --array=1-4
#SBATCH --time=13:00:00   # walltime 
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks) 
#SBATCH --nodes=1   # number of nodes 
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=40G  # memory per CPU core 
#SBATCH --account=def-gmitsis  
#SBATCH --job-name=BN_T1toMNI_reg
#SBATCH --output=out_BN_reg_%j.txt
#SBATCH --mail-user=mary.miedema@mail.mcgill.ca
#SBATCH --mail-type=ALL

module load StdEnv/2023  
module load fsl/6.0.7.7

# get the subject number
cd /home/miedemam/scratch/MGH_derivatives
sub_num=$(sed -n "${SLURM_ARRAY_TASK_ID}p" run_subs_4.txt) 

echo "Processing subject:"
echo $sub_num

# copy over the structural file
anat_folder=/home/miedemam/scratch/bids_MGH_FINAL/sub-${sub_num}/ses-02/anat/
out_folder=/home/miedemam/scratch/MGH_derivatives/transforms/sub-${sub_num}/
mkdir $out_folder
cp ${anat_folder}sub-${sub_num}_ses-02_acq-mprage_T1w.nii.gz ${out_folder}sub-${sub_num}_ses-02_acq-mprage_T1w.nii.gz

# call the Brainstem Navigator script
## cd scripts
## movingimage=sub-${sub_num}_ses-02_acq-mprage_T1w.nii.gz
## fixedimage=/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/CUDA/gcc9/cuda11.0/fsl/6.0.4/fsl/data/standard/MNI152_T1_1mm
## dirIN=$out_folder

module load ants

## ./antsreg_CC_deform_T1wstructural2MNI.sh $movingimage $fixedimage $dirIN

# cannibalize the Brainstem Navigator script -------------

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=8

moving=${out_folder}sub-${sub_num}_ses-02_acq-mprage_T1w.nii.gz
fixed=/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/CUDA/gcc9/cuda11.0/fsl/6.0.4/fsl/data/standard/MNI152_T1_1mm.nii.gz
dirOUT=${out_folder}

out=`imglob ${moving} | xargs basename`_2_`imglob ${fixed} | xargs basename`
#defout=`imglob ${moving} | xargs basename`_Deformed_2_`imglob ${fixed} | xargs basename`.nii.gz
defout=$dirOUT/${out}

TranslationParameters=" -m CC[ ${fixed}, ${moving}, 1, 4 ] -u 1 -t Translation[1] -f 6x4x2x1 -s 4x2x1x0                   -c [10000x10000x0x0, 1.e-8, 10]  -r [ ${fixed}, ${moving}, 1]"

RigidParameters="       -m CC[ ${fixed}, ${moving}, 1, 4 ] -u 1 -t Rigid[1]       -f 6x4x2x1 -s 4x2x1x0 -w [0.005, 0.995] -c [10000x10000x0x0, 1.e-8, 10]  "

AffineParameters="      -m CC[ ${fixed}, ${moving}, 1, 4 ] -u 1 -t Affine[1]      -f 6x4x2x1 -s 4x2x1x0 -w [0.005, 0.995] -c [10000x10000x1500x20, 1.e-8, 10]  "

DeformParameters="      -m CC[ ${fixed}, ${moving}, 1, 4 ] -u 1 -t SyN[0.2,3,0]   -f 6x4x2x1 -s 3x2x1x0 -w [0.005, 0.995] -c [200x200x200x50, 1e-8, 10] "


# -r [ ${fixed}, ${moving}, 1] 1 = center of mass of alignment; 0 alignment of geometric center of the images; 2 alignment of the origin of the images
# -m mattes[ ${fixed}, ${moving}, 1, 32, Regular, 0.05 ]

echo "Translation + Rigid + Affine + Deformable"

antsRegistration -d 3 -o ${defout} ${TranslationParameters} ${RigidParameters} ${AffineParameters} ${DeformParameters} -v 0

antsApplyTransforms -d 3 -e 0 -i ${moving} -r ${fixed} -o ${defout}.nii.gz -t ${defout}1Warp.nii.gz -t ${defout}0GenericAffine.mat -v 0



echo "Finished structural registration!"
echo "This job ran for:"
echo $SECONDS
