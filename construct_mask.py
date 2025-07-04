"""
Implementation of functions to create a 3D volumetric mask based on a spherical ROIs defined in MNI coordinates.

The coordinates provided in this script correspond to regions of the central autonomic network as identified by Beissner et al., 2013

Created on July 4 2025
@author: Mary Miedema, based on https://neurostars.org/t/create-a-10mm-sphere-roi-mask-around-a-given-coordinate/28853/3 
"""

from nilearn import datasets, plotting
from nilearn.masking import _unmask_3d
from nilearn.maskers import nifti_spheres_masker
import nibabel as nib
from nibabel import Nifti1Image

# the brain space can be adapted here as needed
brain_mask = datasets.load_mni152_brain_mask()

# define a triplet for each ROI
coords = [(0, 10, 40), (48, -26, 46), (-20, -8, -12), (-2, 38, -18), (44, 18, -6), (-4, -16, 6), (-44, -36, 42), (-32, -20, 14), (-46, -66, -28), (20, 36, 34), (30, -22, -16), (-20, -6, -18), (-40, 0, 12), (-6, -44, 34), (-56, 6, 8), (50, -24, 2), (44, -38, 14), (-10, -62, -20), (40, 2, 12)]

_, A = nifti_spheres_masker._apply_mask_and_get_affinity(
    seeds=coords,
    niimg=None,
    radius=5, # choose ROI radius
    allow_overlap=False, 
    mask_img=brain_mask)

from nilearn.masking import _unmask_4d
sphere_mask = _unmask_4d(X=A.toarray(),mask=brain_mask.get_fdata().astype(bool))
sphere_mask = Nifti1Image(sphere_mask, brain_mask.affine)

# convert masks from 3-D to 4-D, with separate values for each ROI
data = sphere_mask.get_fdata()
import numpy as np
summed_data = np.sum(mult_data, axis=-1)
mult_data = data;
for i in range(1, sphere_mask.shape[3]):
    print (i)
    mult_data[:,:,:,i]=mult_data[:,:,:,i]+i*data[:,:,:,i]
new_img = nl.image.new_img_like(sphere_mask, summed_data)
nib.save(new_img,"MASKNAME.nii.gz")

####################################################################################################################################
