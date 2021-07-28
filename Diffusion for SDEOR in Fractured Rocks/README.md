#This file provides instructions to run the new diffusion NNC codes in HPC, or a machine with 'Shale' module installed.
#Search for 'Hassan' inside each .m file to see changes.

1- Add "MSdiffNatVars.m" file to Utils folder in Shale module for static implementation.

2- Update "EDFMshalegrid.m" file inside Shale/hfmFns.

3- Update "assembleShaleGlobalGrid.m" file inside Shale/hfmFns.

4- Update "fracturematrixShaleNNC3D.m" file inside Shale/hfmFns.

5- Update "fgridNNCs3D.m" file inside hfm/edfm-hw/Fracture-Fracture Intersection Preprocessing.

6- Update "pMatFracNNCs3D.m" file inside Shale/pEDFM_nnc.

7- Update "setupPEDFMOpsTPFA.m" file inside Shale/pEDFM_nnc.

8- Update "setupEDFMOperatorsTPFA.m" file inside hfm/edfm-hw/Operator Setup.

9- Add "getFaceTransmissibilityShale_Diffusion.m" file to Utils folder in Shale module.

10- Add "computeTransShale_Diffusion.m" file to Utils folder in Shale module.

11- Update "eqnsShaleNaturalVars.m" file inside Shale/compositionalFns.

12- Update "setupOperatorsTPFA.m" file inside ad-core/utils.

13- Add "explicitStencil.m" file to Utils folder in Shale module. #Domain meshing for HF explicit modeling

14- Confirm "GenericNaturalVariablesModel.m" file is added to compositional/models. #To use state functions for fugacity derivative computations

