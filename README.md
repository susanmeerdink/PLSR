# PLSR
Project related to Partial Least Squares Regression for prediction of foliar traits

### INPUTS
---
  * This code works with a spectral database which was collected in the laboratory with an Analytical Spectra Device (ASD) and Nicolet Model 4700 Interferometer Spectrometer, which measures the Visible Short Wave/InfraRed and Thermal InfraRed spectral range respectively. Additionally this code takes a database of plant species leaf traits included cellulose (%), lignin (%), leaf mass per area (g/m^2), nitrogen (%), and water content (%). 
  * For more information about methods and dataset see: 
	* `Meerdink, S.K., Roberts, D.A., King, J.Y., Roth, K.L., Dennison, P., Amaral, C.H., & Hook, S.J. (in revision). Linking seasonal foliar traits to VSWIR-TIR spectroscopy across California ecosystems. Remote Sensing Environment.

### DESCRIPTION OF CODE
---
  * PLSR_HyspIRI_Project_Single_Trait.m: This script selects data for a specific leaf trait by looping through the 6 spectrums (Full, VSWIR, TIR, HyspIRI, AVIRIS, and HyTES) and 5 model categories( Broadleaf, Needleleaf, General, Spring, Summer, and Fall). It calls functions to develop training/validation datasets, run Partial Least Squares Regression (PLSR), and display the results. It outputs the results to a csv and saves the workspace.
  * splitValCalBefore.m: Function that splits the dataset into calibration (model development) and independent validation (used for model assessment).
  * determinefactors.m: Function that determines the number of factors or components that will be used to run PLSR
  * plsr_with_bins_prop.m: Function that predicts leaf trait using PLSR
  * plsfigure.m: Function that displays predicted versus observed plots from the PLSR 
  * loop_through_workspaces_save_plsr_results.m: Script that reads through all traits PLSR results and compiles them into one csv
  * loop_through_workspaces_save_calc_BETA.m: Script that reads through all traits PLSR BETA coefficients and compiles them into one csv
  * Figures_for_Manuscript: This script contains figures that are used in the RSE manuscript.
