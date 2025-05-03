This repository contains codes and materials to reproduce the analyses on the robustness of pEx-SBM in the article "Partially Exchangeable Stochastic Block Models for (Node-Colored) Multilayer Networks" to uninformative layers. In particular, this repository reproduce results discussed in Section S4.2 of the Supplementary Material. 

The analyses are performed with a MacBook Air (M1, 2020), CPU 8–core and 8GB RAM (macOS Monterey, version 12.5), using the R version 4.2.2

--> REPRODUCIBILITY NOTE: A seed is set at the beginning of each code to ensure full reproducibility. Slight changes in the final numbers reported in Tables and Figures (if any) depend on which version of the external R libraries employed has been used in the implementation of the codes. This is due to internal changes of certain functions when versions of some packages have been updated. However, the magnitude of these minor variations (if any) is negligible and does not affect the final conclusions.

------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------
**IMPORTANT: TO RUN SUCH AN ANALYSIS, THE R WORKING DIRECTORY MUST BE SET IN THIS FOLDER**
------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------

The repository contains the following files.

#####################################
### SOURCE CODES ####################
#####################################

• "pex-sbm.R". It contains all the source R functions which are required to perform posterior computation, inference and prediction under the pEx-SBM class (with a particular focus on the three examples studied in the article: H-DP, H-NSP and H-DP with hyperprior).

• "casc.R". It contains the source R functions developed in the GitHub repository (https://github.com/norbertbin/rCASC) to implement the covariate assisted spectral clustering method by Binkiewicz et al. (2017). IMPORTANT NOTE: To implement CASC, first install the packages "rCASC" and "SpecClustPack" following the instructions at https://github.com/norbertbin/rCASC and https://github.com/norbertbin/SpecClustPack, respectively.

#####################################
### DATA ############################
#####################################

• "simulations_perm.RData". It contains the simulated data to assess robustness to uninformative layers (layer division, true underlying partition, block probability matrices in the two scenarios analyzed, simulated adjacency matrices in ten replicated studies under each of the two scenarios). The focus will be on random permutations of the layer division in scenario 1.

######################################
### STEP-BY-STEP CODE AND OUTPUTS ####
######################################

• "simulations_ROBUSTNESS.R". It contains the step-by-step code to reproduce the results on the robustness analysis under pEx-SBM with H-DP, H-NSP, H-DP with hyperprior and CASC
  --- • This code produces the output "output_ROBUSTNESS.RData" with the performance measures for H-DP, H-NSP, H-DP with hyperprior and CASC

• "posterior_analyses_ALL.R". It contains the step-by-step code to reproduce the results discussed in Section S4.2 of the Supplementary Material from the analysis of the output file presented above.

-------------------------------------------------------------
---- INSTRUCTIONS TO REPRODUCE THE RESULTS ------------------
-------------------------------------------------------------

--> 1. Set the R working directory in this folder
--> 2. Open the file "simulations_ROBUSTNESS.R" and run it
--> 3. Step 2. produces the output "output_ROBUSTNESS.RData"
--> 4. Open the file "posterior_analyses_ALL.R" and run it to reproduce the results discussed in Section S4.2 of the Supplementary Material

NOTE: The code in step 2. performs the analysis on ten replicated studies under each of the four methods (H-DP, H-NSP, H-DP with hyperprior and CASC). Hence, although each Gibbs-sampler requires few minutes on a single study, re-running such a sampler over multiple replicated studies increases the time. Hence, since output files are already provided, one can also focus only on steps 1. and 4.

-------------------------------------------------------------
-------------------------------------------------------------
-------------------------------------------------------------
