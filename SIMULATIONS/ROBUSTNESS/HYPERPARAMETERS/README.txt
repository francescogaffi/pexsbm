This repository contains codes and materials to reproduce the analyses on the robustness of pEx-SBM in the article "Partially Exchangeable Stochastic Block Models for (Node-Colored) Multilayer Networks" with respect to the choice of the hyperparameters. In particular, this repository reproduces results in Table S.1 of the Supplementary Material (column "hyperparameters"). 

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

• "numcl.R". It contains useful R functions to study the most remarkable priors (H-DP, H-NSP and H-DP with hyperprior) underlying the proposed pEx-SBM class, with a particular focus on the expected prior number of non-empty clusters.

#####################################
### DATA ############################
#####################################

• "simulations.RData". It contains the simulated data analyzed in the article (layer division, true underlying partition, block probability matrices in the two scenarios analyzed, simulated adjacency matrices in ten replicated studies under each of the two scenarios). 

######################################
### STEP-BY-STEP CODE AND OUTPUTS ####
######################################

• "simulations_HDP.R". It contains the step-by-step code to reproduce the results on the robustness analysis under pEx-SBM with H-DP
  --- • This code produces the output "output_HDP.RData" with the performance measures for H-DP

• "simulations_HNSP.R". It contains the step-by-step code to reproduce the results on the robustness analysis under pEx-SBM with H-NSP
  --- • This code produces the output "output_HNSP.RData" with the performance measures for H-NSP

• "simulations_HDP_hyperprior.R". It contains the step-by-step code to reproduce the results on the robustness analysis under pEx-SBM with H-DP hyperprior
  --- • This code produces the output "output_HDP_hyperprior.RData" with the performance measures for for H-DP with hyperprior 

• "simulations_ESBM.R". It contains the step-by-step code to reproduce the results of the robustness analysis under ESBM by Legramanti et al. (2022)
  --- • This code produces the output "output_ESBM.RData" with the performance measures for ESBM

• "posterior_analyses_ALL.R". It contains the step-by-step code to produce column "hyperparameters" in Table S.1 of the Supplementary Material from the analysis of the output files presented above.

NOTE: Results for ESBM are provided for completeness, but are not reported in the article, since the focus is to assess robustness of the proposed pEx-SBM.

-------------------------------------------------------------
---- INSTRUCTIONS TO REPRODUCE THE RESULTS ------------------
-------------------------------------------------------------

--> 1. Set the R working directory in this folder
--> 2. Run, separately, the codes "simulations_HDP.R", "simulations_HNSP.R", "simulations_HDP_hyperprior.R" and "simulations_ESBM.R"
--> 3. Step 2. produces the outputs "output_HDP.RData", "output_HNSP.RData", "output_HDP_hyperprior.RData" and "output_ESBM.RData"
--> 4. Open the file "posterior_analyses_ALL.R" and run it to obtain column "hyperparameters" in Table S.1 of the Supplementary Material

NOTE 1: each code in step 2. performs the analysis on ten replicated studies in two different scenarios. Hence, although each Gibbs-sampler requires few minutes on a single study, re-running such a sampler over multiple replicated studies increases the time.

NOTE 2: Results for ESBM are provided for completeness, but are not reported in the article, since the focus is to assess robustness of the proposed pEx-SBM.
-------------------------------------------------------------
-------------------------------------------------------------
-------------------------------------------------------------
