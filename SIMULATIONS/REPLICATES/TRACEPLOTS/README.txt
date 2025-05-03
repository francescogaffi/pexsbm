This repository contains codes and materials to reproduce the traceplots in Figure S.1 of the Supplementary Material of the article "Partially Exchangeable Stochastic Block Models for (Node-Colored) Multilayer Networks".

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

#####################################
### DATA ############################
#####################################

• "simulations.RData". It contains the simulated data analyzed in the article (layer division, true underlying partition, block probability matrices in the two scenarios analyzed, simulated adjacency matrices in ten replicated studies under each of the two scenarios). 

######################################
### STEP-BY-STEP CODE AND OUTPUTS ####
######################################

• "traceplots.R". It contains the step-by-step code to reproduce the traceplots in Figure S.1 for H-DP, H-NSP and H-DP with hyperprior in the two scenarios analyzed in Section 4, with a focus on the first of the ten replicated studies.
  --- • This code produces the output "output_TRACEPLOTS.RData" with the log of the beta-binomial likelihood that we monitor in the traceplots 


-------------------------------------------------------------
---- INSTRUCTIONS TO REPRODUCE THE RESULTS ------------------
-------------------------------------------------------------

--> 1. Set the R working directory in this folder
--> 2. Open the file "traceplots.R" and run it to produce the output "output_TRACEPLOTS.RData" and Figure S.1 of the Supplementary Material

-------------------------------------------------------------
-------------------------------------------------------------
-------------------------------------------------------------
