This repository contains codes and materials to reproduce Table 1, Figure 2 and further results in Section 4 of the article "Partially Exchangeable Stochastic Block Models for (Node-Colored) Multilayer Networks". 

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

• "JCDC.cpp". It contains the source R functions developed by Yuan Zhang to implement the community detection method with node feature (JCDC) proposed by Zhang et al., (2016).

#####################################
### DATA ############################
#####################################

• "simulations.RData". It contains the simulated data analyzed in the article (layer division, true underlying partition, block probability matrices in the two scenarios analyzed, simulated adjacency matrices in ten replicated studies under each of the two scenarios). 

######################################
### STEP-BY-STEP CODE AND OUTPUTS ####
######################################

• "simulations_HDP.R". It contains the step-by-step code to reproduce the results of the analysis under pEx-SBM with H-DP
  --- • This code produces the output "output_HDP.RData" with performance measures for H-DP

• "simulations_HNSP.R". It contains the step-by-step code to reproduce the results of the analysis under pEx-SBM with H-NSP
  --- • This code produces the output "output_HNSP.RData" with performance measures for H-NSP

• "simulations_HDP_hyperprior.R". It contains the step-by-step code to reproduce the results of the analysis under pEx-SBM with H-DP hyperprior
  --- • This code produces the output "output_HDP_hyperprior.RData" with performance measures for H-DP with hyperprior 

• "simulations_ESBM.R". It contains the step-by-step code to reproduce the results of the analysis under ESBM (Legramanati et al. 2022)
  --- • This code produces the output "output_ESBM.RData" with performance measures for ESBM (Legramanati et al. 2022)

• "simulations_COMPETITORS.R". It contains the step-by-step code to reproduce the results of the analysis under Louvain, SBM (greed), JCDC and CASC
  --- • This code produces the output "output_COMPETITORS.RData" with performance measures for Louvain, SBM (greed), JCDC and CASC

• "posterior_analyses_ALL.R". It contains the step-by-step code to produce the results in Table 1, Figure 2 and additional analyses reported in Section 4 from the analysis of the output files presented above.

######################################
### OTHER FOLDERS ####################
######################################

• "TRACEPLOTS". It contains code and materials to reproduce the traceplots in Figure S.1 of the Supplementary Material


-------------------------------------------------------------
---- INSTRUCTIONS TO REPRODUCE THE RESULTS ------------------
-------------------------------------------------------------

--> 1. Set the R working directory in this folder
--> 2. Run, separately, the codes "simulations_HDP.R", "simulations_HNSP.R", "simulations_HDP_hyperprior.R", "simulations_ESBM.R" and "simulations_COMPETITORS.R"
--> 3. Step 2. produces the outputs "output_HDP.RData", "output_HNSP.RData", "output_HDP_hyperprior.RData", "output_ESBM.RData" and "output_COMPETITORS.RData"
--> 4. Open the file "posterior_analyses_ALL.R" and run it to produce the results in Table 1, Figure 2 and additional analyses reported in Section 4

NOTE: each code in step 2. performs the analysis on ten replicated studies in two different scenarios. Hence, although each Gibbs-sampler requires few minutes on a single study, re-running such a sampler over multiple replicated studies increases the time.

-------------------------------------------------------------
-------------------------------------------------------------
-------------------------------------------------------------
