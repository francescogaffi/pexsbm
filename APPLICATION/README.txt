This repository contains codes and materials to reproduce the analysis of the "Infinito network", presented in Sections 5 of the article "Partially Exchangeable Stochastic Block Models for (Node-Colored) Multilayer Networks". Results in Section S4.3 of the Supplementary Material are also reproduced.

The analyses are performed with a MacBook Air (M1, 2020), CPU 8–core and 8GB RAM (macOS Monterey, version 12.5), using the R version 4.2.2

--> REPRODUCIBILITY NOTE: A seed is set at the beginning of each code to ensure full reproducibility. Slight changes in the final numbers reported in Tables and Figures (if any) depend on which version of the external R libraries employed has been used in the implementation of the codes. This is due to internal changes of certain functions when versions of some packages have been updated. However, the magnitude of these minor variations (if any) is negligible and does not affect the final conclusions.

------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------
**IMPORTANT: TO RUN SUCH AN ANALYSIS, THE R WORKING DIRECTORY MUST BE SET IN THIS FOLDER**
------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------

The repository contains the following files.

##########################
### SOURCE CODES #########
##########################

• "pex-sbm.R". It contains all the source R functions which are required to perform posterior computation, inference and prediction under the pEx-SBM class (with a particular focus on the three examples studied in the article: H-DP, H-NSP and H-DP with hyperprior).

• "numcl.R". It contains useful R functions to study the most remarkable priors underlying the proposed pEx-SBM class (H-DP, H-NSP and H-DP with hyperprior), with a particular focus on the expected prior number of non-empty clusters.

##########################
### DATA #################
##########################

• "crime_net.RData". It contains the criminal network data analyzed in the article (adjacency matrix, node attributes, and additional useful quantities for graphical representation of the network). The original data are available at "https://sites.google.com/site/ucinetsoftware/datasets/covert-networks/ndrangheta-mafia-2". See also the GitHub repository "https://github.com/danieledurante/ESBM" for information on data pre-processing.

##########################
### STEP-BY-STEP CODE ####
##########################

• "application_CRIME.R". It contains the step-by-step code to reproduce the results in Sections 5 of the main article, and Section S4.3 of the Supplementary Material.


-------------------------------------------------------------
---- INSTRUCTIONS TO REPRODUCE THE RESULTS ------------------
-------------------------------------------------------------

--> 1. Set the R working directory in this folder
--> 2. Open the file "application_CRIME.R" and run the code

-------------------------------------------------------------
-------------------------------------------------------------
-------------------------------------------------------------
