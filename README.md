# pexsbm: Partially exchangeable stochastic block models

This repository is associated with the article [**Partially Exchangeable Stochastic Block Models for (Node-Colored) Multilayer Networks**](https://arxiv.org/abs/2410.10619), and aims at providing detailed materials and codes to implement the general **pEx-SBM** class presented in the article and to reproduce the results presented in Sections 4 and 5, and in the Supplementary Material.

The documentation is organized in two main folders (`SIMULATIONS` and `APPLICATION`, with multiple subfolders) and one source file (`pex-sbm.R`) described below (a `packages.txt` file with further details on the `R` packages employed is also provided).

- `pex-sbm.R`. It contains all the **source** `R` **functions** which are required to perform posterior computation, inference and prediction under the **pEx-SBM** class (with a particular focus on the three examples studied in the article: H-DP, H-NSP and H-DP with hyperprior).

- `SIMULATIONS`. It contains three subfolders to reproduce the results for the **simulation studies** in Section 4:
	- `CONSISTENCY`: Contains codes/materials to reproduce the consistency results in Figure 3 of Section 4.
	- `REPLICATES`: Contains codes/materials to reproduce Table 1, Figure 2 and additional results in Section 4.
	- `ROBUSTNESS`: Contains codes/materials to reproduce the resuts in Section S4.2 (Supplementary Material).
	
- `APPLICATION`. It contains codes and materials to reproduce the **analysis of the *Infinito* network**, presented in Sections 5 of the article. The folder contains also the network studied in the article (see `crime_net.RData`).

- `packages.txt`. It contains detailed information on the version number of the `R` packages employed.

Each folder and subfolder in `SIMULATIONS` and `APPLICATION` contains additional `README.txt` files with details for the implementation.

The analyses are performed with a **MacBook Air (M1, 2020), CPU 8–core and 8GB RAM (macOS Monterey, version 12.5)**, using the `R` version 4.2.2

All the above functions rely on a basic and reproducible `R` implementation, mostly meant to provide a clear understanding of the computational routines associated with the proposed model. Optimized computational routines relying on C++ coding can be easily considered. 

1. **Note on runtime**: Each Gibbs-sampler requires few minutes (< 2 minutes) for a single study on a standard laptop. In the extensive studies within the folders `SIMULATIONS` and `APPLICATION`, such a Gibbs-sampler is run for several replicated studies, under different methods, and in varying scenarios. This leads to a runtime in the order of hours to reproduce in full generality all the results in the article. 

2. **Note on reproducibility**: A seed is set at the beginning of each code to ensure full reproducibility. Slight changes in the final numbers reported in Tables and Figures (if any) depend on which version of the external `R` libraries employed has been used in the implementation of the codes. This is due to internal changes of certain functions when versions of some packages have been updated. However, the magnitude of these minor variations (if any) is negligible and does not affect the final conclusions.

3. **Note on figures**: Step-by-step codes and guidelines to reproduce the figures in the article are provided in the corresponding folders. See, in particular, `application_CRIME.R` (folder `APPLICATION`), `posterior_analyses_ALL.R` (folder `SIMULATIONS/CONSISTENCY`), `posterior_analyses_ALL.R` (folder `SIMULATIONS/REPLICATES`) and `traceplots.R` (folder `SIMULATIONS/REPLICATES/TRACEPLOTS`).
